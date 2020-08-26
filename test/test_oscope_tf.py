"""
TensorFlow 2 OscoNet code
"""
import numpy as np
import pytest
import tensorflow as tf

from OscopeBootstrap import qvalue
from OscopeBootstrap.create_edge_network_represention import create_edge_network_representation
from OscopeBootstrap.oscope_tf import PRECISION_fp, calc_e2, calc_e2_many_genes, find_best_psi_for_each_gene_pair, \
    PRECISION_int, get_permuted_cost, get_pvalues, flatten_upper_triangular, get_symmetric_matrix_from_upper_triangular
from OscopeBootstrap.SyntheticDataset import GetSimISyntheticData, true_adj_matrix
from OscopeBootstrap.oscope_tf import bootstrap_hypothesis_test, get_accuracy, get_metrics_for_different_qvalue_thresholds


def calc_e2_np(X, Y, psi):
    return np.sum(np.square(np.square(X) + np.square(Y) - 2 * X * Y * np.cos(psi) - np.square(np.sin(psi))))


def calc_e2_many_genes_np(X_many_genes: np.ndarray, psi_ng: np.ndarray):
    '''

    :param X_many_genes: G X N tensor of gene expression
    :param psi_ng: G X G tensor of phase shift - should be symmetric
    :return: total cost across all genes
    '''
    G = X_many_genes.shape[0]
    c = 0
    for ix in range(G):
        for iy in range(G):
            c += calc_e2_np(X_many_genes[ix, :], X_many_genes[iy, :], psi_ng[ix, iy])
    return c


def create_single_group_example(N, std_noise, phase_shift):
    t = np.linspace(0, 2 * np.pi, N)
    G = 4
    data = np.zeros((G, N))
    data[0, :] = np.sin(t) + std_noise * np.random.randn(N)
    data[1, :] = np.sin(t + phase_shift) + std_noise * np.random.randn(N)
    data[2, :] = std_noise * np.random.randn(N)
    data[3, :] = std_noise * np.random.randn(N)
    return data


def test_get_symmetric_matrix_from_upper_triangular():
    flatten_vector = np.array([1, 2, 3, 4, 5, 6])
    a = get_symmetric_matrix_from_upper_triangular(4, flatten_vector)
    np.testing.assert_equal(a, a.T)


def test_calc_e2():
    np.random.seed(42)
    N = 10
    X = tf.constant(np.random.randn(N,), dtype=PRECISION_fp)
    Y = tf.constant(np.random.randn(N, ), dtype=PRECISION_fp)
    psi = tf.constant(np.array(3.), dtype=PRECISION_fp)

    assert calc_e2(X, X, tf.constant(np.array(0.), dtype=PRECISION_fp)) == 0, 'must get minimum cost for identical gene with 0 phase'

    e_tf = calc_e2(X, Y, psi)
    e_np = calc_e2_np(X.numpy(), Y.numpy(), psi.numpy())
    np.testing.assert_almost_equal(e_tf, e_np, decimal=1)


def test_calc_e2_many_genes():
    G = 5
    N = 10
    X_many_genes = tf.constant(np.random.randn(G, N), dtype=PRECISION_fp)
    psi_ng = tf.constant(np.random.randn(G, G), dtype=PRECISION_fp)
    # make sure we include 0 as possible phase
    cost = calc_e2_many_genes(X_many_genes, psi_ng)
    cost_np = calc_e2_many_genes_np(X_many_genes.numpy(), psi_ng.numpy())
    # np.testing.assert_almost_equal(cost, cost_np) Big differences!
    assert np.all(cost > 0)


def test_find_best_psi_for_each_gene_pair():
    np.random.seed(42)
    tf.random.set_seed(42)
    # construct example
    phase_shift = np.pi
    N = 10
    G = 4
    data_np = create_single_group_example(N, 0.1, phase_shift=phase_shift)
    data = tf.constant(data_np, dtype=PRECISION_fp)
    # candidate_psi = tf.linspace(0, 2 * tf.constant(np.pi), dtype=PRECISION)
    candidate_psi = tf.constant(np.array([phase_shift, phase_shift/2]), dtype=PRECISION_fp)

    n_permutations = tf.constant(np.array(20), dtype=PRECISION_int)
    psi_ng = tf.Variable(tf.zeros((G, G), dtype=PRECISION_fp) * tf.constant(np.inf, dtype=PRECISION_fp))
    cost_ng = tf.Variable(tf.ones((G, G), dtype=PRECISION_fp) * tf.constant(np.inf, dtype=PRECISION_fp))
    cost_permuted = tf.Variable(tf.ones((G, G, n_permutations), dtype=PRECISION_fp) * tf.constant(np.inf, dtype=PRECISION_fp))
    pvalues = tf.Variable(tf.ones((G, G), dtype=PRECISION_fp) * tf.constant(np.inf, dtype=PRECISION_fp))

    find_best_psi_for_each_gene_pair(psi_ng, cost_ng, data, candidate_psi=candidate_psi)
    assert psi_ng[0, 1] == phase_shift, 'why picked the other phase shift?'

    get_permuted_cost(cost_permuted, data, candidate_psi, n_permutations)

    get_pvalues(pvalues, cost_ng, cost_permuted)

    # then q-values
    # then check we find the right pair
    pvalue_flatten = flatten_upper_triangular(pvalues.numpy())
    qvalues_flatten, pi0 = qvalue.estimate(pvalue_flatten, verbose=True)
    qvalues = get_symmetric_matrix_from_upper_triangular(pvalues.shape[0], qvalues_flatten)

    adjacency_matrix = qvalues < 0.01
    assert adjacency_matrix[0, 1]
    assert adjacency_matrix[1, 0]
    assert adjacency_matrix.sum() == 2, 'Only one significant pair should exist (0, 1)'

    gene_names = [f'gene{i}' for i in range(4)]
    a = create_edge_network_representation(adjacency_matrix, 1/cost_ng.numpy(), gene_names)
    assert a.shape[1] == 3, 'must have gene1, gene2, weight columns'
    assert a.shape[0] == 1, 'only one gene pair is significant'


@pytest.mark.slow
def test_bootstrap():
    # This is a slow test (>10 secs) so need to run with `pytest --runslow -rs`
    np.random.seed(42)
    tf.random.set_seed(42)
    NG = 5
    G = 20
    N = 100
    ngroups = 1
    alpha = 0.01  # significance level for test

    data_df, phaseG, angularSpeed = GetSimISyntheticData(NG=NG, G=G, ngroups=ngroups,
                                                         N=N, noiseLevel=0)

    adjacency_matrix, qvalues, cost = bootstrap_hypothesis_test(n_bootstrap=30,
                                                                data=data_df.values,
                                                                alpha=alpha,
                                                                grid_points_in_search=30)

    assert qvalues.shape == (G, G)
    assert adjacency_matrix.shape == (G, G)
    assert np.all(~np.isnan(qvalues))
    assert np.all(~np.isnan(adjacency_matrix))
    assert cost.shape == (G, G)

    adjacency_matrix_true = true_adj_matrix(G, angularSpeed)
    correct_ratio = get_accuracy(adjacency_matrix, adjacency_matrix_true)
    assert correct_ratio > .98

    TPR, FDR, FPR = get_metrics_for_different_qvalue_thresholds(qvalues,
                                                                adjacency_matrix_true,
                                                                np.array([alpha]))
    # To get appropriate values we need to increase number of bootstrap samples
    assert TPR > 0.75
    assert FDR < 0.3
    assert FPR < 0.1
