"""
TensorFlow 2 OscoNet code
"""
import argparse
import time
from typing import Tuple
from numba import njit,prange

import numpy as np
import scipy.stats
import tensorflow as tf
import tensorflow_probability as tfp

from OscopeBootstrap import qvalue
from OscopeBootstrap.SyntheticDataset import GetSimISyntheticData, true_adj_matrix

PRECISION_fp = tf.float32
PRECISION_int = tf.int32

@njit
def numpy_calc_e2(X, Y, psi):
    ''' numpy calculation of EvalE '''
    return np.sum(
            (X*X + Y*Y - 2 * X * Y * np.cos(psi) - np.sin(psi)**2)**2
    )


def calc_e2(X: tf.Tensor, Y: tf.Tensor, psi: tf.Tensor):
    ''' TensorFlow version of EvalE. '''
    return tf.math.reduce_sum(
        tf.square(
            tf.square(X) + tf.square(Y) - 2 * X * Y * tf.math.cos(psi) - tf.square(tf.math.sin(psi))
        )
    )


def calc_e2_many_genes(X_many_genes: tf.Tensor, psi_ng: tf.Tensor):
    '''

    :param X_many_genes: G X N tensor of gene expression
    :param psi_ng: G X G tensor of phase shift - should be symmetric
    :return: total cost across all genes
    '''
    c = 0
    G = tf.shape(X_many_genes)[0]
    for ix in range(G):
        for iy in range(G):
            c += calc_e2(X_many_genes[ix, :], X_many_genes[iy, :], psi_ng[ix, iy])
    return c


@tf.function
def find_best_psi_for_each_gene_pair(psi_ng: tf.Variable, cost_ng: tf.Variable, X_many_genes: tf.Tensor, candidate_psi: tf.Tensor):
    '''
    Find best phi and cost for each gene pair
    :param X_many_genes: G X N tensor of gene expression
    :param candidate_psi - values to consider to phase shift
    :return: phase and cost across for each gene pair - only upper triangular filled in
    '''
    G = tf.shape(X_many_genes)[0]
    gn = tf.size(candidate_psi)
    for ix in tf.range(G):
        for iy in tf.range(ix+1, G):
            print(ix, iy)
            for i_p in tf.range(gn):
                c = calc_e2(X_many_genes[ix, :], X_many_genes[iy, :], candidate_psi[i_p])
                if c < cost_ng[ix, iy]:
                    cost_ng[ix, iy].assign(c)
                    psi_ng[ix, iy].assign(candidate_psi[i_p])
    return psi_ng, cost_ng

@njit(parallel=True)
def numpy_find_best_psi_for_each_gene_pair(psi_ng, cost_ng, X_many_genes, candidate_psi):
    '''
    Find best phi and cost for each gene pair
    :param X_many_genes: G X N tensor of gene expression
    :param candidate_psi - values to consider to phase shift
    :return: phase and cost across for each gene pair - only upper triangular filled in
    '''
    G = (X_many_genes.shape)[0]
    gn = (candidate_psi).size

    for ix in prange(G):
        for iy in range(ix+1, G):
            for i_p in range(gn):
                c = numpy_calc_e2(X_many_genes[ix, :], X_many_genes[iy, :], candidate_psi[i_p])
                if c < cost_ng[ix, iy]:
                    cost_ng[ix, iy]=c
                    psi_ng[ix, iy]=candidate_psi[i_p]
    return psi_ng, cost_ng


@tf.function
def get_permuted_cost(
        cost_permuted: tf.Variable,
        X_many_genes: tf.Tensor,
        candidate_psi: tf.Tensor,
        n_permutations: tf.Tensor) -> None:
    '''

    :param cost_permuted: output - permuted cost
    :param X_many_genes: G X N tensor of gene expression
    :param candidate_psi: G tensor - psi values for each gene
    :param n_permutations: number of bootstrap permutations
    '''
    G = tf.shape(X_many_genes)[0]
    gn = tf.size(candidate_psi)
    for ix in tf.range(G):
        for iy in tf.range(ix+1, G):
            for n in tf.range(n_permutations):
                y_random = tf.random.shuffle(X_many_genes[iy, :])
                for i_p in tf.range(gn):
                    c = calc_e2(X_many_genes[ix, :], y_random, candidate_psi[i_p])
                    if c < cost_permuted[ix, iy, n]:
                        cost_permuted[ix, iy, n].assign(c)

@njit(parallel=True)
def numpy_get_permuted_cost(
        cost_permuted,
        X_many_genes,
        candidate_psi,
        n_permutations) -> None:
    '''
    :param cost_permuted: output - permuted cost
    :param X_many_genes: G X N tensor of gene expression
    :param candidate_psi: G tensor - psi values for each gene
    :param n_permutations: number of bootstrap permutations
    '''

    G = (X_many_genes.shape)[0]
    gn = (candidate_psi).size

    for ix in prange(G):
        #print("Doing {0} of {1} genes".format(ix,G))
        for iy in range(ix+1, G):
            for n in range(n_permutations):
                y_random = np.random.permutation(X_many_genes[iy, :])
                for i_p in range(gn):
                    c = numpy_calc_e2(X_many_genes[ix, :], y_random, candidate_psi[i_p])
                    if c < cost_permuted[ix, iy, n]:
                        cost_permuted[ix, iy, n] = c

@tf.function
def get_pvalues(pvalues: tf.Variable, cost_unpermuted: tf.Tensor, cost_permuted: tf.Tensor):
    '''
    Get p-values from bootstrap
    :param cost_unpermuted: G X G - only upper triangular used
    :param cost_permuted: G X G X n permutations
    :return: p-value
    '''
    # This should relationship of our previous implementation using scipy and tfp
    # empirical_dist = tfp.distributions.Empirical(samples=tf.constant(np.array([1,2,3,4])))
    # assert empirical_dist.cdf(3).numpy() == scipy.stats.percentileofscore([1,2,3,4], 3)/100.
    # for each gene pair, calculate p-value
    G = tf.shape(cost_unpermuted)[0]

    for ix in tf.range(G):
        for iy in tf.range(ix+1, G):
            empirical_dist = tfp.distributions.Empirical(samples=cost_permuted[ix, iy, :])
            pvalues[ix, iy].assign(empirical_dist.cdf(cost_unpermuted[ix, iy]))

def numpy_get_pvalues(cost_unpermuted, cost_permuted):
    '''
    Get p-values from bootstrap
    :param cost_unpermuted: G X G - only upper triangular used
    :param cost_permuted: G X G X n permutations
    :return: p-value
    '''

    G = cost_unpermuted.shape[0]
    pvalues = np.full((G, G), np.inf)

    for ix in range(G):
        for iy in range(ix+1, G):
            pvalues[ix, iy] = scipy.stats.percentileofscore(cost_permuted[ix,iy, :], cost_unpermuted[ix,iy])/100.0

    return pvalues


def flatten_upper_triangular(upper_triangular: np.ndarray):
    G = upper_triangular.shape[0]
    idx_r, idx_c = np.triu_indices(G, 1)
    a = np.array([upper_triangular[i,j] for i, j in zip(idx_r, idx_c)])
    assert np.all(~np.isinf(a))
    assert np.all(~np.isnan(a))
    assert np.all(a >= 0)
    assert np.all(a <= 1)
    return a


def get_symmetric_matrix_from_upper_triangular(G: int, flatten_vector: np.ndarray):
    assert flatten_vector.size == G*(G-1)/2
    idx = np.triu_indices(G, 1)
    a = np.zeros((G, G))
    a[idx] = flatten_vector
    a += a.T
    a[np.diag_indices(G)] = np.inf  # put infinities on diagonal to protect caller from using these zero entries
    return a


def get_accuracy(adj_a: np.ndarray, adj_b: np.ndarray):
    G = adj_a.shape[0]
    idx = np.triu_indices(G, 1)
    accuracy = (adj_a[idx] == adj_b[idx]).sum() / idx[0].size
    assert 0 <= accuracy <= 1
    return accuracy


def bootstrap_hypothesis_test(n_bootstrap: int, data: np.ndarray, alpha: float, grid_points_in_search: int,
                              ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Perform hypothesis test using bootstrap
    :param n_bootstrap: number of bootstrap samples
    :param data: expressions data. Assumed to be G X N where G number of genes and N number of cells
    :param alpha: significance level for test.
    :param grid_points_in_search: number of grid points when searching over shift. Higher values increase runtime
    :return: adjacency matrix of size G X G with binary entries. =1 significant pair, 0 otherwise
        Also return qvalue matrix (G X G) and cost matric (G X G)
    """
    G, N = data.shape
    n_permutations = tf.constant(np.array(n_bootstrap), dtype=PRECISION_int)
    psi_ng = np.zeros((G, G)) * np.inf
    cost_ng = np.ones((G, G)) * np.inf
    cost_permuted = tf.Variable(tf.ones((G, G, n_permutations), dtype=PRECISION_fp) * tf.constant(np.inf, dtype=PRECISION_fp))
    pvalues = tf.Variable(tf.ones((G, G), dtype=PRECISION_fp) * tf.constant(np.inf, dtype=PRECISION_fp))

    candidate_psi = np.linspace(0, 2 * np.pi, grid_points_in_search)

    t=time.time()
    numpy_find_best_psi_for_each_gene_pair(psi_ng, cost_ng, data, candidate_psi=candidate_psi)
    print(f'find_best_psi_for_each_gene_pair {time.time()-t:.0f} secs')

    cost_permuted = cost_permuted.numpy()
    #data_tf = data_tf.numpy()
    #candidate_psi = candidate_psi.numpy()
    n_permutations = n_permutations.numpy()

    t = time.time()
    numpy_get_permuted_cost(cost_permuted, data, candidate_psi, n_permutations)
    print(f'get_permuted_cost {time.time()-t:.0f} secs')
     
    #cost_ng = tf.convert_to_tensor(cost_ng,dtype=tf.float32)
    #cost_permuted = tf.convert_to_tensor(cost_permuted,dtype=tf.float32)
    t = time.time()
    pvalues = numpy_get_pvalues(cost_ng,cost_permuted)
    #get_pvalues(pvalues, cost_ng, cost_permuted)
    print(f'get_pvalues {time.time()-t:.0f} secs')
    # then q-values
    # then check we find the right pair
    pvalue_flatten = flatten_upper_triangular(pvalues)
    assert pvalue_flatten.size == G * (G-1) / 2, 'only upper triangular should be save'
    qvalues_flatten, pi0 = qvalue.estimate(pvalue_flatten, verbose=True)
    qvalues = get_symmetric_matrix_from_upper_triangular(pvalues.shape[0], qvalues_flatten)
    np.allclose(qvalues, qvalues.T, atol=1e-7)
    adjacency_matrix = qvalues < alpha

    return adjacency_matrix, qvalues, cost_ng


def _calculate_error_rates(adjMatrixEstimated, adjMatrixTrue):
    TP = np.logical_and(adjMatrixEstimated, adjMatrixTrue).sum()  # true positive
    TN = np.logical_and(np.logical_not(adjMatrixEstimated), np.logical_not(adjMatrixTrue)).sum()
    FP = np.logical_and(adjMatrixEstimated, np.logical_not(adjMatrixTrue)).sum()  # false positive
    FN = np.logical_and(np.logical_not(adjMatrixEstimated), adjMatrixTrue).sum()
    # True positive rate -TP divided by number of true oscillating pairs
    TPR = TP / float(adjMatrixTrue.sum())
    # False discovery rate -FP divided by number of gene pairs reported
    FPR = FP / (FP+TN)
    if(adjMatrixEstimated.sum() == 0):
        # No co-osc genes found - so FDR is 0
        FDR = 0.
        FPR = 0.
    else:
        FDR = FP / float(adjMatrixEstimated.sum())
    assert np.all(~np.isnan(TPR))
    assert np.all(~np.isnan(FDR))
    return TPR, FPR, FDR, TP, TN, FP, FN


def get_metrics_for_different_qvalue_thresholds(qvalues: np.ndarray,
                                                adjMatrix_true: np.ndarray,
                                                alpha_values: np.ndarray
                                                ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Get true positive, false discovery and false positive rates
    :param qvalues: qvalue G X G matrix where G=number of genes
    :param adjMatrix_true: GXG binary adjacency matrix of true gene pairs
    :param alpha_values: one dimensional vector for threshold values
    :return: TPR, FDR, FPR vector of same size as alpha_values
    """
    # FDR and TPR for different thresholds
    G = adjMatrix_true.shape[0]
    TPR, FDR, FPR = np.zeros(alpha_values.size), np.zeros(alpha_values.size), np.zeros(alpha_values.size)
    TP, TN, FP, FN = np.zeros(alpha_values.size), np.zeros(alpha_values.size), np.zeros(alpha_values.size), \
                     np.zeros(alpha_values.size)
    for iaT, aT in enumerate(alpha_values):
        adjMatrixBootstrapQvalue = np.zeros((G, G), dtype=bool)
        adjMatrixBootstrapQvalue[qvalues < aT] = True
        TPR[iaT], FPR[iaT], FDR[iaT], TP[iaT], TN[iaT], FP[iaT], FN[iaT] = _calculate_error_rates(
            adjMatrixBootstrapQvalue, adjMatrix_true)
    return TPR, FDR, FPR


if __name__ == '__main__':
    t_start = time.time()
    parser = argparse.ArgumentParser(description="My parser")
    parser.add_argument('--test_mode', dest='test_mode', default=False, action='store_true')
    args = parser.parse_args()

    if args.test_mode:
        # fast example
        NG = 5
        G = 20
        N = 1000
        n_bootstrap = 100
        ngroups = 1
        grid_points_in_search = 10
    else:
        # full example
        NG = 30
        G = 1000
        N = 100
        n_bootstrap = 1000
        ngroups = 3
        grid_points_in_search = 30
    # start
    alpha = 0.001  # significance level
    data_df, phaseG, angularSpeed = GetSimISyntheticData(NG=NG, G=G, ngroups=ngroups, N=N, noiseLevel=0)

    adjacency_matrix, qvalues = bootstrap_hypothesis_test(n_bootstrap, data_df.values, alpha=alpha,
                                                          grid_points_in_search=grid_points_in_search)

    adjacency_matrix_true = true_adj_matrix(G, angularSpeed)
    correct_ratio = get_accuracy(adjacency_matrix, adjacency_matrix_true)
    print(f'Ratio of correctly identified pairs {correct_ratio:.2f}')

    TPR, FDR, _ = get_metrics_for_different_qvalue_thresholds(qvalues, adjacency_matrix_true, np.array([alpha]))
    print(f'True positive rate {float(TPR):.2f}, False discovery rate {float(FDR):.2f}')
    print(f'Total required time {time.time()-t_start:.1f} seconds')
