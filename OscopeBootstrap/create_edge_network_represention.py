import numpy as np
import pandas as pd


def create_edge_network_representation(adjMatrixBootstrap, weight_matrix, gene_names):
    """
        CreateEdgeNetwork - Create Edge file.
        This is needed before hypothesis test q-value derived adjacency matrix
        can be consumed by R network analysis code.
        Return a pandas dataframe with 3 columns, two gene names for the gene-pair and the cost value
    """
    assert np.all(adjMatrixBootstrap.shape == weight_matrix.shape)
    # we remove significant pairs that are not symmetric
    assert np.allclose(adjMatrixBootstrap, adjMatrixBootstrap.T), 'not symmetric'
    G = weight_matrix.shape[0]
    nt = G*(G-1)  # number of tests without diagonal
    print('Sparseness %f' % (adjMatrixBootstrap.sum() / float(nt)))
    # Get gene names
    assert(len(gene_names) == G)

    # Create edge representation
    # G_i, G_j, cost for all significant genes
    nSignificantPairs = adjMatrixBootstrap.sum() / 2.  # symmetric matrix
    assert(nSignificantPairs.is_integer())
    edgeNetwork = []  # np.empty((int(nSignificantPairs), 3), dtype='string, string, float64')
    iterC = 0
    for i in range(G):
        for j in range(i+1, G):
            if(adjMatrixBootstrap[i, j] == 1):
                edgeNetwork.append([gene_names[i], gene_names[j], weight_matrix[i, j]])
                iterC += 1
    a = pd.DataFrame(data=edgeNetwork, columns=['gene1', 'gene2', 'weight'])

    return a
