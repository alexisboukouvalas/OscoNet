import numpy as np
from numpy.matlib import repmat
from matplotlib import pyplot as plt
import pandas as pd

from OscopeBootstrap.pseudotime import estimate_pseudotime_using_spectral_embedding, calculate_metrics, \
    plot_latent_space, plot_gene_fits, plot_correspondence_of_peaktime_and_times


def reproduce_figures(configuration: int, path_to_results : str = 'Results'):
    if configuration == 1:
        resultsTraining = 'cluster1'
        resultsTest = 'cluster1'
        n_neighbors = 33
    else:
        resultsTraining = 'cluster3'
        resultsTest = 'cluster1'
        n_neighbors = 31
    assert resultsTraining in ['cluster1', 'cluster3']
    assert resultsTest in ['cluster1', 'cluster3']

    np.random.seed(42)

    plt.ion()
    plt.close('all')
    data = pd.read_csv(f'{path_to_results}/filt0.9Data.txt', index_col=[0], sep=None)
    commData = pd.read_csv(f'{path_to_results}/witfildCommG50.txt', sep=None)

    # get cluster 265 genes, 239 are CC
    commId = 124 # community to look at
    comm = commData[commData.CommunityID == commId][['CloneID', 'GeneSymbol']]
    assert comm.shape[0] == 265, comm.shape

    # Calculate pseudotime time periodic for truth
    ptTruePeriod = repmat(np.arange(15), 1, 4).flatten()[:48]
    assert ptTruePeriod.size == 48, 'number of time points must be 48'

    ####################################
    # Load ENI results
    ####################################
    enid = '_30k'
    clusterENI = {}
    for i in range(1, 4):
        cloneids = pd.read_csv(f'{path_to_results}/OscopeClusterC%i%s.csv' % (i, enid), index_col=[0])
        eniIndex = pd.read_csv(f'{path_to_results}/OscopeClusterC%iENI%s.csv' % (i, enid), index_col=[0])
        eni=(eniIndex-1) / (eniIndex-1).max()
        print('%g: Nk=%g, eni:%g' % (i, cloneids.shape[0], eni.shape[0]), cloneids.head(), eni.head())
        clusterENI['cluster%g' % i] = {'pt': eni.values.flatten(), 'cloneid': cloneids}

    ####################################
    # Data and parameters
    ####################################
    cloneidlist = clusterENI[resultsTraining]['cloneid'].values.flatten()
    # Training data
    training_data = data.loc[cloneidlist, :].values.T
    assert training_data.shape[1] == cloneidlist.size
    assert training_data.shape[0] == 48

    cloneidlistEval = clusterENI[resultsTest]['cloneid'].values.flatten()  # use just cluster 1 genes for evaluation

    ####################################
    # Results in paper
    ####################################

    geneplot = ['PLK', 'CCNE1', 'CDC6']

    # look at fits using the true time (ptTruePeriod)
    metrics_true = calculate_metrics(data, commData, cloneidlistEval, ptTruePeriod, ptTruePeriod)
    _ = plot_gene_fits(geneplot, data, cloneidlistEval, commData, ptTruePeriod, ptTruePeriod)

    # spectral method
    pt_spectral, latent_space_2d = estimate_pseudotime_using_spectral_embedding(training_data, n_neighbors)
    metrics_spectral = calculate_metrics(data, commData, cloneidlistEval, pt_spectral, ptTruePeriod)
    _ = plot_latent_space(latent_space_2d, ptTruePeriod)
    _ = plot_gene_fits(geneplot, data, cloneidlistEval, commData, pt_spectral, ptTruePeriod)
    _ = plot_correspondence_of_peaktime_and_times(metrics_spectral['peakTimes'], metrics_true['peakTimes'], ptTruePeriod, pt_spectral)

    # eni fit (ptENI)
    ptENI = clusterENI[resultsTraining]['pt']
    metrics_eni = calculate_metrics(data, commData, cloneidlistEval, ptENI, ptTruePeriod)
    _ = plot_gene_fits(geneplot, data, cloneidlistEval, commData, ptENI, ptTruePeriod)
    _ = plot_correspondence_of_peaktime_and_times(metrics_eni['peakTimes'], metrics_true['peakTimes'], ptTruePeriod, ptENI)

if __name__ == "__main__":
    reproduce_figures()
