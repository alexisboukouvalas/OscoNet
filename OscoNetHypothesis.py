from __future__ import print_function
import numpy as np
import pandas as pd
import numpy as np
import OscopeVectorize
import pandas as pd
from joblib import Parallel, delayed
import time
import scipy.stats
import pickle
import multiprocessing
import os


def GetSimISyntheticData(fPlot=False, NG=15, G=1000, N=100, noiseLevel=0, fReturnTruth=False, ngroups=3):
    '''
    Sim I: Oscope paper supplementary
    1,000 genes and 100 cells.  90 out of the 1,000 genes
    were  simulated  as  oscillators.
    The  90  oscillators  were  simulated  in  3  frequency  groups,  each
    group contains 30 genes.

     Group 1 and 3 following the same  order,  while  genes  in  group  2  following  another  order.   In
    Sim  I ,  the  relative  speeds  of the  3  groups  are  proportional  to  2:3:6.

    Within each frequency group, genes were further simulated with strong and weak signals.
    Half of the oscillatory genes were simulated as strong oscillators with sigma_g = sigm_str . The other half
    were simulated as weak oscillators with sigma_g = sigma_wk = 2sigma_str .

    Starting phase phi_g varies in different genes within a frequency group.

    The remaining genes except the oscillators are called noise
    genes. Noise genes were simulated as random Gaussian noise. The noise level was adjusted to
    be comparable to the average noise signal among all oscillators.

    Simulation study
    the sigma_str varies from 0.05 to 0.4 in 5 steps.
    '''
    assert ngroups <= 3, 'Only 3 groups implemented'
    # Construct oscillatory groups
    sigma_strLevel = [0.05, 0.1, 0.2, 0.3, 0.4, 0.6]
    sigma_str = sigma_strLevel[noiseLevel]

    # two different orders
    t1 = np.linspace(0, 2*np.pi, N)
    t2 = np.random.permutation(t1)

    data = np.zeros((G, N))
    data[:] = np.nan
    # genes per weak/strong oscillatory group
    # Group 1
    cellName = []
    for i in range(N):
        cellName.append('C'+str(i))

    geneName = []

    phaseG = np.zeros((G))
    angularSpeed = np.zeros((G))

    for i in range(NG):  # strong oscillators
        startingPhase = np.random.uniform(0, 2*np.pi)
        phaseG[i] = startingPhase
        angularSpeed[i] = 2
        data[i, :] = np.sin(2*t1 + startingPhase) + sigma_str*np.random.randn(N)
        geneName.append('G1SO'+str(i))
    for i in range(NG, 2*NG):  # weak oscillators
        startingPhase = np.random.uniform(0, 2*np.pi)
        phaseG[i] = startingPhase
        angularSpeed[i] = 2
        data[i, :] = np.sin(2*t1 + startingPhase) + 2*sigma_str*np.random.randn(N)
        geneName.append('G1WO'+str(i))
    if(ngroups >= 2):
        # Group 2
        for i in range(2*NG, 3*NG):  # strong oscillators
            startingPhase = np.random.uniform(0, 2*np.pi)
            phaseG[i] = startingPhase
            angularSpeed[i] = 3
            data[i, :] = np.sin(3*t2 + startingPhase) + sigma_str*np.random.randn(N)
            geneName.append('G2SO'+str(i))
        for i in range(3*NG, 4*NG):  # weak oscillators
            startingPhase = np.random.uniform(0, 2*np.pi)
            phaseG[i] = startingPhase
            angularSpeed[i] = 3
            data[i, :] = np.sin(3*t2 + startingPhase) + 2*sigma_str*np.random.randn(N)
            geneName.append('G2WO'+str(i))
    if(ngroups >= 3):
        # Group 3
        for i in range(4*NG, 5*NG):  # strong oscillators
            startingPhase = np.random.uniform(0, 2*np.pi)
            phaseG[i] = startingPhase
            angularSpeed[i] = 6
            data[i, :] = np.sin(6*t1 + startingPhase) + sigma_str*np.random.randn(N)
            geneName.append('G3SO'+str(i))
        for i in range(5*NG, 6*NG):  # weak oscillators
            startingPhase = np.random.uniform(0, 2*np.pi)
            phaseG[i] = startingPhase
            angularSpeed[i] = 6
            data[i, :] = np.sin(6*t1 + startingPhase) + 2*sigma_str*np.random.randn(N)
            geneName.append('G3WO'+str(i))

    # white noise genes
    for w in range(i+1, G):  # use i index from above where it stopped
        phaseG[w] = np.nan
        angularSpeed[w] = np.nan
        data[w, :] = np.max([3/2 * sigma_str, 1]) * np.random.randn(N)
        geneName.append('R'+str(w))

    assert np.all(~ np.isnan(data)), 'Entries with nans!'
    assert len(geneName) == G
    assert len(cellName) == N

    if(fPlot):
        from matplotlib import pyplot as plt
        plt.ion()
        _, axList = plt.subplots(7, sharex=True, sharey=True)
        for c, i in enumerate([1, 15, 30, 45, 60, 75, 90]):
            if i == 30 or i == 45:
                t = t2
            else:
                t = t1
            axList.flatten()[c].plot(t, data[i, :], 'bo')

    df = pd.DataFrame(data, index=geneName, columns=cellName)

    if(fReturnTruth):
        return df, phaseG, angularSpeed  # return GXN matrix
    else:
        return df


def CreateEdgeNetwork(adjMatrixBootstrap, cost, psi, geneNames):
    '''
        CreateEdgeNetwork - Create Edge file
    '''
    assert np.all(adjMatrixBootstrap.shape == cost.shape)
    # we remove significant pairs that are not symmetric
    print('Is adjacency matrix symmetric?', np.allclose(adjMatrixBootstrap, adjMatrixBootstrap.T))
    adjMatrixBootstrap = np.logical_and(adjMatrixBootstrap, adjMatrixBootstrap.T)
    assert np.allclose(adjMatrixBootstrap, adjMatrixBootstrap.T)
    print('Adjacency matrix now symmetric after applying AND transformation')
    G = cost.shape[0]
    nt = G*(G-1)  # number of tests without diagonal
    print('Sparseness %f' % (adjMatrixBootstrap.sum() / float(nt)))

    # Create edge representation
    # G_i, G_j, cost for all significant genes
    nSignificantPairs = adjMatrixBootstrap.sum() / 2.  # symmetric matrix
    assert(nSignificantPairs.is_integer())
    edgeNetwork = []  # np.empty((int(nSignificantPairs), 3), dtype='string, string, float64')
    iterC = 0
    for i in range(G):
        for j in range(i+1, G):
            if(adjMatrixBootstrap[i, j] == 1):
                assert np.allclose(psi[i, j], psi[j, i])
                assert np.allclose(cost[i, j], cost[j, i])
                edgeNetwork.append([geneNames[i], geneNames[j], cost[i, j], psi[i, j]])
                iterC += 1
    a = pd.DataFrame(data=edgeNetwork, columns=['gene1', 'gene2', 'cost', 'psi'])
    # then rank by cost - distance, the bigger the worse - try to minimize, so order is ascending
    b = a.sort_values('cost')
    return b


def ensure_dir(f):
    ''' from http://stackoverflow.com/questions/273192/how-to-check-if-a-directory-exists-and-create-it-if-necessary '''
    assert len(f) > 0, 'given bad file name %s' % f
    d = os.path.dirname(f)
    if not os.path.exists(d):
        os.makedirs(d)

# If we don't want to use rpy to call R Oscope
def OscopePreprocessing(p):
    '''
    Replaces outliers with quantile (5-95%) and rescales to [-1, 1]
    '''
    assert isinstance(p, pd.DataFrame), 'Must be panda df as Oscope requires column (cell)  and row (gene) names'
    # Remove genes with unique values len == 0 - they are not oscillating
    gf = np.array([np.unique(p.loc[g,:]).size > 1 for g in p.index])
    print('Removing genes', p.index[gf == False])
    data = p.loc[gf, :].copy()

    Q5 = data.quantile(q=0.01, axis=1)[:, None]
    Q95 = data.quantile(q=0.99, axis=1)[:, None]
    Rg = Q95 - Q5
    if np.any(Rg == 0):
        gl = data.index[np.flatnonzero(Rg==0)]
        print('Unique values for genes with singular quantile range:')
        for g in gl:
            print(g, np.unique(data.loc[g, :]))
        raise NameError('Range must be bigger than 0 - the following genes have 0 range %s.' % str(gl))
    rescaledData = ((data - np.tile(Q5,data.shape[1])) * 2 / np.tile(Rg, data.shape[1])) - 1
    rescaledData[rescaledData < (-1)] = -1
    rescaledData[rescaledData > 1] = 1

    assert np.all(rescaledData >= -1)
    assert np.all(rescaledData <= 1)

    return rescaledData


def RunBootstrapPercentile(data, NumPerm=2000, gN=30, useTF=True,
                           n_cores= multiprocessing.cpu_count(),strSave='test', fSavefile=True):
    # Uses globals
    G, N = data.shape
    assert np.isnan(data.values).sum() == 0
    t0 = time.time()
    pvalues, psi, cost = DoTheBootstrapFast(data, pWorkers=n_cores, NumPerm=NumPerm, gN=gN, useTF=useTF, fDebug=True)
    print(('Bootstrap parallel time with %g workers took %.1f secs. Bootstrap samples %g. Use TF %g.' %
          (n_cores, time.time()-t0, NumPerm, useTF)))
    saveDict = {'pvalues': pvalues, 'n_cores': n_cores, 'NumPerm': NumPerm,
                'psi': psi, 'cost': cost, 'G': G, 'N': N, 'gN': gN, 'useTF': useTF}
    if(fSavefile):
        f = '%s_SummaryPartition.p' % strSave
        pickle.dump(saveDict, open(f, "wb"))
    return saveDict



def DoTheBootstrap(data, pWorkers=2, NumPerm=200, percentileTest=99, fDebug=False, listOfTargetGenes=None):
        assert isinstance(data, pd.DataFrame)
        assert percentileTest > 0 and percentileTest < 100
        if(listOfTargetGenes is None):
            listOfTargetGenes = data.index  # use all genes
        print('DoTheBootstrap:: Running on %g  cores with %g target genes.' % (pWorkers, len(listOfTargetGenes)))
        # Calculate similarity target-all candidate genes
        G, N = data.shape
        psi, cost = OscopeVectorize.OptimizeByGridSearch_np(data.values, gN=50)
        idxSignificantGenes = np.zeros((G, G), dtype=np.bool_)
        if(pWorkers == 1):
            for t in range(len(listOfTargetGenes)):
                idxSignificantGenes[t, :] = SingleHypothesisTF(m, t, psi, cost, data.values, NumPerm, percentileTest)
        else:
            idx = Parallel(n_jobs=pWorkers)(delayed(SingleHypothesisNP)(t, psi, cost, data.values, NumPerm, percentileTest)
                                            for t in range(len(listOfTargetGenes)))
            idxSignificantGenes = np.array(idx)  # convert list
        return idxSignificantGenes  # Return G index with boolean index of significant genes


def SingleHypothesisTF(m, t, psi, cost, dataMatrix, NumPerm=200, percentileTest=99):
    '''
    t is gene of interest
    Permute all other genes.
    '''
    G = dataMatrix.shape[0]
    assert t >= 0 and t <= G  # target gene index
    # perform hypothesis test by permuting cells
    bootstrapCost = np.zeros((G, NumPerm))
    for ip in range(NumPerm):
        # shuffle all the cells except for the target gene - for this we need to make a copy of the array
        # Pass is N X G, this will permute first axis, i.e. cells, creates a copy
        permuted = np.random.permutation(dataMatrix.T).T
        permuted[t, :] = dataMatrix[t, :]  # overwrite true order
        costP = m.EvalE(permuted, psi)
        # costP is a G X G matrix - we just want  the target gene
        bootstrapCost[:, ip] = costP[t, :]  # permuted cost
    # Hypothesis test
    acSim = -np.log10(cost[t, :])  # actual similarity
    bootSim = -np.log10(bootstrapCost)
    qBoot = np.percentile(bootSim, percentileTest, 1)  # along all cells
    return qBoot < acSim  # return boolean index


def SingleHypothesisNP(t, psi, cost, dataMatrix, NumPerm=200, percentileTest=99):
    '''
    idxSignificantGenes is output
    '''
    G = dataMatrix.shape[0]
    assert t >= 0 and t <= G  # target gene index
    # perform hypothesis test by permuting cells
    bootstrapCost = np.zeros((G, NumPerm))
    for ip in range(NumPerm):
        # shuffle all the cells except for the target gene - for this we need to make a copy of the array
        # Pass is N X G, this will permute first axis, i.e. cells, creates a copy
        permuted = np.random.permutation(dataMatrix.T).T
        permuted[t, :] = dataMatrix[t, :]  # overwrite true order
        costP = OscopeVectorize.EvalE2_np(permuted, psi)
        bootstrapCost[:, ip] = costP[t, :]  # permuted cost
    # Hypothesis test
    acSim = -np.log10(cost[t, :])  # actual similarity
    bootSim = -np.log10(bootstrapCost)
    qBoot = np.percentile(bootSim, percentileTest, 1)  # along all bootstrap samples
    return qBoot < acSim


def GeneratePvalues(X, gN=50):
    '''
    X is a G X N matrix, randomly permute each column independently.
    '''
    G, N = X.shape
    XPermuted = np.empty_like(X)
    for g in range(G):
        idxCell = np.random.permutation(N)  # permuted cell order
        XPermuted[g, :] = X[g, idxCell]

    _, cost = OscopeVectorize.OptimizeByGridSearch_np(XPermuted, gN=gN)
    return cost


def DoTheBootstrapFast(data, pWorkers=2, NumPerm=200, fDebug=False,
                       listOfTargetGenes=None, gN=50, fSinglePassRandom=False, useTF=True):
        assert isinstance(data, pd.DataFrame)
        if(listOfTargetGenes is None):
            listOfTargetGenes = data.index  # use all genes
        print('DoTheBootstrap:: Running on %g  cores with %g target genes.' % (pWorkers, len(listOfTargetGenes)))
        # Calculate similarity target-all candidate genes
        G, N = data.shape
        if(useTF):
            if(fDebug):
                print('DoTheBootstrapFast: Initial search Using TF')
            m = OscopeVectorize.OscopeTF(G, N)  # setup TF evaluation function
            psi, cost = m.OptimizeByGridSearch(data.values, gN=gN)
        else:
            if(fDebug):
                print('DoTheBootstrapFast: Initial search Using NP')
            psi, cost = OscopeVectorize.OptimizeByGridSearch_np(data.values, gN=gN)

        if(fSinglePassRandom):
            if(fDebug):
                print('DoTheBootstrapFast: fSinglePassRandom')
            # Generate column wide random reorderings
            nullDistributionList = Parallel(n_jobs=pWorkers)(delayed(GeneratePvalues)(data.values, gN=gN)
                                                             for _ in range(NumPerm))
            nullDistribution = np.array(nullDistributionList).ravel()
            pvalueList = Parallel(n_jobs=pWorkers, max_nbytes=1e6)(delayed(SingleHypothesisCommonNull)(t, cost, nullDistribution)
                                                                   for t in range(len(listOfTargetGenes)))
#             pvalueList = Parallel(n_jobs=pWorkers, max_nbytes=1e6)(delayed(SingleHypothesisNPFast)(t, psi, cost, data.values, NumPerm, gN=gN,
#                                                                                    nullDistributionForGene=nullDistribution[:, t, :])
#                                                    for t in range(len(listOfTargetGenes)))
        else:
            if(fDebug):
                print('DoTheBootstrapFast: multiple passes per gene')
            pvalueList = Parallel(n_jobs=pWorkers, max_nbytes=1e6)(delayed(SingleHypothesisNPFast)(t, psi, cost, data.values, NumPerm, gN=gN)
                                                                   for t in range(len(listOfTargetGenes)))

        pvalues = np.array(pvalueList)  # convert list
        return pvalues, psi, cost  # Return G index with boolean index of significant genes


def SingleHypothesisCommonNull(t, cost, nullDistribution):
    '''
    Common null distribution for all
    '''
    # Hypothesis test
    pvalue = np.ones((cost.shape[0]))  # G vector of p-values
    G = cost.shape[0]
    for i in range(G):
        pvalue[i] = scipy.stats.percentileofscore(nullDistribution, cost[t, i], 'mean') / 100.
    return pvalue


def SingleHypothesisNPFast(t, psi, cost, X, NumPerm=200, gN=50, nullDistributionForGene=None):
    '''
    dataMatrix is G X N
    idxSignificantGenes is output
    '''
    G, N = X.shape
    assert t >= 0 and t <= G  # target gene index
    if nullDistributionForGene is None:
        # perform hypothesis test by permuting cells
        bootstrapCost = np.zeros((G, NumPerm))
        Y = X  # the other genes, include target gene for simplicity, no copy needed
        for ip in range(NumPerm):
            # shuffle all the cells except for the target gene - for this we need to make a copy of the array
            idxCell = np.random.permutation(N)  # permuted cell order
            costP = OscopeVectorize.EvalSingle_MultiplePsi(np.expand_dims(X[t, :], 0), Y[:, idxCell],
                                                           np.expand_dims(psi[t, :], 1), gN=gN)
            bootstrapCost[:, ip] = costP  # permuted cost
    else:
        # Already randomly permuted cost computed
        bootstrapCost = np.reshape(nullDistributionForGene, (NumPerm, G)).T

    # Hypothesis test
    pvalue = np.zeros((bootstrapCost.shape[0]))
    pvalue[:] = np.nan
    for i in range(G):
        pvalue[i] = scipy.stats.percentileofscore(bootstrapCost[i, :], cost[t, i], 'mean') / 100.

    return pvalue

def CalculateErrorRates(adjMatrixEstimated, adjMatrixTrue):
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


if __name__ == '__main__':
    data, phaseG, angularSpeed = GetSimISyntheticData(fPlot=False, NG=10, G=100, N=100, noiseLevel=0, ngroups=2, fReturnTruth=True)

    G, N = data.shape
    print('G=%g, N=%g' % (G, N))

    testddirect = OscopePreprocessing(data)
    print('postprocessed data size', testddirect.shape)
    t = time.time()
    r = RunBootstrapPercentile(testddirect, NumPerm=100, gN=G, useTF=False,
                                           n_cores=1, fSavefile=False)
    print('Completed in %g seconds' % (time.time()-t))

    ####################### Hypothesis test ############################################
    import qvalue
    alpha = 1e-6  # significance level for hypothesis test. Can vary this here.
    dj = np.diag_indices(G, 2)  # index to diagonal of a G X G matrix
    # qvalue calculation should be similar to Bonferonni
    qvalues, pi0 = qvalue.estimate(r['pvalues'], verbose=True)
    adjMatrixBootstrapQ = np.zeros((G, G), dtype=bool)
    adjMatrixBootstrapQ[qvalues < alpha] = True
    adjMatrixBootstrapQ[dj] = False  # remove self-transitions

    adjMatrixTrue = np.zeros((G, G), dtype=bool)
    for i in range(G):
        for j in range(G):
            adjMatrixTrue[i, j] = (angularSpeed[i] == angularSpeed[j])
    adjMatrixTrue[dj] = False  # have to do it, because white noise diagonal not true

    trueDF = pd.DataFrame(adjMatrixTrue, columns=data.index, index=data.index)
    trueDF.to_csv('trueClusteringAdjMatrix.csv')

    ####################### Create edge representation ############################################
    edgeNetwork = CreateEdgeNetwork(adjMatrixBootstrap=adjMatrixBootstrapQ, cost=r['cost'],
                                                    psi=r['psi'], geneNames=testddirect.index)
    print('Writing edge network to', 'toy.csv', 'number of edges', edgeNetwork.shape[0])
    edgeNetwork.to_csv('toy.csv')

    # Calculate error rates
    TPR, FPR, FDR, TP, TN, FP, FN = CalculateErrorRates(adjMatrixBootstrapQ, adjMatrixTrue)
    print('False positive rate', FPR, 'True positive rate', TPR, 'False discrovery rate', FDR)
    print('True positives=', TP, 'True negatives=', TN, 'False positives=', FP, 'False negatives', FN)
