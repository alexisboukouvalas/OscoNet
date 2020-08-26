from __future__ import print_function

from typing import Tuple

import numpy as np
import pandas as pd



def GetSimISyntheticData(NG: int = 15, G: int = 1000, N: int = 100, noiseLevel: int = 0, ngroups: int = 3) \
        -> Tuple[pd.DataFrame, np.ndarray, np.ndarray]:

    '''
    Generate synthetic data.
   :param NG: half-size of each gene group. For example =3, each group will have 6 co-oscillating genes. 3 of them
   will be strong oscillators and 3 weak (double the noise)
   :param G: number of total genes
   :param N: number of cells
   :param noiseLevel: noise level index (0 to 5 index)
   :param ngroups: number of groups
   :return: dataframe of data of shape (G X G), numpy arrays of true phase and angular
   speed, each a vector of size G

    Reproduced from the original Oscope paper
    https://www.nature.com/articles/nmeth.3549

    We include below the original description from the supplementary material of that paper:
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
        print('I am in group 2')
        print('NG = {0}'.format(NG))
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
    print('White noise from {0}'.format(i+1))
    # white noise genes
    for w in range(i+1, G):  # use i index from above where it stopped
        phaseG[w] = np.nan
        angularSpeed[w] = np.nan
        data[w, :] = np.max([3/2 * sigma_str, 1]) * np.random.randn(N)
        geneName.append('R'+str(w))

    assert np.all(~ np.isnan(data)), 'Entries with nans!'
    assert len(geneName) == G
    assert len(cellName) == N
    df = pd.DataFrame(data, index=geneName, columns=cellName)
    return df, phaseG, angularSpeed  # return GXN matrix


def true_adj_matrix(G: int, angularSpeed: np.ndarray):
    # Count significant pairs - create the real adjacency matrix
    dj = np.diag_indices(G, 2)  # index to diagonal of a G X G matrix
    adjMatrixTrue = np.zeros((G, G), dtype=bool)
    for i in range(G):
        for j in range(G):
            adjMatrixTrue[i, j] = (angularSpeed[i] == angularSpeed[j])
    adjMatrixTrue[dj] = False  # have to do it, because white noise diagonal not true
    return adjMatrixTrue
