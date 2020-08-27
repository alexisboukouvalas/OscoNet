Nose levels are 
[0.05, 0.1, 0.2, 0.3, 0.4]


The adjMatrixTrue is the true adjacency matrix
The adjMatrixQvalue has been estimated using the bootstrap oscope method at alpha = 1e-4 significance level.

Code
=====================================================
To generate the data yourself see OscopeBootstrap code SyntheticRun.py which calls
data, phaseG, angularSpeed = SyntheticDataset.GetSimISyntheticData(fPlot=False, NG=NG, G=G,
                                                                           N=100, noiseLevel=noiselevel, fReturnTruth=True)

Description of the synthetic data
============================================
    1,000 genes and 100 cells.  90 out of the 1,000 genes
    were  simulated  as  oscillators.
    The  90  oscillators  were  simulated  in  3  frequency  groups,  each
    group contains 60 genes.

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

