from matplotlib import pyplot as plt
import pickle as pickle


def plot(i, strTitle, strExperiment):
    plt.figure()
    plt.title(strTitle + ' ' + strExperiment)
    linList = []
    lgl = []
    for ip in range(pairwiseMetrics.shape[1]):
        if ip == pairwiseMetrics.shape[1]-1:
            leg = 'Oscope'
        else:
            leg = str(nbootsamples[ip])
        l, = plt.plot(noiseLevelSigma, pairwiseMetrics[:, ip, i],  '--o')
        linList.append(l)
        lgl.append(leg)
    plt.legend(lgl)
    plt.savefig(strDir+'/'+strTitle+'.png', bbox_inches='tight')


# TODO: plotting?
# def plot_batch_experiment(results, q_lower=1-0.999999, q_upper=0.999999, n_q=20):
#     import matplotlib.pyplot as plt
#     qtTry = np.linspace(q_lower, q_upper, n_q)
#     f, ax = plt.subplots(1, 1)
#     ax.set_xlabel('threshold', fontsize=20)
#     ax.plot([0, 1], [0, 1], '--k', lw=5)
#     ax.set_xlim(0.0, 0.1)
#     ax.set_ylim(0.0, 0.1)
#
#     for noiselevel, result in enumerate(results):
#         adjMatrix_true = true_adj_matrix(result['G'], result['angularSpeed'])
#         print('**** Noise level ', noiselevel)
#         for q in qtTry:
#             adj = get_adjacency_matrices(result['G'], result['angularSpeed'], result['pvalues'], alpha=q,
#                                          fBonferonni=True)
#             incorrect_ratio = calculate_inaccuracy(adj['adjMatrixBootstrapQvalue'], adj['adjMatrixTrue'])
#             print('Noise level', noiselevel, ' q=', q, 'accuracy', 1-incorrect_ratio)
#         TPR, FDR, FPR = get_metrics_for_different_qvalue_thresholds(result['qvalues'], adjMatrix_true, qtTry)
#         ax.plot(qtTry, FDR, ':', label=noiselevel, lw=5)
#     ax.legend(loc='lower right', ncol=1)
#     plt.show(block=True)


if __name__ == '__main__':
    plt.close("all")
    strDir = '/home/mqbssaby/Periodicity/Oscope/data/bootstrapRunCPU16'
    print('Plotting results after loading processcluster file from ' + strDir)
    saveDict = pickle.load(open(strDir + '/processClusterResultingPlotting.p', "rb"))
    pairwiseMetrics = saveDict['pairwiseMetrics']
    fBonferonni = saveDict['fBonferonni']
    alpha = saveDict['alpha']
    noiseLevelSigma = saveDict['noiseLevelSigma']
    nbootsamples = saveDict['nbootsamples']

    strExperiment = 'Bonferonni=%g, alpha=%g ' % (fBonferonni, alpha)
    plt.ion()

    plot(0, 'TPR', strExperiment)
    plot(1, 'FDR', strExperiment)
#     plot(2, 'TP', strExperiment)
#     plot(3, 'TN', strExperiment)
#     plot(4, 'FP', strExperiment)
#     plot(5, 'FN', strExperiment)
