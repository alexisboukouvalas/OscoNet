"""
Estimate pseudotime using spectral embedding methods.
We also include plotting and error calculation functions .
"""
import numpy as np
from matplotlib import pyplot as plt
from scipy import optimize
from scipy.interpolate import UnivariateSpline
from scipy.stats import spearmanr
from sklearn import manifold
from sklearn import preprocessing


def estimate_pseudotime_using_spectral_embedding(data, n_neighbors: int):
    """
    :param data: gene data
    :param n_neighbors: parameter for SpectralEmbedding
    :return: pseudotime and 2-D latent space
    """
    X_init_embed = manifold.SpectralEmbedding(n_components=2, n_neighbors=n_neighbors).fit_transform(data)
    X_init = preprocessing.scale(X_init_embed)  # 0 mean unit variance
    # Fit circle
    r = _get_circle(X_init)
    X_init -= r['center']  # remove mean
    ptRadians = np.arctan2(X_init[:, 1], X_init[:, 0])
    pt = (ptRadians - ptRadians.min()) / (ptRadians.max() - ptRadians.min())  # convert to [0, 1]
    return pt, X_init


def calculate_metrics(data, commData, cloneidlistEval, pt, ptTruePeriod):
    # Measure fit
    peakTimes = np.zeros(cloneidlistEval.size)
    roughness = np.zeros(cloneidlistEval.size)

    peakTimes_true = np.zeros(cloneidlistEval.size)
    roughness_true = np.zeros(cloneidlistEval.size)

    for ic, c in enumerate(cloneidlistEval):
        a = commData[commData.CloneID == c]
        peakTimes[ic] = _get_peak_time(data, a.CloneID, pt)
        peakTimes_true[ic] = _get_peak_time(data, a.CloneID, ptTruePeriod)
        roughness[ic] = _calc_roughness(data.loc[a.CloneID, :].values.flatten(), pt)
        roughness_true[ic] = _calc_roughness(data.loc[a.CloneID, :].values.flatten(), ptTruePeriod)
    # KS test to uniformity
    from scipy import stats
    ptKS = stats.kstest(pt, 'uniform')[0]

    ks, pseudoCorr, peakCorr = list(), list(), list()
    ks.append(ptKS)
    pseudoCorr.append(spearmanr(pt, ptTruePeriod)[0])
    peakCorr.append(spearmanr(peakTimes, peakTimes_true)[0])
    print('peak R=%.2f' % peakCorr[-1],
          'pseudotime - true roughness %.2f' % (np.median(roughness - np.median(roughness_true))),
          'corr pseudotime %.2f, ks=%.2f' % (pseudoCorr[-1], ptKS))
    return dict(peakTimes=peakTimes, roughness=roughness, ptKS=ptKS, ks=ks,
                pseudoCorr=pseudoCorr, peakCorr=peakCorr)


def _calc_R(xc, yc, dataC):
    """ calculate the distance of each 2D points from the center (xc, yc) """
    return np.sqrt((dataC[:, 0]-xc)**2 + (dataC[:, 1]-yc)**2)


def _f_2(c, dataC):
    """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
    Ri = _calc_R(*c, dataC)
    return Ri - Ri.mean()


def _get_circle(dataC):
    ''' Calculate circle using lease circles'''
    center_estimate = np.mean(dataC[:, 0]), np.mean(dataC[:, 1])
    center_2, ier = optimize.leastsq(_f_2, center_estimate, args=dataC)
    Ri_2       = _calc_R(*center_2, dataC)
    R_2        = Ri_2.mean()
    residu_2   = sum((Ri_2 - R_2)**2)
    return {'center': center_2, 'radius': R_2, 'residual': residu_2}


def _get_peak_time(data, cloneid, pt, baseCycleTime=None):
    d = data.loc[cloneid, :].values.flatten()
    if baseCycleTime is None:
        baseCycle = d
    else:
        baseCycle = d[baseCycleTime[0]:baseCycleTime[1]]
    return pt[np.argmax(baseCycle)]


def _calc_roughness(x, pt):
    ''' evaluate roughness of pseudotime
    This metric measures the smoothness of the gene expression profile by looking at the differences
    of consecutive measurements.
    Smaller values indicate a smoother response. '''
    x=np.atleast_2d(x)
    i=np.argsort(pt)
    x = x[:,i]
    N = x.shape[1]
    assert(N > 0)
    S = x.std(axis=1)
    return np.sqrt(np.sum(np.square(x[:,0:(N-1)] - x[:, 1:N]),1) / (N-1)) / S

####################
# Plotting
####################
def plot_latent_space(x_latent, ptTruePeriod):
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    ax.scatter(x_latent[:, 0], x_latent[:, 1], c=ptTruePeriod, s=35)
    r = _get_circle(x_latent)
    circle1 = plt.Circle(r['center'], r['radius'], color='r', alpha=0.1)
    ax.add_artist(circle1)
    return fig, ax


def plot_gene_fits(geneplot, data, cloneidlistEval, commData, pt, ptTruePeriod):
    # Plot gene fits
    if len(geneplot) < 5:
        f, ax = plt.subplots(int(len(geneplot)), 1, figsize=(7, 7), sharey=False, sharex=False)
    else:
        f, ax = plt.subplots(int(np.sqrt(len(geneplot))), int(np.sqrt(len(geneplot))),
                             figsize=(20, 20), sharey=False, sharex=False)
    ax = ax.flatten()
    iterC = 0
    for ic, c in enumerate(cloneidlistEval):
        a = commData[commData.CloneID == c]
        assert a.shape[0] == 1, 'must find exactly one match for cloneid %s' % c
        g = a.GeneSymbol.values[0]
        if g in geneplot:
            _plotgene(data, ptTruePeriod, a.CloneID, ax[iterC], pt)
            ax[iterC].set_title(g)
            # add cubic interpolator
            d = data.loc[a.CloneID, :].values.flatten()
            idx = np.argsort(pt)
            fspline = UnivariateSpline(pt[idx], d[idx], check_finite=True)
            xtest = np.linspace(np.min(pt), np.max(pt), 100)
            fs = fspline(xtest)
            assert not np.any(np.isnan(fs)), 'spline %s-%s' % (xtest, f)
            ax[iterC].plot(xtest, fs, color='m', lw=3)
            iterC += 1
    return f, ax


def plot_correspondence_of_peaktime_and_times(peakTimes_estimate, peakTimes_true, ptTruePeriod, pt):
    f, ax = plt.subplots(1, 1, figsize=(7, 7))
    _plot_correspondance(ax, peakTimes_estimate, peakTimes_true)
    ax.set_ylabel('True peak time')
    ax.set_xlabel('Estimated Peak time')
    f, ax = plt.subplots(1, 1, figsize=(7, 7))
    _plot_correspondance(ax, pt, ptTruePeriod, c=ptTruePeriod)
    ax.set_ylabel('True time')
    ax.set_xlabel('Estimated time')
    return f, ax


def _plotgene(data, ptTruePeriod, cloneid, ax, pt, baseCycleTime=None, fFocus=False, f=None):
    d = data.loc[cloneid, :].values.flatten()
    peakTime = _get_peak_time(data, cloneid, pt, baseCycleTime)
    if baseCycleTime is not None:
        ax.axvspan(baseCycleTime[0], baseCycleTime[1], color='blue', alpha=0.5)
        if fFocus:
            ax.set_xlim(baseCycleTime)
    h = ax.scatter(pt, d, label=cloneid, s=35, c=ptTruePeriod)  # plot gene expression
    ax.axvline(peakTime, lw=5, color='green')
    if f is not None:
        f.colorbar(h, ticks=np.arange(15))
    return peakTime


def _plot_correspondance(ax, a, b, c=None):
    ax.scatter(a, b, s=100, c=c)
    ax.set_title('Spearmanr correlation %.2f' % spearmanr(a, b)[0])
    v = ax.axis()
    ax.plot( [np.min(v[:2]), np.max(v[:2])], [np.min(v[-2:]), np.max(v[-2:])], '-k', lw=3, alpha=0.5)
    ax.plot( [np.max(v[:2]), np.min(v[:2])], [np.min(v[-2:]), np.max(v[-2:])], '-k', lw=3, alpha=0.5)
