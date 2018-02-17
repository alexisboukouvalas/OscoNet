from __future__ import print_function
import numpy as np

def OptimizeByGridSearch_np(X, gN, fDebug=False):
    ''' Find Psi for every combination by grid search '''
    assert gN % 2 == 0, 'Pass even number of grid point'
    cphi = np.linspace(0, 2*np.pi, gN)
    G = X.shape[0]
    Psi = np.zeros((G, G))
    currentBestCost = np.zeros((G, G))
    currentBestCost[:] = np.inf
    currentBestPsi = -1*np.ones((G, G))
    for g in range(0, gN, 2):  # scales with number of candidate points, not number of genes
        lt = np.tril_indices(G, -1)
        ut = np.triu_indices(G, 1)
        Psi[lt] = cphi[g]
        Psi[ut] = cphi[g+1]
        c = EvalE2_np(X, Psi)  # return G X G matrix
        b = np.stack([c, currentBestCost])  # creates a 2 X G X G matrix
        currentBestCost = np.min(b, 0)
        currentBestPsi = np.choose(np.argmin(b, 0), [Psi, currentBestPsi])
    # Post-processing step - collapse to symmetric matrix by picking best overall
    for ig in range(G):
        for jg in range(ig+1, G):
            a = currentBestCost[[ig, jg], [jg, ig]]
            currentBestCost[[ig, jg], [jg, ig]] = np.min(a)
            b = currentBestPsi[[ig, jg], [jg, ig]]
            currentBestPsi[[ig, jg], [jg, ig]] = np.choose(np.argmin(a), b)
    if fDebug:
        assert np.allclose(currentBestCost, currentBestCost.T), 'Cost matrix  must be symmetric'
        assert np.allclose(currentBestPsi, currentBestPsi.T), 'Psi matrix  must be symmetric'
    return currentBestPsi, currentBestCost


def EvalE2_np(X, Psi):
    ''' Numpy version of EvalE '''
    N = X.shape[1]
    X2 = np.square(np.expand_dims(X, 1))  # now G x 1 x N
    Y2 = np.square(np.expand_dims(X, 0))  # now 1 x G x N
    Xb = np.expand_dims(X, 1)  # now G x 1 x N
    Yb = np.expand_dims(X, 0)  # now 1 x G x N
    cosb = np.expand_dims(np.cos(Psi), 2)  # G X G X 1
    cosbt = np.tile(cosb, [1, 1, N])
    XYb = Xb * Yb  # G X G X N
    sinb = np.expand_dims(np.square(np.sin(Psi)), 2)  # G X G X 1
    sinbt = np.tile(sinb, [1, 1, N])
    return np.sum(np.abs(X2+Y2 - 2 * XYb * cosbt - sinbt), 2)  # G X G


def EvalE2_npSingle(X, Y, Psi):
    ''' Numpy version of EvalE. Has target X single gene, Y is set of genes G'''
    # X is 1 X N
    # Y is G X N
    # X*Y is GXN
    # so is X-Y
    return np.sum(np.abs(np.square(X) + np.square(Y) - 2 * X * Y * np.cos(Psi) - np.square(np.sin(Psi))), 1)  # G X 1


def EvalSingle_MultiplePsi(X, Y, Psi, gN=5):
    cphi = np.linspace(0, 2*np.pi, gN)
    f = EvalE2_npSingle(X, Y, Psi)  # best guess is to take unpermuted values
    for i in range(gN):
        fi = EvalE2_npSingle(X, Y, cphi[i])
        f = np.where(fi < f, fi, f)  # u pdate with minimum value
    return f

if __name__ == '__main__':
    np.random.seed(0)
