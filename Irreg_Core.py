import numpy as np


def IrregAllanPrep(times, values, taumin=-1.):
    diff = np.diff(times)
    if taumin <= 0.:
        taumin = float("%.2e" % (np.min(diff)))  # Find a better way to round
    else:
        taumin = float(taumin)
    splitlist = [[] for i in xrange(int(np.ceil(times[-1] / taumin)) + 1)]
    for i, stamp in enumerate(times):
        splitlist[int(np.floor(stamp / taumin))].insert(1, [times[i], values[i]])
    return splitlist, taumin

