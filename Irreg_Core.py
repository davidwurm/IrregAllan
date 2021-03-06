import numpy as np
import allantools

#Provieds the taumin based on the most frequent
def optimal_tau_min(times):
    diff = np.diff(times)
    hist,border=np.histogram(diff,bins=len(times))
    tau_opt=border[np.argmax(hist)]
    return tau_opt

def IrregAllanPrep(times, values, taumin=-1.):
    diff = np.diff(times)
    if taumin <= 0.:
        taumin = float("%.2e" % (np.min(diff)))  # Find a better way to round
    else:
        taumin = float(taumin)
    splitlist = [[] for i in xrange(int(np.ceil(times[-1] / taumin)) + 1)]
    for i, stamp in enumerate(times):
        splitlist[1+int(np.floor((stamp-taumin/2.) / taumin))].insert(1, [times[i], values[i]])
    return np.array(splitlist), taumin


def singleAllan(splitlist=[], tau0=-1., m=1, maxpnts=np.inf):
    if m <= 0:
        return None
    squarediffs = []
    weights = []
    number_of_bins=len(splitlist)
    #Setup of start indeces
    if maxpnts >= number_of_bins:
        start_indices=range(number_of_bins-2*m)
    else:
        start_indices=range(0,number_of_bins-2*m,(number_of_bins-2*m)/maxpnts)
    for start_index in start_indices:
        #This includes one elem of squarediffs etc.
        blockdiffs = []
        blockmeans = []
        for n in xrange(start_index, number_of_bins - m - start_index, m):  # List of all start points from l till end in m-sized steps
            times = []
            values = []
            for i in splitlist[n:n + m]:  # collecting of subelemtns of one m-sized block
                if i != []:
                    values.append(i[0][1])
                    times.append(i[0][0])
            times.append((n + m) * tau0)
            times.insert(0, n * tau0)
            blockdiffs.append(
                (np.sum([i ** 2 for i in np.diff(np.array(times))]) / (m * tau0)) ** 0.5)  # Get better idea of penalty
            blockmeans.append(np.nanmean(np.array(values)))
        # print splitlist
        # print blockdiffs
        # print blockmeans
        if len(blockdiffs) >= 2:
            for k in range(len(blockdiffs) - 1):
                if not np.isnan(blockmeans[k:k + 2]).any():
                    squarediffs.append((blockmeans[k] - blockmeans[k + 1]) ** 2)
                    weights.append(1. / np.max([blockdiffs[k], blockdiffs[k + 1]]))
    # print squarediffs
    # print weights
    wSum = float(np.sum(weights))
    oadev = np.sqrt(0.5 * np.sum([squarediffs[i] * w for i, w in enumerate(weights)]) / wSum)
    eadv = oadev / np.sqrt(wSum)
    return oadev, eadv, tau0 * m


m = 1
taus = [m]
nop = 5
data = [ i for i in range(nop)]

times=np.array([i for i in xrange(nop)])
rate = (len(times) - 1) / float(times[-1] - times[0])
print "rate %f" % (rate)
def altertheo(m):
    return m * float(1) / np.sqrt(2.)


print "allan theo: %f" % (altertheo(m))
t2, ad, ade, adn = allantools.oadev(data, rate=rate, data_type='freq', taus=taus)
print t2
print "official: %f" % (ad)
print "error official %f" %(ade)
quick = np.sqrt(0.5 * np.mean([(np.mean(data[i:i + m]) - np.mean(data[i + m:i + 2 * m])) ** 2
                               for i in range(0, len(data) - 2 * m)]))
print "quickvalue: %f" % (quick)


splitlist, tauminout = IrregAllanPrep(times, data, taumin=1.)
oadev, eadev, taueff = singleAllan(splitlist, tauminout, m=m, maxpnts=np.inf)
print "error irreg: %f" %(eadev)

print "Irreg Oadev: %f" % (oadev)

