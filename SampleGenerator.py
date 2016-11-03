import numpy as np
import allantools


def Drift_Generator(slope=1.,taus=[],points=1000,fs=1.0):
    try:
        taus=np.array(taus)
        slope=float(slope)
        points=int(points)
        fs=float(fs)
    except:
        print "Input parsing failed for Drift Generator"
        raise TypeError
    allanlist = np.array([2. * (tau * slope) ** 2 for tau in taus])
    outpoints = np.array([slope*point/fs for point in xrange(points)])
    return allanlist,outpoints

def White_Generator(h_white=1.,taus=[],points=1000,fs=1.0):
    try:
        taus=np.array(taus)
        h_white=float(h_white)
        points=np.array(points)
    except:
        raise TypeError
    allanlist = np.array([np.sqrt(float(h_white)/(2.*tau)) for tau in taus])
    outpoints=allantools.noise.white(N=points,b0=h_white,fs=fs)
    return allanlist,outpoints

