import pytest
import numpy as np
import Irreg_Core

def test_evenness():
    times=np.array([i for i in xrange(10)])
    values=np.array([1. for i in xrange(10)])
    splitlist,taumin = Irreg_Core.IrregAllanPrep(times,values,taumin=1.0)
    out=Irreg_Core.optimal_tau_min(times)
    assert splitlist.shape==(10,1,2)

