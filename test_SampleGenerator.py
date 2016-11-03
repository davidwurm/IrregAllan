import pytest
import numpy as np
from SampleGenerator import Drift_Generator


def test_Drift_Generator():
    allan_Drift_theo, points_Drift_list = Drift_Generator(taus=[2,3,4], slope=1.,points=7,fs=1.0)
    assert (points_Drift_list == np.array([0,1,2,3,4,5,6])).all
