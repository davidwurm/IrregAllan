import pytest
import numpy as np
from MAIN import AllanI as ai

def test_default_input():
    test_obj=ai()
    test_obj.importData()
    assert test_obj.status=="Init"

def test_std_rate_input():
    test_obj=ai()
    test_obj.importData(rate=23.,data=[3,4,5,6,7,8,6,4,3,4,6,2])
    assert test_obj.status=="Rate Data loaded"

def test_segmented_date_rate_input():
    test_obj=ai()
    test_obj.importData(rate=23.,data=[[3,4,5,6,7,8],[6,4,3,4,6,2]])
    assert test_obj.status=="Init"

def test_tau_list():
    test_obj=ai()
    test_obj.importData(rate=2., data=[3, 4, 5, 6, 7, 8, 9,10,43,43,4,4,43,3,7])
    test_obj.setup_taus(tau_in=[0.0007,3., 5.4,6.,1234.])
    assert (test_obj.tau_list==np.array([3,5,6],dtype=long)).all


def test_tau_default():
    test_obj=ai()
    test_obj.importData(rate=1., data=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    test_obj.setup_taus(no_of_taus=100)
    assert (test_obj.tau_list == np.array([1,2,3,4,5,6], dtype=long)).all

