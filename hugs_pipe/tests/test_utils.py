"""
Unit test for hugs.utils.
"""

import numpy as np
from numpy.testing import assert_array_equal
from ..utils import annuli, get_psf_sigma, get_test_exp

exp = get_test_exp()


def test_annuli():
    arr = np.zeros((10, 10))
    row_idx, col_idx = annuli(6, 6, 0, 5, arr.shape)
    arr[row_idx, col_idx] = 1
    _arr = np.array(
        [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
         [0, 0, 0, 0, 1, 1, 1, 1, 1, 0],
         [0, 0, 0, 1, 1, 1, 1, 1, 1, 1],
         [0, 0, 1, 1, 1, 1, 1, 1, 1, 1],
         [0, 0, 1, 1, 1, 1, 1, 1, 1, 1],
         [0, 0, 1, 1, 1, 1, 1, 1, 1, 1],
         [0, 0, 1, 1, 1, 1, 1, 1, 1, 1],
         [0, 0, 1, 1, 1, 1, 1, 1, 1, 1],
         [0, 0, 0, 1, 1, 1, 1, 1, 1, 1]])
    assert_array_equal(arr, _arr)


def test_get_psf_sigma():
    sigma = get_psf_sigma(exp)
    assert sigma > 0.0
