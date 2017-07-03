"""
Unit tests for hugs.imtools
"""
from __future__ import division, print_function

import numpy as np
from ...utils import get_test_exp
from ..core import *

exposure = get_test_exp()

def test_get_cutout():
    center = (13687, 28321)
    size = 300
    cutout = get_cutout(center, size, exp=exposure)
    img = cutout.getMaskedImage().getImage().getArray()
    assert img.shape == (435, 419)


def test_smooth_gauss():
    mi = exposure.getMaskedImage()
    mi_smooth = smooth_gauss(mi, 10.0)
    arr = mi.getImage().getArray()
    arr_smooth = mi_smooth.getImage().getArray()
    assert np.nanstd(arr) > np.nanstd(arr_smooth)


def test_rmedian():
    data = exposure.getMaskedImage().getImage().getArray()
    filtered_data = rmedian(exposure, 4, 6)
    assert type(filtered_data) == np.ndarray
    assert filtered_data.std() < data.std()


def test_unsharp_mask():
    data = exposure.getMaskedImage().getImage().getArray()
    unsharp = unsharp_mask(exposure, 15)
    assert type(unsharp) == np.ndarray
    assert unsharp.std() < data.std()
