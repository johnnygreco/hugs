"""
Unit tests for hugs.imtools
"""
from __future__ import division, print_function

from ...utils import get_test_exp
from ..core import get_cutout

exposure = get_test_exp()

def test_get_cutout():
    center = (13687, 28321)
    size = 300
    cutout = get_cutout(center, size, exp=exposure)
    img = cutout.getMaskedImage().getImage().getArray()
    assert img.shape == (435, 419)
