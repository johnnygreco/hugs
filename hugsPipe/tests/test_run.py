"""
Test hugsPipe.run
"""

import os
import pytest
import lsst.afw.image
from ..run import run

dataDIR = os.environ.get('TEST_DATA_DIR')
fn = os.path.join(dataDIR, 'test_exposure.fits')
exposure = lsst.afw.image.ExposureF(fn)


def test_run():
    exp = run(
        fn, thresh={'high':20.0, 'det':4.0, 'low':3.0},
        npix={'det':50}, assoc={'r_in':0, 'min_pix': '6 psf sigma'},
        kern_fwhm=4.0, visualize=False)
    mask = exp.getMaskedImage().getMask()
    planes = list(mask.getMaskPlaneDict().keys())
    assert 'DETECTED' in planes
    assert 'THRESH_LOW' in planes
    assert 'THRESH_HIGH' in planes
