"""
Tests for hugsPipe.primitives
"""
from __future__ import division, print_function

import os
import lsst.afw.image
from .. import primitives as prim
from .. import utils

dataDIR = os.environ.get('TEST_DATA_DIR')
fn = os.path.join(dataDIR, 'test_exposure.fits')
exposure = lsst.afw.image.ExposureF(fn)
masked_image = exposure.getMaskedImage().clone()
mask = masked_image.getMask().clone()


def test_image_threshold():
    fpset_1 = prim.image_threshold(masked_image, 20)
    assert len(fpset_1.getFootprints()) > 0
    fpset_2 = prim.image_threshold(masked_image, 20, npix=50)
    assert len(fpset_2.getFootprints()) > 0
    assert len(fpset_2.getFootprints()) < len(fpset_1.getFootprints())


def test_associate():
    fpset = prim.image_threshold(masked_image, 30, 
                                 plane_name='THRESH_HIGH', 
                                 mask=mask)
    assoc = prim.associate(mask, fpset)
    assert assoc.shape == mask.getArray().shape
    assert assoc.sum() > 0


def test_deblend_stamps():
    exp_wcs = utils.get_astropy_wcs(fn)
    table = prim.deblend_stamps(exposure, wcs=exp_wcs)
    assert len(table) > 0
    assert table['ra_icrs_centroid'][0] is not None 
    assert table['dec_icrs_centroid'][0] is not None 
