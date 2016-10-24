"""
Tests for hugs_pipe.primitives
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
sources = prim.deblend_stamps(exposure)


def test_associate():
    fpset = prim.image_threshold(masked_image, 30, 
                                 plane_name='THRESH_HIGH', 
                                 mask=mask)
    assoc = prim.associate(mask, fpset)
    assert assoc.shape == mask.getArray().shape
    assert assoc.sum() > 0


def test_deblend_stamps():
    assert len(sources) > 0
    assert sources['ra'][0] is not None 
    assert sources['dec'][0] is not None 


def test_image_threshold():
    fpset_1 = prim.image_threshold(masked_image, 20)
    assert len(fpset_1.getFootprints()) > 0
    fpset_2 = prim.image_threshold(masked_image, 20, npix=50)
    assert len(fpset_2.getFootprints()) > 0
    assert len(fpset_2.getFootprints()) < len(fpset_1.getFootprints())


def test_photometry():
    img = exposure.getMaskedImage().getImage().getArray()
    prim.photometry(img, sources)
    assert 'mag_ell' in sources.colnames
    assert 'mag_circ_3' in sources.colnames
