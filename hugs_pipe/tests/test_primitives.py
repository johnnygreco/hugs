"""
Unit tests for hugs_pipe.primitives
"""
from __future__ import division, print_function

from .. import primitives as prim
from .. import utils
from .. import imtools

exposure = utils.get_test_exp()
masked_image = exposure.getMaskedImage()
mask = masked_image.getMask()
utils.remove_mask_planes(mask, ['CR', 'CROSSTALK', 'DETECTED_NEGATIVE'])
mask.addMaskPlane('BLEND')
sources = prim.measure_sources(exposure)

def test_clean():
    utils.remove_mask_planes(mask, ['BLEND'])
    fpset_hi = prim.image_threshold(masked_image, 30, 
                                 plane_name='THRESH_HIGH', 
                                 mask=mask)
    fpset_lo = prim.image_threshold(masked_image, 5, 
                                 plane_name='THRESH_LOW', 
                                 mask=mask)
    exp_clean = prim.clean(exposure, fpset_lo)
    img_clean = exp_clean.getMaskedImage().getImage().getArray()
    img = exposure.getMaskedImage().getImage().getArray()
    assert 'CLEANED' in mask.getMaskPlaneDict().keys()
    assert img.std() > img_clean.std()


def test_find_blends():
    utils.remove_mask_planes(mask, ['CLEANED'])
    mi_smooth = imtools.smooth_gauss(masked_image, 8.0)
    fpset_det = prim.image_threshold(mi_smooth, 5, 
                                     plane_name='DETECTED', 
                                     mask=mask)
    prim.find_blends(exposure, fpset_det)
    assert 'BLEND' in mask.getMaskPlaneDict().keys()


def test_measure_sources():
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
