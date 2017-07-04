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


def test_image_threshold():
    fpset_1 = prim.image_threshold(masked_image, 20)
    assert len(fpset_1.getFootprints()) > 0
    fpset_2 = prim.image_threshold(masked_image, 20, npix=50)
    assert len(fpset_2.getFootprints()) > 0
    assert len(fpset_2.getFootprints()) < len(fpset_1.getFootprints())
