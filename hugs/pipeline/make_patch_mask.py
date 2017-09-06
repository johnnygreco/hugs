from __future__ import division, print_function

import os
import numpy as np
from astropy.table import Table, hstack
from lsst.pipe.base import Struct
from .. import utils
from .. import imtools
from .. import primitives as prim
from ..cattools import xmatch
from ..randoms import get_mask_array
import lsstutils

__all__ = ['run']

def run(cfg):
    """
    Run hugs pipeline to make patch mask for a randoms catalog. 

    Parameters
    ----------
    cfg : hugs_pipe.Config 
        Configuration object which stores all params 
        as well as the exposure object. 

    Returns
    -------
    mask_array : ndarray
        The mask array.
    """

    assert cfg.tract and cfg.patch, 'No patch id given!'
    cfg.timer # start timer

    ############################################################
    # Get masked image and check if we have enough good data
    ############################################################

    mi = cfg.exp[cfg.band_detect].getMaskedImage()
    mask = mi.getMask()

    if cfg.exp.patch_meta.good_data_frac < 0.2:
        cfg.logger.warning('***** not enough data!!! ****')
        return Struct(mask_array=None, ra_lim=None, dec_lim=None)

    ############################################################
    # Image thesholding at low and high thresholds. In both 
    # cases, the image is smoothed at the psf scale.
    ############################################################
        
    mi_smooth = imtools.smooth_gauss(mi, cfg.psf_sigma)
    cfg.logger.info('performing low threshold at '
                    '{} sigma'.format(cfg.thresh_low['thresh']))
    fpset_low = prim.image_threshold(
        mi_smooth, mask=mask, plane_name='THRESH_LOW', **cfg.thresh_low)
    cfg.logger.info('performing high threshold at '
                    '{} sigma'.format(cfg.thresh_high['thresh']))
    fpset_high = prim.image_threshold(
        mi_smooth, mask=mask, plane_name='THRESH_HIGH', **cfg.thresh_high)

    ############################################################
    # Get "cleaned" image, with noise replacement
    ############################################################

    cfg.logger.info('generating cleaned exposure')
    exp_clean = prim.clean(cfg.exp[cfg.band_detect], fpset_low, **cfg.clean)

    mask_array = get_mask_array(exp_clean)
    cfg.reset_mask_planes()

    corners = lsstutils.bbox_to_radec(cfg.exp[cfg.band_detect])
    ra_lim = [corners[:, 0].min(), corners[:, 0].max()]
    dec_lim = [corners[:, 1].min(), corners[:, 1].max()]

    cfg.logger.info('task completed in {:.2f} min'.format(cfg.timer))

    results = Struct(mask_array=mask_array, 
                     ra_lim=ra_lim, 
                     dec_lim=dec_lim)
    
    return results
