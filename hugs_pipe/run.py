from __future__ import division, print_function

import os
import numpy as np
import lsst.pipe.base
from . import utils
from . import imtools
from . import primitives as prim

__all__ = ['run']


def run(cfg, debug_return=False):
    """
    Run hugs pipeline.

    Parameters
    ----------
    cfg : hugs_pipe.Config 
        Configuration object which stores all params 
        as well as the exposure object. 
    debug_return : bool, optional
        If True, return struct with outputs from 
        every step of pipeline.

    Returns
    -------
    sources : astropy.table.Table 
        Source catalog.
    """

    ############################################################
    # Image thesholding at low and high thresholds. In both 
    # cases, the image is smoothed at the psf scale.
    ############################################################
    
    mi_smooth = imtools.smooth_gauss(cfg.mi, cfg.psf_sigma)
    cfg.logger.info('performing low threshold at '
                    '{} sigma'.format(cfg.thresh_low['thresh']))
    fp_low = prim.image_threshold(mi_smooth, mask=cfg.mask, 
                                  plane_name='THRESH_LOW', **cfg.thresh_low)
    cfg.logger.info('performing high threshold at '
                    '{} sigma'.format(cfg.thresh_high['thresh']))
    fp_high = prim.image_threshold(mi_smooth, mask=cfg.mask, 
                                   plane_name='THRESH_HIGH', **cfg.thresh_high)

    ############################################################
    # Generate noise array, which we will use to replace 
    # unwanted sources with noise. 
    ############################################################

    shape = cfg.mask.getArray().shape
    back_rms = cfg.mi.getImage().getArray()[cfg.mask.getArray()==0].std()
    noise_array = back_rms*np.random.randn(shape[0], shape[1])

    ############################################################
    # Get "association" segmentation image; will be non-zero for 
    # low-thresh footprints that are associated with high-thresh
    # footprints. Then, replace these sources with noise.
    ############################################################

    cfg.logger.info('associating low & high thresh objects')
    assoc = prim.associate(cfg.mask, fp_low, **cfg.assoc)
        
    exp_clean = cfg.exp.clone()
    mi_clean = exp_clean.getMaskedImage()
    mi_clean.getImage().getArray()[assoc!=0] = noise_array[assoc!=0]
    mask_clean = mi_clean.getMask()

    ############################################################
    # Smooth with large kernel for detection.
    ############################################################

    kern_fwhm = cfg.thresh_det.pop('kern_fwhm')
    cfg.logger.info('gaussian smooth for detection with '
                    'with fhwm = {} arcsec'.format(kern_fwhm))
    fwhm = kern_fwhm/utils.pixscale # pixels
    sigma = fwhm/(2*np.sqrt(2*np.log(2)))
    mi_clean_smooth = imtools.smooth_gauss(mi_clean, sigma)

    ############################################################
    # Image thresholding at final detection threshold 
    ############################################################

    cfg.logger.info('performing detection threshold at '
                    '{} sigma'.format(cfg.thresh_det['thresh']))
    fp_det = prim.image_threshold(mi_clean_smooth, plane_name='DETECTED',
                                  mask=mask_clean, **cfg.thresh_det)
    fp_det.setMask(cfg.mask, 'DETECTED')

    ############################################################
    # Deblend sources in 'detected' footprints
    ############################################################

    cfg.logger.info('building source catalog')
    sources = prim.deblend_stamps(exp_clean, **cfg.deblend_stamps)
    img = cfg.mi.getImage().getArray()

    cfg.logger.info('measuring aperture magnitudes')
    prim.photometry(img, sources, **cfg.photometry)

    if type(cfg.data_id)==str:
        utils.add_cat_params(sources)
    else:
        tract, patch = cfg.data_id['tract'], cfg.data_id['patch']
        utils.add_cat_params(sources, tract, patch)

    cfg.logger.info('task complete')
        
    if debug_return:
        return lsst.pipe.base.Struct(sources=sources,
                                     exposure=cfg.exp,
                                     exp_clean=exp_clean,
                                     mi_clean_smooth=mi_clean_smooth,
                                     fp_low=fp_low,
                                     fp_high=fp_high,
                                     fp_det=fp_det)
    else:
        return sources
