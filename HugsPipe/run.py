from __future__ import division, print_function

import os
import numpy as np
import lsst.pipe.base
from . import utils
from . import primitives as prim

__all__ = ['run']


def run(cfg, debug_return=False):
    """
    Run hugs pipeline.

    Parameters
    ----------
    cfg : HugsPipe.Config 
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
    
    mi_smooth = utils.smooth_gauss(cfg.mi, cfg.psf_sigma)
    fp_low = prim.image_threshold(mi_smooth, mask=cfg.mask, 
                                  plane_name='THRESH_LOW', **cfg.thresh_low)
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

    assoc = prim.associate(cfg.mask, fp_low, **cfg.assoc)
        
    exp_clean = cfg.exp.clone()
    mi_clean = exp_clean.getMaskedImage()
    mi_clean.getImage().getArray()[assoc!=0] = noise_array[assoc!=0]

    ############################################################
    # Smooth with large kernel for detection.
    ############################################################

    kern_fwhm = cfg.thresh_det.pop('kern_fwhm')
    fwhm = kern_fwhm/utils.pixscale # pixels
    sigma = fwhm/(2*np.sqrt(2*np.log(2)))
    mi_clean_smooth = utils.smooth_gauss(mi_clean, sigma)

    ############################################################
    # Image thresholding at final detection threshold 
    ############################################################

    fp_det = prim.image_threshold(mi_clean_smooth, plane_name='DETECTED',
                                  mask=cfg.mask, **cfg.thresh_det)

    ############################################################
    # Deblend sources in 'detected' footprints
    ############################################################

    sources = prim.deblend_stamps(cfg.exp, **cfg.deblend_stamps)
        
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
