from __future__ import division, print_function

import os
import numpy as np
import lsst.pipe.base
from . import utils
from . import imtools
from . import primitives as prim

__all__ = ['run']


def run(cfg, debug_return=False, synth_factory=None):
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
    synth_factory : hugs_pipe.SynthFactory
        Synth injecting object. If None, no synths 
        will be injected.

    Returns
    -------
    results : astropy.table.Table 
        Source catalog.
    """

    assert cfg.data_id is not None, 'No data id!'
    cfg.timer # start timer

    ############################################################
    # If desired, inject synthetic galaxies 
    ############################################################
    
    if synth_factory:
        num_synths = synth_factory.num_synths
        cfg.logger.warning('**** injecting {} synths ****'.format(num_synths))
        synth_factory.inject(cfg.exp, band='i')
        if cfg.phot_colors:
            for band in cfg.color_data.keys():
                synth_factory.inject(cfg.color_data[band], band=band)

    ############################################################
    # Image thesholding at low and high thresholds. In both 
    # cases, the image is smoothed at the psf scale.
    ############################################################
    
    mi_smooth = imtools.smooth_gauss(cfg.mi, cfg.psf_sigma)
    cfg.logger.info('performing low threshold at '
                    '{} sigma'.format(cfg.thresh_low['thresh']))
    fpset_low = prim.image_threshold(
        mi_smooth, mask=cfg.mask, plane_name='THRESH_LOW', **cfg.thresh_low)
    cfg.logger.info('performing high threshold at '
                    '{} sigma'.format(cfg.thresh_high['thresh']))
    fpset_high = prim.image_threshold(
        mi_smooth, mask=cfg.mask, plane_name='THRESH_HIGH', **cfg.thresh_high)

    ############################################################
    # Get "cleaned" image, with noise replacement
    ############################################################

    cfg.logger.info('generating cleaned exposure')
    exp_clean = prim.clean(cfg.exp, fpset_low, **cfg.clean)
    mi_clean = exp_clean.getMaskedImage()
    mask_clean = mi_clean.getMask()

    ############################################################
    # Smooth with large kernel for detection.
    ############################################################

    kern_fwhm = cfg.thresh_det.pop('kern_fwhm')
    cfg.logger.info('gaussian smooth for detection with '
                    'with fhwm = {} arcsec'.format(kern_fwhm))
    fwhm = kern_fwhm/utils.pixscale # pixels
    sigma = fwhm/(2*np.sqrt(2*np.log(2)))
    mi_clean_smooth = imtools.smooth_gauss(mi_clean, sigma, use_scipy=True)

    ############################################################
    # Image thresholding at final detection threshold 
    ############################################################

    cfg.logger.info('performing detection threshold at '
                    '{} sigma'.format(cfg.thresh_det['thresh']))
    fpset_det = prim.image_threshold(mi_clean_smooth, plane_name='DETECTED',
                                     mask=mask_clean, **cfg.thresh_det)
    fpset_det.setMask(cfg.mask, 'DETECTED')

    ############################################################
    # Find and remove obvious blends and deblend remaining
    # sources in 'detected' footprints
    ############################################################

    cfg.logger.info('finding obvious blends')
    prim.find_blends(exp_clean, fpset_det, **cfg.find_blends)

    cfg.logger.info('building source catalog')
    sources = prim.measure_sources(
        exp_clean, logger=cfg.logger, sf=synth_factory, **cfg.measure_sources)
    img_data = cfg.mi.getImage().getArray()

    ############################################################
    # Perform aperture photometry
    ############################################################

    cfg.logger.info('performing aperture photometry')
    if cfg.phot_colors:
        cfg.color_data['I'] = img_data
        prim.photometry(cfg.color_data, sources, **cfg.photometry)
    else:
        prim.photometry(img_data, sources, **cfg.photometry)

    if type(cfg.data_id)==dict:
        tract, patch = cfg.data_id['tract'], cfg.data_id['patch']
        utils.add_cat_params(sources, tract, patch)
    else:
        utils.add_cat_params(sources)

    mask_fracs = utils.calc_mask_bit_fracs(exp_clean)

    cfg.logger.info('task completed in {:.2f} min'.format(cfg.timer))
        
    if debug_return:
        # detection plane was modified by find_blends
        fpset_det = utils.get_fpset(mask_clean, 'DETECTED')
        results = lsst.pipe.base.Struct(sources=sources,
                                        exposure=cfg.exp,
                                        exp_clean=exp_clean,
                                        mi_clean_smooth=mi_clean_smooth,
                                        fpset_low=fpset_low,
                                        fpset_high=fpset_high,
                                        fpset_det=fpset_det,
                                        mask_fracs=mask_fracs)
    else:
        results = lsst.pipe.base.Struct(sources=sources, 
                                        mask_fracs=mask_fracs)
        cfg.reset_mask_planes()

    return results
