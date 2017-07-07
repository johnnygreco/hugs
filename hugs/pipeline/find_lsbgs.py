from __future__ import division, print_function

import os
import numpy as np
from astropy.table import Table, hstack
from lsst.pipe.base import Struct
from .. import utils
from .. import imtools
from .. import primitives as prim
from ..cattools import xmatch

__all__ = ['run']

def run(cfg):
    """
    Run hugs pipeline using SExtractor for the final detection 
    and photometry.

    Parameters
    ----------
    cfg : hugs_pipe.Config 
        Configuration object which stores all params 
        as well as the exposure object. 

    Returns
    -------
    results : astropy.table.Table 
        Source catalog.
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
        results = _null_return(cfg)
        return results

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
    mi_clean = exp_clean.getMaskedImage()
    mask_clean = mi_clean.getMask()

    ############################################################
    # Detect sources and measure props with SExtractor
    ############################################################

    cfg.logger.info('detecting in {}-band'.format(cfg.band_detect))
    label = '{}-{}-{}'.format(cfg.tract, cfg.patch[0], cfg.patch[-1])

    cfg.logger.info('cleaning non-detection bands')
    replace = cfg.exp.get_mask_array(cfg.band_detect)
    for band in cfg.bands:
        if band!=cfg.band_detect:
            mi_band = cfg.exp[band].getMaskedImage()
            noise_array = utils.make_noise_image(mi_band, cfg.rng)
            mi_band.getImage().getArray()[replace] = noise_array[replace]

    sources = Table()

    for band in cfg.bands:
        cfg.logger.info('measuring in {}-band'.format(band))
        dual_exp = None if band==cfg.band_detect else cfg.exp[band]
        sources_band = prim.detect_sources(
            exp_clean, cfg.sex_config, cfg.sex_io_dir, label=label, 
            dual_exp=dual_exp, delete_created_files=cfg.delete_created_files) 
        if len(sources_band)>0:
            sources = hstack([sources, sources_band])
        else:
            cfg.logger.warn('**** no sources found by sextractor ****')
            results = _null_return(cfg, exp_clean)
            return results

    ############################################################
    # Verify detections in other bands using SExtractor
    ############################################################

    all_detections = sources.copy()

    for band in cfg.band_verify:
        cfg.logger.info('verifying dection in {}-band'.format(band))
        sources_verify = prim.detect_sources(
            cfg.exp[band], cfg.sex_config, cfg.sex_io_dir,
            label=label, delete_created_files=cfg.delete_created_files)
        if len(sources_verify)>0:
            match_masks, _ = xmatch(
                sources, sources_verify, max_sep=cfg.verify_max_sep)
            txt = 'cuts: {} out of {} objects detected in {}-band'.format(
                len(match_masks[0]), len(sources), band)
            cfg.logger.info(txt)
            if len(match_masks[0])==0:
                cfg.logger.warn('**** no matched sources with '+band+' ****')
                results = _null_return(cfg, exp_clean)
                return results
            sources = sources[match_masks[0]]
        else:
            cfg.logger.warn('**** no sources detected in '+band+' ****')
            results = _null_return(cfg, exp_clean)
            return results

    mask_fracs = utils.calc_mask_bit_fracs(exp_clean)
    cfg.exp.patch_meta.cleaned_frac = mask_fracs['cleaned_frac']
    cfg.exp.patch_meta.bright_obj_frac = mask_fracs['bright_object_frac']

    cfg.logger.info('task completed in {:.2f} min'.format(cfg.timer))
    results = Struct(all_detections=all_detections,
                     sources=sources,
                     exp=cfg.exp,
                     exp_clean=exp_clean,
                     success=True)

    cfg.reset_mask_planes()
    return results


def _null_return(config, exp_clean=None):
    config.reset_mask_planes()
    return Struct(all_detections=None,
                  sources=None,
                  exp=config.exp,
                  exp_clean=exp_clean,
                  success=False)
