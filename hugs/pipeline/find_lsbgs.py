from __future__ import division, print_function

import os
import numpy as np
import lsst.pipe.base
from astropy.table import Table, hstack
from .. import utils
from .. import imtools
from .. import primitives as prim
from .. import randoms 
from ..cattools import xmatch

__all__ = ['run']

def run(cfg, randoms_only=False):
    """
    Run hugs pipeline using SExtractor for the final detection 
    and photometry.

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

    if cfg.exp.good_data_fraction() < 0.2:
        cfg.logger.warning('***** not enough data!!! ****')
        results = _null_return(cfg.exp)
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
    # Find randoms in footprint that are not masked
    ############################################################

    randoms_results = None
    if cfg.randoms_db_fn is not None:
        cfg.logger.info('finding detectable randoms in patch')
        randoms_df, randoms_db = randoms.find_randoms_in_footprint(
            cfg.randoms_db_fn, exp_clean, return_db=True)
        randoms_results = lsst.pipe.base.Struct(df=randoms_df, db=randoms_db)

    if randoms_only:
        cfg.logger.info('task completed in {:.2f} min'.format(cfg.timer))
        results = _null_return(cfg.exp, exp_clean, randoms_results)
        return results

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
            exp_clean, cfg.sex_config, cfg.sex_params, cfg.sex_io_dir, 
            label=label, dual_exp=dual_exp, 
            delete_created_files=cfg.delete_created_files) 
        if len(sources_band)>0:
            sources = hstack([sources, sources_band])
        else:
            cfg.logger.warn('**** no sources found by sextractor ****')
            results = _null_return(cfg.exp, exp_clean, randoms_results)
            return results

    ############################################################
    # Verify detections in other bands using SExtractor
    ############################################################

    all_detections = sources.copy()

    for band in cfg.band_verify:
        cfg.logger.info('verifying dection in {}-band'.format(band))
        sources_verify = prim.detect_sources(
            cfg.exp[band], cfg.sex_config, cfg.sex_params, cfg.sex_io_dir,
            label=label, delete_created_files=cfg.delete_created_files)
        if len(sources_verify)>0:
            match_masks, _ = xmatch(
                sources, sources_verify, max_sep=cfg.verify_max_sep)
            txt = 'cuts: {} out of {} objects detected in {}-band'.format(
                len(match_masks[0]), len(sources), band)
            cfg.logger.info(txt)
            if len(match_masks[0])==0:
                cfg.logger.warn('**** no matched sources with '+band+' ****')
                results = _null_return(cfg.exp, exp_clean, randoms_results)
                return results
            sources = sources[match_masks[0]]
        else:
            cfg.logger.warn('**** no sources detected in '+band+' ****')
            results = _null_return(cfg.exp, exp_clean, randoms_results)
            return results

    mask_fracs = utils.calc_mask_bit_fracs(exp_clean)
    cfg.logger.info('task completed in {:.2f} min'.format(cfg.timer))
    results = lsst.pipe.base.Struct(all_detections=all_detections,
                                    sources=sources,
                                    exp=cfg.exp,
                                    exp_clean=exp_clean,
                                    mask_fracs=mask_fracs, 
                                    randoms_results=randoms_results)
    return results


def _null_return(exp, exp_clean=None, randoms_results=None):
    if exp_clean is not None:
        mask_fracs = utils.calc_mask_bit_fracs(exp_clean)
    else:
        mask_fracs = {}
    return lsst.pipe.base.Struct(all_detections=Table(),
                                 sources=Table(), 
                                 candy=Table(),
                                 mask_fracs=mask_fracs,
                                 randoms_results=randoms_results,
                                 exp=exp,
                                 exp_clean=exp_clean)