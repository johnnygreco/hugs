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

def run(cfg, debug_return=False, synth_factory=None, randoms_only=False):
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
    synth_factory : hugs_pipe.SynthFactory
        Synth injecting object. If None, no synths 
        will be injected.

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

    nodata = mask.getArray() & mask.getPlaneBitMask('NO_DATA') != 0
    nodata = nodata.astype(float)
    no_data_frac = nodata.sum()/nodata.size
    if no_data_frac > 0.8:
        cfg.logger.warning('***** not enough data!!! ****')
        _exp = cfg.exp[cfg.band_detect]
        debug_exp = {'exp':cfg.exp, 'exp_clean':_exp} if debug_return else {}
        results = _null_return(None, debug_exp)
        cfg.reset_mask_planes()
        return results

    ############################################################
    # If desired, inject synthetic galaxies 
    ############################################################
    
    if synth_factory:
        num_synths = synth_factory.num_synths
        cfg.logger.warning('**** injecting {} synths ****'.format(num_synths))
        for band in cfg.bands:
            synth_factory.inject(cfg.exp[band], band=band)

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

    del mi_smooth

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

    if cfg.randoms_db_fn is not None:
        cfg.logger.info('finding detectable randoms in patch')
        randoms_df, randoms_db = randoms.find_randoms_in_footprint(
            cfg.randoms_db_fn, exp_clean, return_db=True)
        randoms_results = lsst.pipe.base.Struct(df=randoms_df, db=randoms_db)
    else:
        randoms_results = None

    if randoms_only:
        debug_exp = {'exp':cfg.exp, 'exp_clean':exp_clean} if debug_return else {}
        results = _null_return(randoms_results, debug_exp)
        if not debug_return:
            cfg.reset_mask_planes()
        cfg.logger.info('task completed in {:.2f} min'.format(cfg.timer))
        return results

    ############################################################
    # Detect sources and measure props with SExtractor
    ############################################################

    cfg.logger.info('detecting in {}-band'.format(cfg.band_detect))
    label = '{}-{}-{}'.format(cfg.tract, cfg.patch[0], cfg.patch[-1])
    if cfg.group_id:
        label = str(cfg.group_id)+'-'+label

    cfg.sex_measure['sf'] = synth_factory
    num_aps = len(cfg.sex_measure['apertures'])

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
        add_params = band==cfg.band_detect
        dual_exp = None if band==cfg.band_detect else cfg.exp[band]
        sources_band = prim.sex_measure(
            exp_clean, label=label+'-'+cfg.band_detect, add_params=add_params, 
            dual_exp=dual_exp, **cfg.sex_measure)
        if len(sources_band)>0:
            utils.add_band_to_name(sources_band, band, num_aps)
            if band!=cfg.band_detect:
                cols = [c for c in sources_band.colnames if '('+band+')' in c]
                sources = hstack([sources, sources_band[cols]])
            else:
                sources = hstack([sources, sources_band])
        else:
            break

    if len(sources)>0:

        ############################################################
        # Verify detections in other bands using SExtractor
        ############################################################

        all_detections = sources.copy()

        for band in cfg.band_verify:
            cfg.logger.info('verifying dection in {}-band'.format(band))
            sources_verify = prim.sex_measure(
                cfg.exp[band], label=label+'-'+band, 
                add_params=False, **cfg.sex_measure)
            if len(sources_verify)>0:
                match_masks, _ = xmatch(
                    sources, sources_verify, max_sep=cfg.verify_max_sep)
                txt = 'cuts: {} out of {} objects detected in {}-band'.format(
                    len(match_masks[0]), len(sources), band)
                cfg.logger.info(txt)
                if len(match_masks[0])==0:
                    sources = Table()
                    break
                sources = sources[match_masks[0]]
            else:
                sources = Table()
                break

    else: 
        cfg.logger.warn('**** no sources found by sextractor ****')
        all_detections = Table()
        candy = Table()

    mask_fracs = utils.calc_mask_bit_fracs(exp_clean)

    cfg.logger.info('task completed in {:.2f} min'.format(cfg.timer))

    if debug_return:
        results = lsst.pipe.base.Struct(all_detections=all_detections,
                                        sources=sources,
                                        exp=cfg.exp,
                                        exp_clean=exp_clean,
                                        mask_fracs=mask_fracs, 
                                        randoms_results=randoms_results)
    else:
        results = lsst.pipe.base.Struct(all_detections=all_detections,
                                        sources=sources, 
                                        mask_fracs=mask_fracs,
                                        randoms_results=randoms_results)
        cfg.reset_mask_planes()

    return results


def _null_return(randoms_results, debug_exp={}):
    return lsst.pipe.base.Struct(all_detections=Table(),
                                 sources=Table(), 
                                 candy=Table(),
                                 mask_fracs={},
                                 randoms_results=randoms_results,
                                 **debug_exp)
