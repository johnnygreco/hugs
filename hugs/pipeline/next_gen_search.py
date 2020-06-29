from __future__ import division, print_function

import os
import numpy as np
from astropy.table import Table, hstack
from lsst.pipe.base import Struct
from ..sep_stepper import SepLsstStepper, sep_ellipse_mask
from ..stats import get_clipped_sig_task
from ..utils import pixscale, zpt
from .. import utils
from .. import imtools
from .. import primitives as prim
from ..cattools import xmatch

__all__ = ['run']

def run(cfg, reset_mask_planes=False):
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
    results : lsst.pipe.base.Struct
        Object containing results:
        results.all_detections : catalog of all detections
        results.sources : catalog of sources we are keeping
        results.exp : exposure object for this run
        results.exp_clean : cleaned exposure object for this run
        results.success : boolean flag of run status 
    """

    assert cfg.tract and cfg.patch, 'No patch id given!'
    cfg.timer # start timer

    ############################################################
    # Get masked image and check if we have enough good data
    ############################################################

    try:

        mi = cfg.exp[cfg.band_detect].getMaskedImage()
        mask = mi.getMask()
        if cfg.band_mask != cfg.band_detect:
            cfg.logger.info('making mask using {}-band'.format(cfg.band_mask))
            mi_band_mask = cfg.exp[cfg.band_mask].getMaskedImage().clone()
        else:
            mi_band_mask = mi
        stat_task = get_clipped_sig_task()

        cfg.logger.info('good data fraction = {:.2f}'.\
            format(cfg.exp.patch_meta.good_data_frac))

        if cfg.exp.patch_meta.good_data_frac < cfg.min_good_data_frac:
            msg = '***** not enough data in {} {} {}-band!!! ****'
            cfg.logger.warning(msg.format(cfg.tract, cfg.patch, cfg.band_mask))
            results = _null_return(cfg)
            return results

        if cfg.band_mask != cfg.band_detect :
            min_good = cfg.min_good_data_frac
            if cfg.exp.good_data_fraction(cfg.band_mask) < min_good:
                msg = '***** not enough data in {} {} {}-band!!! ****'
                cfg.logger.warning(
                    msg.format(cfg.tract, cfg.patch, cfg.band_mask))
                results = _null_return(cfg)
                return results

        ############################################################
        # Image thesholding at low and high thresholds. In both 
        # cases, the image is smoothed at the psf scale.
        ############################################################
            
        #mi_smooth = imtools.smooth_gauss(mi, cfg.psf_sigma)
        stats = stat_task.run(mi_band_mask)

        if stats.stdev <= 0.0 or np.isnan(stats.stdev):
            msg = '***** {} | {} -- stddev = {} !!! ****'
            cfg.logger.warning(msg.format(cfg.tract, cfg.patch, stats.stdev))
            results = _null_return(cfg)
            return results

        if cfg.thresh_type.lower() == 'sb':
            cfg.logger.info('thresh type set to ' + cfg.thresh_type)
            low_th = cfg.thresh_low['thresh'] 
            flux_th = 10**(0.4 * (zpt - low_th)) * pixscale**2
            cfg.thresh_low['thresh'] = flux_th / stats.stdev    
            high_th = cfg.thresh_high['thresh']
            flux_th = 10**(0.4 * (zpt - high_th)) * pixscale**2
            cfg.thresh_high['thresh'] = flux_th / stats.stdev    
        elif cfg.thresh_type.lower() == 'stddev':
            cfg.logger.info('thresh type set to ' + cfg.thresh_type)
        else:
            raise Exception('invalid threshold type')

        cfg.logger.info('performing low threshold at '
                        '{:.2f} sigma'.format(cfg.thresh_low['thresh']))
        fpset_low = prim.image_threshold(
            mi_band_mask, mask=mask, plane_name='THRESH_LOW', **cfg.thresh_low)
        cfg.logger.info('performing high threshold at '
                        '{:.2f} sigma'.format(cfg.thresh_high['thresh']))
        fpset_high = prim.image_threshold(
            mi_band_mask, mask=mask, plane_name='THRESH_HIGH', 
            **cfg.thresh_high)

        ############################################################
        # Get "cleaned" image, with noise replacement
        ############################################################

        cfg.logger.info('generating cleaned exposure')
        exp_clean = prim.clean(cfg.exp[cfg.band_detect], 
                               fpset_low, 
                               **cfg.clean)

        ############################################################
        # Remove small sources using HSC pipeline detection map
        ############################################################

        if cfg.hsc_small_sources_r_max is not None:
            if cfg.hsc_small_sources_r_max > 0:
                cfg.logger.info('removed small sources with HSC detection map')
                exp_clean = prim.remove_small_sources_thresholding(
                    exp_clean, cfg.hsc_small_sources_r_max, cfg.rng)
        
        mi_clean = exp_clean.getMaskedImage()
        mask_clean = mi_clean.getMask()

        ############################################################
        # use sep to find and mask point-like sources
        ############################################################
        
        if cfg.sep_steps is not None:
            sep_stepper = SepLsstStepper(config=cfg.sep_steps)
            sep_stepper.setup_image(exp_clean, cfg.rng)

            step_mask = cfg.exp.get_mask_array(
                planes=['BRIGHT_OBJECT', 'NO_DATA', 'SAT'])
            sep_sources, _ = sep_stepper.run('sep_point_sources', 
                                             mask=step_mask)

            cfg.logger.info('generating and applying sep ellipse mask')
            r_min = cfg.sep_min_radius
            sep_sources = sep_sources[sep_sources['flux_radius'] < r_min]
            ell_msk = sep_ellipse_mask(
                sep_sources, sep_stepper.image.shape, cfg.sep_mask_grow)
            nimage_replace = sep_stepper.noise_image[ell_msk]
            mi_clean.getImage().getArray()[ell_msk] = nimage_replace
            mask_clean.addMaskPlane('SMALL')
            mask_clean.getArray()[ell_msk] += mask_clean.getPlaneBitMask('SMALL')

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
                dual_exp=dual_exp, 
                delete_created_files=cfg.delete_created_files, 
                original_fn=cfg.exp.fn[cfg.band_detect]) 
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
                label=label, delete_created_files=cfg.delete_created_files, 
                original_fn=cfg.exp.fn[band])
            if len(sources_verify)>0:
                match_masks, _ = xmatch(
                    sources, sources_verify, max_sep=cfg.verify_max_sep)
                txt = 'cuts: {} out of {} objects detected in {}-band'.format(
                    len(match_masks[0]), len(sources), band)
                cfg.logger.info(txt)
                if len(match_masks[0])==0:
                    cfg.logger.warn(
                        '**** no matched sources with '+band+' ****')
                    results = _null_return(cfg, exp_clean)
                    return results
                sources = sources[match_masks[0]]
            else:
                cfg.logger.warn('**** no sources detected in '+band+' ****')
                results = _null_return(cfg, exp_clean)
                return results

        mask_fracs = utils.calc_mask_bit_fracs(exp_clean)
        cfg.exp.patch_meta.small_frac = mask_fracs['small_frac']
        cfg.exp.patch_meta.cleaned_frac = mask_fracs['cleaned_frac']
        cfg.exp.patch_meta.bright_obj_frac = mask_fracs['bright_object_frac']

        cfg.logger.info('measuring mophology metrics')
        prim.measure_morphology_metrics(exp_clean.getImage().getArray(), 
                                        sources)

        cfg.logger.info('task completed in {:.2f} min'.format(cfg.timer))

        results = Struct(all_detections=all_detections,
                         sources=sources,
                         hugs_exp=cfg.exp,
                         exp_clean=exp_clean,
                         success=True, 
                         synths=cfg.exp.synths)

        if reset_mask_planes:
            cfg.reset_mask_planes()

        return results

    except:
        cfg.logger.critical(
            'tract - patch {} - {} failed'.format(cfg.tract, cfg.patch))
        results = _null_return(cfg)
        return results


def _null_return(config, exp_clean=None):
    config.reset_mask_planes()
    return Struct(all_detections=None,
                  sources=None,
                  hugs_exp=config.exp,
                  exp_clean=exp_clean,
                  success=False, 
                  synths=None)
