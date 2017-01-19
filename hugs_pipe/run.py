from __future__ import division, print_function

import os
import numpy as np
import lsst.pipe.base
from astropy.table import Table, hstack
from . import utils
from . import imtools
from . import primitives as prim
from .cattools import cutter, xmatch

__all__ = ['run']

sex_io = os.environ.get('SEX_IO_DIR')

def run(cfg, debug_return=False, synth_factory=None):
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
        
    mi = cfg.exp[cfg.band_detect].getMaskedImage()
    mask = mi.getMask()
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

    if cfg.remove_fpt_blends:

        ############################################################
        # Smooth with large kernel for footprint blend detection
        ############################################################

        kern_fwhm = cfg.thresh_det.pop('kern_fwhm')
        cfg.logger.info('gaussian smooth for blend detection with '
                        'with fhwm = {} arcsec'.format(kern_fwhm))
        fwhm = kern_fwhm/utils.pixscale # pixels
        sigma = fwhm/(2*np.sqrt(2*np.log(2)))
        mi_clean_smooth = imtools.smooth_gauss(mi_clean, sigma, use_scipy=True)

        ############################################################
        # Image thresholding 
        ############################################################

        cfg.logger.info('performing detection threshold at '
                        '{} sigma'.format(cfg.thresh_det['thresh']))
        fpset_smooth = prim.image_threshold(
            mi_clean_smooth, plane_name='SMOOTHED', 
            mask=mask_clean, **cfg.thresh_det)
        fpset_smooth.setMask(mask, 'SMOOTHED')

        ############################################################
        # Find and remove "obvious" blends 
        ############################################################

        cfg.logger.info('finding obvious blends')
        prim.find_blends(exp_clean, fpset_smooth, 
                         plane_name='SMOOTHED', **cfg.find_blends)

        noise_array = utils.make_noise_image(mi_clean, cfg.rng)
        replace = mask_clean.getArray() &\
                  mask_clean.getPlaneBitMask('BLEND') != 0
        mi_clean.getImage().getArray()[replace] = noise_array[replace]

    ############################################################
    # Detect sources and measure props with SExtractor
    ############################################################

    cfg.logger.info('detecting in {}-band'.format(cfg.band_detect))
    label = '{}-{}-{}'.format(cfg.tract, cfg.patch[0], cfg.patch[-1])
    sources = prim.sex_measure(
        exp_clean, label=label+'-'+cfg.band_detect, 
        add_params=True, **cfg.sex_measure)
    num_aps = len(cfg.sex_measure['apertures'])
    utils.add_band_to_name(sources, cfg.band_detect, num_aps)
    all_detections = sources.copy()

    ############################################################
    # Verify detections in other bands using SExtractor
    ############################################################

    replace = cfg.exp.get_mask_array(cfg.band_detect)
    for band in cfg.band_verify:
        mi_band = cfg.exp[band].getMaskedImage()
        cfg.logger.info('verifying dection in {}-band'.format(band))
        noise_array = utils.make_noise_image(mi_band, cfg.rng)
        mi_band.getImage().getArray()[replace] = noise_array[replace]
        sources_verify = prim.sex_measure(
            cfg.exp[band], label=label+'-'+band, 
            add_params=False, **cfg.sex_measure)
        match_masks, _ = xmatch(
            sources, sources_verify, max_sep=cfg.verify_max_sep)
        sources = sources[match_masks[0]]
        sources_verify = sources_verify[match_masks[1]]
        utils.add_band_to_name(sources_verify, band, num_aps)
        cols = [c for c in sources_verify.colnames if '('+band+')' in c]
        sources = hstack([sources, sources_verify[cols]])

    ############################################################
    # Cut catalog and run imfit on cutouts
    ############################################################

    if len(sources)>0:
        cfg.logger.info('cutting catalog')
        candy = cutter(sources.to_pandas(), 
                       min_cuts=cfg.min_cuts, 
                       max_cuts=cfg.max_cuts)
        candy = Table.from_pandas(candy)
        if len(candy)>0:
            cfg.logger.info(
                'running imfit on '+cfg.band_detect+'-band cutouts')
            candy = prim.run_imfit(
                exp_clean, candy, label=label, **cfg.run_imfit)
            for band in cfg.band_meas:
                cfg.logger.info('running imfit on '+band+'-band cutouts')
                candy = prim.run_imfit(
                    cfg.exp[band], candy, label=label, 
                    master_band=cfg.band_detect, **cfg.run_imfit)
        else:
            cfg.logger.warning('**** no sources satisfy cuts! ****')
    else:
        cfg.logger.warning('**** no sources found! ****')
        candy = Table()

    mask_fracs = utils.calc_mask_bit_fracs(exp_clean)
    cfg.logger.info('task completed in {:.2f} min'.format(cfg.timer))

    if debug_return:
        results = lsst.pipe.base.Struct(all_detections=all_detections,
                                        sources=sources,
                                        candy=candy,
                                        exp=cfg.exp,
                                        exp_clean=exp_clean,
                                        mask_fracs=mask_fracs)
    else:
        results = lsst.pipe.base.Struct(all_detections=all_detections,
                                        sources=sources, 
                                        candy=candy,
                                        mask_fracs=mask_fracs)
        cfg.reset_mask_planes()

    return results
