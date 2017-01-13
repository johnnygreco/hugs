from __future__ import division, print_function

import os
import numpy as np
import lsst.pipe.base
from astropy.table import Table
from . import utils
from . import imtools
from . import primitives as prim
from .cattools import cutter

__all__ = ['run']

sex_io = os.environ.get('SEX_IO_DIR')

def run_use_sex(cfg, debug_return=False, synth_factory=None):
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

    assert cfg.data_id is not None, 'No data id!'
    cfg.timer # start timer

    # don't need this mask for sextractor run
    #utils.remove_mask_planes(cfg.mask, ['DETECTED'])

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
    fpset_smooth = prim.image_threshold(mi_clean_smooth, plane_name='SMOOTHED',
                                     mask=mask_clean, **cfg.thresh_det)
    fpset_smooth.setMask(cfg.mask, 'SMOOTHED')

    ############################################################
    # Find and remove "obvious" blends 
    ############################################################

    cfg.logger.info('finding obvious blends')
    prim.find_blends(exp_clean, fpset_smooth, 
                     plane_name='SMOOTHED', **cfg.find_blends)

    back_rms = mi_clean.getImage().getArray()[mask_clean.getArray()==0].std()
    shape = mask_clean.getArray().shape
    noise_array = back_rms*cfg.rng.randn(shape[0], shape[1])
    replace = mask_clean.getArray() & mask_clean.getPlaneBitMask('BLEND') != 0
    mi_clean.getImage().getArray()[replace] = noise_array[replace]

    ############################################################
    # Detect sources and measure props with SExtractor
    ############################################################

    cfg.logger.info('detecting final sources with sextractor')
    p1, p2 = cfg.data_id['patch'][0], cfg.data_id['patch'][-1]
    label = '{}-{}-{}'.format(cfg.data_id['tract'], p1, p2)
    sources = prim.sex_measure(exp_clean, label)

    ############################################################
    # Cut catalog and run imfit on cutouts
    ############################################################

    cfg.logger.info('cutting catalog and running imfit on cutouts')
    sources_big = cutter(sources.to_pandas(), min_cuts={'FWHM_IMAGE':25})
    sources_big = Table.from_pandas(sources_big)

    sources_big = prim.run_imfit(exp_clean, sources_big, label=label)

    cfg.logger.info('task completed in {:.2f} min'.format(cfg.timer))

    mask_fracs = utils.calc_mask_bit_fracs(exp_clean)
        
    if debug_return:
        # detection plane was modified by find_blends
        results = lsst.pipe.base.Struct(sources=sources,
                                        candy=sources_big,
                                        exposure=cfg.exp,
                                        exp_clean=exp_clean,
                                        mi_clean_smooth=mi_clean_smooth,
                                        mask_fracs=mask_fracs)
    else:
        results = lsst.pipe.base.Struct(sources=sources, 
                                        candy=sources_big,
                                        mask_fracs=mask_fracs)
        cfg.reset_mask_planes()

    return results
