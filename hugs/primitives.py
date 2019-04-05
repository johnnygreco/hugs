"""
Primitive image processing functions for the hugs pipeline.
"""
from __future__ import division, print_function

import os
import numpy as np
import lsst.afw.detection as afwDet
from .utils import pixscale
from . import utils
from . import sextractor
from .morphology import Morphology

import sep
from scipy import ndimage
from astropy.io import fits
from astropy.table import Table
from astropy.nddata import Cutout2D
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_sigma_to_fwhm, gaussian_fwhm_to_sigma

dustmap = utils.get_dust_map()

__all__ = [
    'image_threshold', 
    'clean', 
    'clean_use_hsc_mask',
    'remove_small_sources_thresholding', 
    'detect_sources',
    'measure_morphology_metrics'
]


def image_threshold(masked_image, thresh=3.0, thresh_type='stdev', npix=1, 
                    rgrow=None, isogrow=False, plane_name='', mask=None,
                    clear_mask=True):
    """
    Image thresholding. A bit mask will be set with name 'plane_name'.

    Parameters
    ----------
    masked_image : lsst.afw.image.MaskedImageF
        A masked image object.
    thresh : float
        Threshold value.
    thresh_type : string, optional
        Threshold type: stdev, pixel_stdev, bitmask, value, or variance.
    npix : int, optional
        Minimum number of touching pixels in an object.
    rgrow : int, optional
        Number of pixels to grow footprints.
    isogrow : bool, optional
        If True, use (expensive) isotropic grow. 
    plane_name : string, optional
        Name of bit plane.
    mask : lsst.afw.image.MaskU, optional
        Mask to set if not same as in masked_imaged
    clear_mask : bool, optional
        If True, clear the bit plane before thresholding

    Returns
    -------
    fpset : lsst.afw.detection.FootprintSet
        Footprints associated with detected objects.
    """

    mask = masked_image.getMask() if mask is None else mask
    thresh_type = getattr(afwDet.Threshold, thresh_type.upper())
    thresh = afwDet.Threshold(thresh, thresh_type)
    fpset = afwDet.FootprintSet(masked_image, thresh, plane_name, npix)
    if rgrow is not None:
        fpset = afwDet.FootprintSet(fpset, rgrow, isogrow)
    if plane_name:
        mask.addMaskPlane(plane_name)
        if clear_mask:
            mask.clearMaskPlane(mask.getMaskPlane(plane_name))
        fpset.setMask(mask, plane_name)
    return fpset


def clean(exposure, fpset_low, name_high='THRESH_HIGH', 
          max_frac_high_thresh=0.15, rgrow=None, random_state=None, 
          bright_object_mask=True):
    """
    Clean image of bright sources and associated diffuse regions by 
    replacing them with sky noise. Also, remove low-threshold regions 
    that are smaller than min_pix_low_thresh.

    Parameters
    ----------
    exposure : lsst.afw.image.ExposureF
        Exposure object with masks from hugsPipe.run.
    fpset_low : lsst.afw.detection.FootprintSet 
        Low threshold footprints.
    name_high : string, optional
        The name of the high-threshold bit plane.
    max_frac_high_thresh : float, optional
        Maximum fraction of high-thresh pixels to keep footprint. 
    rgrow : int, optional
        Number of pixels to grow footprints.
    random_state : int, list of ints, RandomState instance, or None, optional 
        If int or list of ints, random_state is the rng seed.
        If RandomState instance, random_state is the rng.
        If None, the rng is the RandomState instance used by np.random.

    Returns
    -------
    exp_clean : lsst.afw.ExposureF 
        The cleaned exposure object. 
    """
    
    # generate array of gaussian noise
    mi = exposure.getMaskedImage()
    mask = mi.getMask()
    noise_array = utils.make_noise_image(mi, random_state)

    # associate high thresh with low thresh and find small fps
    fpset_replace = afwDet.FootprintSet(mi.getBBox())
    fp_list = []
    for fp in fpset_low.getFootprints():
        hfp = afwDet.HeavyFootprintF(fp, mi)
        pix = hfp.getMaskArray()
        bits_hi = (pix & mask.getPlaneBitMask(name_high)!=0).sum()
        if bits_hi>0:
            ratio = bits_hi/float(fp.getArea())
            if ratio > max_frac_high_thresh:
                fp_list.append(fp)
    fpset_replace.setFootprints(fp_list)
    if rgrow:
        fpset_replace = afwDet.FootprintSet(fpset_replace, rgrow, True)
    mask.addMaskPlane('CLEANED')
    fpset_replace.setMask(mask, 'CLEANED')

    # create new exposure and replace footprints with noise
    exp_clean = exposure.clone()
    mi_clean = exp_clean.getMaskedImage()
    replace = mask.getArray() & mask.getPlaneBitMask('CLEANED') != 0
    if bright_object_mask:
        replace |= mask.getArray() & mask.getPlaneBitMask('BRIGHT_OBJECT') != 0
    mi_clean.getImage().getArray()[replace] = noise_array[replace]

    return exp_clean


def clean_use_hsc_mask(exposure, ref_plane='THRESH_HIGH', rgrow=None, 
                       random_state=None, bright_object_mask=True):

    # generate array of gaussian noise
    mi = exposure.getMaskedImage()
    mask = mi.getMask()
    noise_array = utils.make_noise_image(mi, random_state)
    
    threshold = afwDet.Threshold(mask.getPlaneBitMask(['DETECTED']))
    fp_det = afwDet.FootprintSet(mask, threshold, afwDet.Threshold.BITMASK)

    fp_list = []
    for fp in fp_det.getFootprints():
        hfp = afwDet.HeavyFootprintF(fp, mi)
        pix = hfp.getMaskArray()
        check = (pix & mask.getPlaneBitMask(ref_plane)!=0).sum()
        if check > 0:
            fp_list.append(fp)     
    fpset_replace = afwDet.FootprintSet(mi.getBBox())
    fpset_replace.setFootprints(fp_list)
    if rgrow:
        fpset_replace = afwDet.FootprintSet(fpset_replace, rgrow, True)
    mask.addMaskPlane('CLEANED')
    fpset_replace.setMask(mask, 'CLEANED')

    exp_clean = exposure.clone()
    mi_clean = exp_clean.getMaskedImage()
    replace = mask.getArray() & mask.getPlaneBitMask('CLEANED') != 0
    if bright_object_mask:
        replace |= mask.getArray() & mask.getPlaneBitMask('BRIGHT_OBJECT') != 0
    mi_clean.getImage().getArray()[replace] = noise_array[replace]

    return exp_clean


def remove_small_sources_thresholding(exposure, min_radius_arcsec, 
                                      random_state=None):

    mi = exposure.getMaskedImage()
    mask = mi.getMask()
    noise_array = utils.make_noise_image(mi, random_state)

    threshold = afwDet.Threshold(mask.getPlaneBitMask(['DETECTED']))
    fp_det = afwDet.FootprintSet(mask, threshold, afwDet.Threshold.BITMASK)
    area_min = np.pi * (min_radius_arcsec / utils.pixscale)**2 

    fp_list = []
    for fp in fp_det.getFootprints():
        if fp.getArea() < area_min:
            fp_list.append(fp)     
    fp_small = afwDet.FootprintSet(mi.getBBox())
    fp_small.setFootprints(fp_list)
    mask.addMaskPlane('SMALL')
    fp_small.setMask(mask, 'SMALL')

    exp_clean = exposure.clone()
    mi_clean = exp_clean.getMaskedImage()
    replace = mask.getArray() & mask.getPlaneBitMask('SMALL') != 0
    mi_clean = exp_clean.getMaskedImage()
    mi_clean.getImage().getArray()[replace] = noise_array[replace]

    return exp_clean


def detect_sources(exp, sex_config, sex_io_dir, dual_exp=None, 
                   delete_created_files=True, label='hugs',
                   original_fn=None):
    """
    Source detection using SExtractor.

    Parameters
    ----------
    exp : lsst.afw.image.ExposureF
        Exposure object to run sextractor on.
    sex_config : dict
        SExtractor configuration that is is different for its default.
    sex_io_dir : str
        Input/output directory for files needed to run SExtractor.
    dual_exp : lsst.afw.image.ExposureF, optional
        Dual exposure for forced photometry. 
    delete_created_files : bool, optional
        If True, remove all files created for and by SExtractor.
    label : str
        Label for this run.

    Returns
    -------
    cat : astropy.table.Table
        SExtractor catalog.
    """

    sw = sextractor.Wrapper(sex_config, sex_io_dir)

    #########################################################
    # write exposure for sextractor input and run
    #########################################################

    detect_band = exp.getFilter().getName().lower()
    # some bands have numbers --> get the relevant letter
    detect_band = [b for b in detect_band if b in 'gri'][0] 
    exp_fn = sw.get_io_dir('exp-{}-{}.fits'.format(label, detect_band))

    # HACK: work around strange new bug related to SEXtractor
    if original_fn is None:
        exp.writeFits(exp_fn)
    else:
        header = fits.getheader(original_fn, ext=1)
        fits.writeto(exp_fn, exp.getImage().getArray(), header, overwrite=True)

    if dual_exp is not None:
        meas_band = dual_exp.getFilter().getName().lower()
        meas_band = [b for b in meas_band if b in 'gri'][0]
        dual_fn = sw.get_io_dir('exp-{}-{}.fits'.format(label, meas_band))

        # HACK: work around strange new bug related to SEXtractor
        if original_fn is None:
            dual_exp.writeFits(dual_fn)
        else:
            fn = original_fn.replace('HSC-'+detect_band.upper(), 
                                     'HSC-'+meas_band.upper())
            header = fits.getheader(fn, ext=1)
            fits.writeto(dual_fn, dual_exp.getImage().getArray(), header, 
                         overwrite=True)

        run_fn = exp_fn+'[1],'+dual_fn+'[1]'
        cat_label = 'sex-{}-{}-{}'.format(label, detect_band, meas_band)
    else:
        meas_band = detect_band
        cat_label = 'sex-{}-{}'.format(label, detect_band)
        run_fn = exp_fn+'[1]'

    cat_fn = sw.get_io_dir(cat_label+'.cat')

    #########################################################
    # run SExtactor and get catalog
    #########################################################

    sw.run(run_fn, cat_fn=cat_fn)
    cat = sextractor.read_cat(sw.get_io_dir(cat_fn))

    if len(cat)>0:

        #########################################################
        # only save positions from the primary detection band
        #########################################################

        detect_band_only = [
            'X_IMAGE', 'Y_IMAGE', 'ALPHA_J2000', 'DELTA_J2000', 'FLAGS', 
            'PETRO_RADIUS', 'THETA_IMAGE', 'A_IMAGE', 'B_IMAGE', 
            'ELLIPTICITY', 'KRON_RADIUS'
        ]

        ebv = dustmap.ebv(cat['ALPHA_J2000'], cat['DELTA_J2000'])

        if meas_band==detect_band:
            x0, y0 = exp.getXY0()
            cat['X_IMAGE'] -= 1
            cat['Y_IMAGE'] -= 1
            cat['X_HSC'] = cat['X_IMAGE'] + x0 
            cat['Y_HSC'] = cat['Y_IMAGE'] + y0 
            detect_band_only.append('X_HSC')
            detect_band_only.append('Y_HSC')
        else:
            cat.remove_columns(detect_band_only)

        #########################################################
        # rename columns, change units of flux_radius and 
        # fwhm_image to arcsec, add extinction params
        #########################################################

        cat.rename_column('MAG_APER', 'MAG_APER_0')
        cat.rename_column('MAGERR_APER', 'MAGERR_APER_0')
        cat.rename_column('FLUX_RADIUS', 'FLUX_RADIUS_0')

        for i, diam in enumerate(sex_config['PHOT_APERTURES'].split(',')):
            cat.rename_column('MAG_APER_'+str(i), 'mag_ap'+str(i))
            cat.rename_column('MAGERR_APER_'+str(i), 'magerr_ap'+str(i))
        
        for i, frac in enumerate(sex_config['PHOT_FLUXFRAC'].split(',')):
            frac = str(int(100 * float(frac)))
            cat.rename_column('FLUX_RADIUS_'+str(i), 'flux_radius_'+frac)
            cat['flux_radius_'+frac] *= utils.pixscale

        cat['FWHM_IMAGE'] = cat['FWHM_IMAGE']*utils.pixscale
        cat.rename_column('FWHM_IMAGE', 'FWHM')

        for name in cat.colnames:
            if name not in detect_band_only: 
                cat.rename_column(name, name.lower()+'_'+meas_band)
            else:
                cat.rename_column(name, name.lower())
        if meas_band==detect_band:
            cat.rename_column('alpha_j2000', 'ra')
            cat.rename_column('delta_j2000', 'dec')
            cat['ebv'] = ebv
        cat['A_'+meas_band] = ebv*getattr(utils.ext_coeff, meas_band)

    #########################################################
    # delete files created by and for sextractor
    #########################################################

    if delete_created_files:
        if dual_exp is not None:
            os.remove(dual_fn)
        os.remove(exp_fn)
        os.remove(cat_fn)

    return cat


def measure_morphology_metrics(image, sources, scale_size=3):
    gini_full = []
    gini_1 = []
    gini_1p5 = []
    gini_2 = []
    gini_1p5_circ = []
    gini_2_circ = []

    acorr_peak = []
    acorr_bkgd = []

    for src in sources:
        try:
            centroid = [src['x_image'], src['y_image']]
            cutout = Cutout2D(image, 
                              centroid, 
                              src['flux_radius_65_g'] * scale_size / pixscale) 
            x_shift = cutout.center_original[0] - cutout.center_cutout[0]
            y_shift = cutout.center_original[1] - cutout.center_cutout[1]
            centroid = [src['x_image'] - x_shift, src['y_image'] - y_shift]

            morph = Morphology(cutout.data, 
                               centroid, 
                               src['a_image'], 
                               src['b_image'], 
                               src['ellipticity'], 
                               src['theta_image'])

            _, circ_mean, ann_mean = morph.autocorr(1.5)
            acorr_peak.append(circ_mean)
            acorr_bkgd.append(ann_mean)

            gini_full.append(morph.gini(0))
            gini_1.append(morph.gini(1))
            gini_1p5.append(morph.gini(1.5))
            gini_2.append(morph.gini(2))
            gini_1p5_circ.append(morph.gini(1.5, 3))
            gini_2_circ.append(morph.gini(2, 3))

        except:
            acorr_peak.append(np.nan)
            acorr_bkgd.append(np.nan)
            gini_full.append(np.nan)
            gini_1.append(np.nan)
            gini_1p5.append(np.nan)
            gini_2.append(np.nan)
            gini_1p5_circ.append(np.nan)
            gini_2_circ.append(np.nan)

    sources['gini_full'] = gini_full
    sources['gini_1'] = gini_1
    sources['gini_1p5'] = gini_1p5
    sources['gini_2'] = gini_2
    sources['gini_1p5_circ'] = gini_1p5_circ
    sources['gini_2_circ'] = gini_2_circ

    sources['acorr_peak'] = acorr_peak
    sources['acorr_bkgd'] = acorr_bkgd
    acorr_bkgd =  np.array(acorr_bkgd)
    acorr_bkgd[acorr_bkgd==0] = np.nan
    sources['acorr_ratio'] = np.array(acorr_peak) / acorr_bkgd
