from __future__ import division, print_function

import os
import numpy as np
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
from astropy.io import fits
from astropy.table import Table, vstack, hstack
from scipy.spatial.distance import cdist
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
import sep
from . import utils
from . import imtools
from . import sextractor

__all__ = [
    'clean', 
    'image_threshold', 
    'detect_sources',
]

def clean(exposure, fpset_low, min_pix_low_thresh=100, name_high='THRESH_HIGH', 
          max_frac_high_thresh=0.3, rgrow=None, random_state=None):
    """
    Clean image of bright sources and associated diffuse regions by 
    replacing them with sky noise. Also, remove low-threshold regions 
    that are small. 

    Parameters
    ----------
    exposure : lsst.afw.image.ExposureF
        Exposure object with masks from hugsPipe.run.
    fpset_low : lsst.afw.detection.FootprintSet 
        Low threshold footprints.
    min_pix_low_thresh : int, optional
        Minimum number of pixels for a low-thresh footprint. 
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
        else:
            if fp.getArea() < min_pix_low_thresh:
                fp_list.append(fp)
    fpset_replace.setFootprints(fp_list)
    if rgrow:
        fpset_replace = afwDet.FootprintSet(fpset_replace, rgrow, True)
    mask.addMaskPlane('CLEANED')
    fpset_replace.setMask(mask, 'CLEANED')

    # create new exposure and replace footprints with noise
    exp_clean = exposure.clone()
    mi_clean = exp_clean.getMaskedImage()
    replace = mask.getArray() & mask.getPlaneBitMask('BRIGHT_OBJECT') != 0
    replace |= mask.getArray() & mask.getPlaneBitMask('CLEANED') != 0
    mi_clean.getImage().getArray()[replace] = noise_array[replace]
    return exp_clean


def image_threshold(masked_image, thresh=3.0, thresh_type='stdev', npix=1, 
                    rgrow=None, isogrow=False, plane_name='', mask=None,
                    clear_mask=True):
    """
    Image thresholding. A bit mask will be set with name 'plane_name'.

    Parameters
    ----------
    masked_image : lsst.afw.image.imageLib.MaskedImageF
        A masked image object.
    thresh : float
        Threshold value.
    thresh_type : string, optional
        Threshold type: stdev, pixel_stdev, bitmask, value,
        or variance.
    npix : int, optional
        Minimum number of touching pixels in an object.
    rgrow : int, optional
        Number of pixels to grow footprints.
    isogrow : bool, optional
        If True, use (expensive) isotropic grow. 
    plane_name : string, optional
        Name of bit plane.
    mask : lsst.afw.image.imageLib.MaskU, optional
        Mask to set if not same as in masked_imaged
    clear_mask : bool, optional
        If True, clear the bit plane before thresholding

    Returns
    -------
    fpset : lsst.afw.detection.detectionLib.FootprintSet
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


def detect_sources(exp, sex_config, sex_params, sex_io_dir, dual_exp=None, 
                   delete_created_files=True, label='hugs'):
    """
    Source detection using SExtractor.

    Parameters
    ----------
    exp : lsst.afw.image.ExposureF
        Exposure object to run sextractor on.
    label : str
        Label for this run.

    Returns
    -------
    cat : astropy.table.Table
        Sextractor catalog.
    """

    sw = sextractor.Wrapper(sex_config, sex_params, sex_io_dir)

    #########################################################
    # write exposure for sextractor input and run
    #########################################################

    detect_band = exp.getFilter().getName().lower()
    # some bands have numbers --> get the relevant letter
    detect_band = [b for b in detect_band if b in 'gri'][0] 
    exp_fn = sw.get_io_dir('exp-{}-{}.fits'.format(label, detect_band))
    exp.writeFits(exp_fn)

    if dual_exp is not None:
        meas_band = dual_exp.getFilter().getName().lower()
        meas_band = [b for b in meas_band if b in 'gri'][0]
        dual_fn = sw.get_io_dir('exp-{}-{}.fits'.format(label, meas_band))
        dual_exp.writeFits(dual_fn)
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
        # change mag param names and add SB within each aperture
        #########################################################

        if 'MAG_APER' in cat.colnames:
            cat.rename_column('MAG_APER', 'MAG_APER_0')
            for i, diam in enumerate(sex_config['PHOT_APERTURES'].split(',')):
                r = utils.pixscale*float(diam)/2 # arcsec
                sb = cat['MAG_APER_'+str(i)] + 2.5*np.log10(np.pi*r**2)
                cat['mu_aper_'+str(i)] = sb
      
        #########################################################
        # add band to non-position and shape parameter names
        #########################################################

        detect_band_only = [
            'X_IMAGE', 'Y_IMAGE', 'ALPHA_J2000', 'DELTA_J2000', 
            'THETA_IMAGE', 'A_IMAGE', 'B_IMAGE', 'ELLIPTICITY', 'KRON_RADIUS']

        for name in cat.colnames:
            if name not in detect_band_only: 
                cat.rename_column(name, name+'({})'.format(meas_band))

        #########################################################
        # only save positions from the primary detection band
        #########################################################

        if meas_band==detect_band:
            x0, y0 = exp.getXY0()
            cat['x_img'] = cat['X_IMAGE'] - 1
            cat['y_img'] = cat['Y_IMAGE'] - 1
            cat['x_hsc'] = cat['X_IMAGE'] + x0 - 1
            cat['y_hsc'] = cat['Y_IMAGE'] + y0 - 1
        else:
            cat.remove_columns(detect_band_only)

    #########################################################
    # delete files created by and for sextractor
    #########################################################

    if delete_created_files:
        if dual_exp is not None:
            os.remove(dual_fn)
        os.remove(exp_fn)
        os.remove(cat_fn)
        os.remove(sw.config['PARAMETERS_NAME'])

    return cat
