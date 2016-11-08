from __future__ import division, print_function

import numpy as np
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
from astropy.table import Table, vstack
from astropy.convolution import Gaussian2DKernel
import photutils as phut
from . import utils
from . import imtools

__all__ = ['clean', 
           'find_blends',
           'image_threshold', 
           'measure_sources', 
           'photometry']


def clean(exposure, fpset_low, min_pix_low_thresh=100, name_high='THRESH_HIGH', 
          max_frac_high_thresh=0.3, rgrow=None):
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

    Returns
    -------
    exp_clean : lsst.afw.ExposureF 
        The cleaned exposure object. 
    """
    
    # generate array of gaussian noise
    mi = exposure.getMaskedImage()
    mask = mi.getMask()
    shape = mask.getArray().shape
    back_rms = mi.getImage().getArray()[mask.getArray()==0].std()
    noise_array = back_rms*np.random.randn(shape[0], shape[1])

    # associate high thresh with low thresh and find small fps
    fpset_replace = afwDet.FootprintSet(mi.getBBox())
    fp_list = afwDet.FootprintList()
    for fp in fpset_low.getFootprints():
        hfp = afwDet.HeavyFootprintF(fp, mi)
        pix = hfp.getMaskArray()
        bits_hi = (pix & mask.getPlaneBitMask(name_high) != 0).sum()
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


def find_blends(exp, fpset_det, name_low='THRESH_LOW', num_low_fps=2,
                min_fp_area=100, min_parent_area=0.3, min_parent_flux=0.8):
    """
    Search footprints for obvious blends of multiple small sources. 
    The search is carried out by counting the number of low-thresh 
    footprints within each detection footprint and looking at properties
    of the low-thresh footprints (e.g., size and flux).

    Parameters
    ----------
    exp : lsst.afw.image.ExposureF
        Exposure object. 
    fpset_det : lsst.afw.FootprintSet
        Detection footprint set.
    name_low : string, optional
        Name of the low-thresh bit plane.
    num_low_fps : int, optional
        Number of low-thresh fps to be considered blend.
    min_fp_area : int, optional
        Minimum number of pixels to be considered an object.
    min_parent_area : float, optional
        Minimum fraction of the footprint occupied by parent.
    min_parent_flux : float, optional
        Minimum fraction of the total flux that must 
        be contributed by the parent. 

    Notes
    -----
    A new bit plane called BLEND will be added to the mask 
    of the input exposure object, and DETECTED bits will be 
    turned off for blended sources. 
    """

    mask = exp.getMaskedImage().getMask()
    plane_low = mask.getPlaneBitMask(name_low)
    fp_list = afwDet.FootprintList()
    
    for fp_det in fpset_det.getFootprints():

        # get bbox of detection footprint
        bbox = fp_det.getBBox()
        exp_fp = exp.Factory(exp, bbox, afwImage.PARENT) 
        mask_fp = exp_fp.getMaskedImage().getMask()
        mi_fp = exp_fp.getMaskedImage()
        total_area = float(fp_det.getArea())

        # find low-thresh footprints in bbox
	fpset_low = afwDet.FootprintSet(
            mask_fp, afwDet.Threshold(plane_low, afwDet.Threshold.BITMASK))

        # find blends and save footprint
        if len(fpset_low.getFootprints()) >= num_low_fps:
            areas = np.array([fp_.getArea() for fp_ in
                              fpset_low.getFootprints()])
            biggest = areas.argmax()
            if (areas > min_fp_area).sum() >= num_low_fps:
                hfps = [afwDet.HeavyFootprintF(fp_, mi_fp)\
                        for fp_ in fpset_low.getFootprints()]
                fluxes = np.array(
                    [hfp_.getImageArray().sum()for hfp_ in hfps])
                is_blend = areas[biggest]/total_area < min_parent_area
                is_blend  |= fluxes[biggest]/fluxes.sum() < min_parent_flux
                if is_blend:
                    fp_list.append(fp_det)
                    mask_fp.clearMaskPlane(mask_fp.getMaskPlane('DETECTED'))
            elif (areas < min_fp_area).sum() == areas.size:
                fp_list.append(fp_det)
                mask_fp.clearMaskPlane(mask_fp.getMaskPlane('DETECTED'))


    # create footprint set and set mask plane
    fpset_blends = afwDet.FootprintSet(exp.getBBox())
    fpset_blends.setFootprints(fp_list)
    mask.addMaskPlane('BLEND')
    fpset_blends.setMask(mask, 'BLEND')


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


def measure_sources(exposure, npix=5, thresh_snr=0.5, kern_sig_pix=3, 
                   grow_stamps=None, detect_kwargs={}, deblend_kwargs={}):
    """
    Use photutils to measure sources within "detection" footprints.

    Parameters
    ----------
    exposure : lsst.afw.ExposureF
        Exposure object with masks from hugsPipe.run.
    thresh_snr : float, optional
        The signal-to-noise ratio per pixel above the background 
        for which to consider a pixel as possibly being part of a source.
    kern_sig_pix : float
        Sigma (in pixels) of Gaussian kernel used to smooth detection image.
    npix : int, optional
        Number of pixels above required to be an object.
    grow_stamps : int
        Number of pixels to grow postage stamp in all directions.
    detect_kwargs :  dict, optional
        Kwargs for photutils.detect_sources.
    deblend_kwargs :  dict, optional
        Kwargs for photutils.deblend_sources.

    Returns
    -------
    table : astropy.table.Table
        Source measurements. 
    """

    # get detected footprints
    mask = exposure.getMaskedImage().getMask()
    plane = mask.getPlaneBitMask('DETECTED')
    fpset = afwDet.FootprintSet(
        mask, afwDet.Threshold(plane, afwDet.Threshold.BITMASK))

    # loop over footprints
    table = Table()
    kern = Gaussian2DKernel(kern_sig_pix)
    kern.normalize()
    fp_id = 1
    for fp in fpset.getFootprints():

        # get exposure in bbox of footprint 
        bbox = fp.getBBox()
        if grow_stamps:
            bbox.grow(grow_stamps)
            bbox.clip(exposure.getBBox())
        exp = exposure.Factory(exposure, bbox, afwImage.PARENT)

        # detect sources in footprint
        img = exp.getMaskedImage().getImage().getArray().copy()
        x0, y0 = exp.getXY0()
        thresh = phut.detect_threshold(img, 
                                       snr=thresh_snr, 
                                       **detect_kwargs)
        seg = phut.detect_sources(img, thresh, npixels=npix, 
                                  filter_kernel=kern)

        # go to next footprint if no sources are detected
        if seg.nlabels==0: 
            continue

        # deblend sources
        seg_db = phut.deblend_sources(
            img, seg, npixels=npix, 
            filter_kernel=kern, **deblend_kwargs)

        # measure source properties
        props = phut.source_properties(img, seg_db)
        props = phut.properties_table(props)
        props['x_hsc'] = props['xcentroid'] + x0
        props['y_hsc'] = props['ycentroid'] + y0
        props['fp_id'] = [fp_id]*len(props)
        fp_id += 1
        if len(props)>1:
            props['is_deblended'] = [True]*len(props)
        else:
            props['is_deblended'] = False
        table = vstack([table, props])
    
    # clean up table columns
    table['id'] = np.arange(1, len(table)+1)
    remove_cols = [
        'ra_icrs_centroid', 'dec_icrs_centroid', 'source_sum_err', 
        'background_sum', 'background_mean', 'background_at_centroid']
    table.remove_columns(remove_cols)

    # use wcs to add ra & dec to table
    ra_list = []
    dec_list = []
    wcs = exposure.getWcs()
    for x, y in table['x_hsc', 'y_hsc']:
        ra = wcs.pixelToSky(x, y).getLongitude().asDegrees()
        dec = wcs.pixelToSky(x, y).getLatitude().asDegrees()
        ra_list.append(ra)
        dec_list.append(dec)
    table['ra'] = ra_list
    table['dec'] = dec_list
    
    # calculate primary array coordinates
    table['x_img'] = table['x_hsc'] - exposure.getX0()
    table['y_img'] = table['y_hsc'] - exposure.getY0()
    return table


def photometry(img, sources, zpt_mag=27.0, ell_nsig=4.0, 
               circ_radii=[3, 6, 9]):
    """
    Do basic aperture photometry within circular apertures and 
    an elliptical aperture. Results will be saved to the 
    input table.

    Parameters
    ----------
    img : 2D ndarray
        The image data.
    sources : astropy.table.Table
        Output table from measure_sources.
    zpt_mag : float, optional
        The zero point magnitude.
    ell_nsig : float, optional
        Number of sigma for major/minor axis sizes.
    circ_radii : list, optional
        Radii within which to circular aperture photometry.
    """

    pos = [(x, y) for x,y in sources['x_img', 'y_img']] 

    mag = []
    for idx in range(len(sources)):
        x, y = pos[idx]
        a, b = sources['semimajor_axis_sigma', 'semiminor_axis_sigma'][idx]
        theta = sources['orientation'][idx]
        aperture = phut.EllipticalAperture((x,y), ell_nsig*a, ell_nsig*b, theta)
        flux = phut.aperture_photometry(img, aperture)['aperture_sum'][0]
        mag.append(zpt_mag - 2.5*np.log10(flux))
    sources['mag_ell'] = mag

    for r in circ_radii:
        r_arcsec = r*utils.pixscale
        apertures = phut.CircularAperture(pos, r=r)
        flux = phut.aperture_photometry(img, apertures)['aperture_sum']
        mag = zpt_mag - 2.5*np.log10(flux)
        sources['mag_circ_{}'.format(r)] = mag
        mu = mag + 2.5*np.log10(np.pi*r_arcsec**2)
        sources['mu_{}'.format(r)] = mu
