from __future__ import division, print_function

import os
import numpy as np
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
from astropy.io import fits
from astropy.table import Table, vstack, hstack
from astropy.convolution import Gaussian2DKernel
from scipy.spatial.distance import cdist
import photutils as phut
from . import utils
from . import imtools

__all__ = [
    'clean', 
    'find_blends',
    'image_threshold', 
    'measure_sources', 
    'photometry', 
    'run_imfit',
    'sex_measure'
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
    fp_list = afwDet.FootprintList()
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


def find_blends(exp, fpset_det, name_low='THRESH_LOW', num_low_fps=2,
                min_fp_area=100, min_parent_area=0.3, min_parent_flux=0.8, 
                plane_name='DETECTED'):
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
    plane_name : str
        Name of bit plane to search for blends.

    Notes
    -----
    A new bit plane called BLEND will be added to the mask 
    of the input exposure object, and DETECTED bits will be 
    turned off for blended sources. 
    """

    mi = exp.getMaskedImage()
    mask = mi.getMask()
    msk_arr = mask.getArray()
    plane_low = mask.getPlaneBitMask(name_low)
    fp_list = afwDet.FootprintList()
    x0, y0 = exp.getXY0()
    
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
                    for s in fp_det.getSpans():
                        y = s.getY() - y0
                        for x in range(s.getX0()-x0, s.getX1()-x0+1):
                            msk_arr[y][x] -= mask.getPlaneBitMask(plane_name)
            elif (areas < min_fp_area).sum() == areas.size:
                fp_list.append(fp_det)
                for s in fp_det.getSpans():
                    y = s.getY() - y0
                    for x in range(s.getX0()-x0, s.getX1()-x0+1):
                        msk_arr[y][x] -= mask.getPlaneBitMask(plane_name)

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
                    grow_stamps=None, detect_kwargs={}, deblend_kwargs={},
                    logger=None, sf=None):
    """
    Use photutils to measure sources within "detection" footprints.

    Parameters
    ----------
    exposure : lsst.afw.ExposureF
        Exposure object with masks from hugs_pipe.run.
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
    sf : hugs_pipe.synths.SynthFactory, optional
        synth factory object for synth flag.

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

    # don't need these photutils columns
    exclude_cols = [
        'ra_icrs_centroid', 'dec_icrs_centroid', 'source_sum_err', 
        'background_sum', 'background_mean', 'background_at_centroid',
        'xmin', 'xmax', 'ymin', 'ymax', 'min_value',
        'max_value', 'minval_xpos', 'minval_ypos', 'maxval_xpos',
        'maxval_ypos'
    ]

    for fp_id, fp in enumerate(fpset.getFootprints()):

        # get exposure in bbox of footprint 
        bbox = fp.getBBox()
        if grow_stamps:
            bbox.grow(grow_stamps)
            bbox.clip(exposure.getBBox())
        exp_fp = exposure.Factory(exposure, bbox, afwImage.PARENT)
        mi_fp = exp_fp.getMaskedImage()
        msk_fp = mi_fp.getMask()
        hfp = afwDet.HeavyFootprintF(fp, mi_fp)
        pix = hfp.getMaskArray()

        # detect sources in footprint
        img = mi_fp.getImage().getArray().copy()
        x0, y0 = exp_fp.getXY0()
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
        props = phut.properties_table(props, exclude_columns=exclude_cols)

        # calculate different coordinates
        props['x_hsc'] = props['xcentroid'] + x0
        props['y_hsc'] = props['ycentroid'] + y0
        props['x_img'] = props['x_hsc'] - exposure.getX0()
        props['y_img'] = props['y_hsc'] - exposure.getY0()

        # generate ids and flags
        props['fp_id'] = [fp_id]*len(props)
        props['nchild'] = [len(props)]*len(props)
        num_edge = (pix & msk_fp.getPlaneBitMask('EDGE')!=0).sum()
        num_bright = (pix & msk_fp.getPlaneBitMask('BRIGHT_OBJECT')!=0).sum()
        props['num_edge_pix'] = [num_edge]*len(props)
        props['num_bright_pix'] = [num_bright]*len(props)
        props['fp_det_area'] = [fp.getArea()]*len(props)

        # find synths
        max_npix_dist = 5
        if 'SYNTH' in msk_fp.getMaskPlaneDict().keys() and sf is not None:
            props['synth_id'] = [-111]*len(props)
            props['synth_fail'] = [False]*len(props)
            if (pix & msk_fp.getPlaneBitMask('SYNTH')!=0).sum() > 0:
                synth_coords = sf.get_psets()[['X0', 'Y0']]
                props_coords = props['x_img', 'y_img'].to_pandas()
                dist = cdist(synth_coords, props_coords)
                min_idx = np.unravel_index(dist.argmin(), dist.shape)
                if dist.min() < max_npix_dist:
                    props['synth_id'][min_idx[1]] = min_idx[0]
                    num_pass = (np.min(dist, axis=1) < max_npix_dist).sum()
                    if len(props)>1 and num_pass>1:
                        props['synth_fail'] = [True]*len(props)
                else:
                    props['synth_fail'] = [True]*len(props)

        table = vstack([table, props])
    
    if len(table)>0:
        table['id'] = np.arange(0, len(table))
        table.remove_columns(['xcentroid', 'ycentroid'])

        # use wcs to add ra & dec to table
        ra_list = []
        dec_list = []
        wcs = exposure.getWcs()
        if wcs is None:
            warn = 'no wcs with exposure'
            if logger:
                logger.warning(warn)
            else:
                print(warn)
        else:
            for x, y in table['x_hsc', 'y_hsc']:
                ra = wcs.pixelToSky(x, y).getLongitude().asDegrees()
                dec = wcs.pixelToSky(x, y).getLatitude().asDegrees()
                ra_list.append(ra)
                dec_list.append(dec)
            table['ra'] = ra_list
            table['dec'] = dec_list

        # masked columns make to_pandas convert ints to floats
        table = Table(table, masked=False)
    else:
        warn = '**** no sources found! ****'
        if logger:
            logger.warning(warn)
        else:
            print(warn)

    return table


def photometry(img_data, sources, zpt_mag=27.0, ell_nsig=5.0, 
               circ_radii=[3, 6, 9]):
    """
    Do basic aperture photometry within circular apertures and 
    an elliptical aperture. Results will be saved to the 
    input table.

    Parameters
    ----------
    img_data : 2D ndarray, or dict
        The image data. If a dict is given, the keys
        must be the photometric band
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

    if type(img_data)==np.ndarray:
        img_data = {'i': img_data}

    for band, img in img_data.items():
        mag = []
        for idx in range(len(sources)):
            x, y = pos[idx]
            a, b = sources['semimajor_axis_sigma', 'semiminor_axis_sigma'][idx]
            theta = sources['orientation'][idx]
            aperture = phut.EllipticalAperture((x,y), ell_nsig*a, ell_nsig*b, theta)
            flux = phut.aperture_photometry(img, aperture)['aperture_sum'][0]
            mag.append(zpt_mag - 2.5*np.log10(flux))
        sources['mag_ell_'+band.lower()] = mag

        for r in circ_radii:
            r_arcsec = r*utils.pixscale
            apertures = phut.CircularAperture(pos, r=r)
            flux = phut.aperture_photometry(img, apertures)['aperture_sum']
            mag = zpt_mag - 2.5*np.log10(flux)
            sources['mag_circ_{}_{}'.format(r, band.lower())] = mag
            mu = mag + 2.5*np.log10(np.pi*r_arcsec**2)
            sources['mu_{}_{}'.format(r, band.lower())] = mu


def run_imfit(exp, cat, label='run', clean=True, viz_fit=False,
              psf_convolve=False):
    """
    Run imfit on postage stamps to make cuts on sample.

    Parameters
    ----------
    exp : lsst.afw.image.ExposureF
        Exposure object. 
    cat : astropy.table.Table
        Source catalog (from sextractor run)
    label : string
        Label for this run. 
    clean : bool
        If True, delete files created by this function.
    psf_convolve : bool
        If True, convolve image with psf.

    Returns
    -------
    results : astropy.table.Table
        Source catalog with columns added for fit params.
    """
    from hugs.tasks import sersic_fit
    from hugs.imfit.viz import img_mod_res
    if viz_fit:
        import matplotlib.pyplot as plt

    size = 120
    band = exp.getFilter().getName().lower()
    wcs = exp.getWcs()
    prefix = os.path.join(utils.io, 'temp-io')

    results = Table()

    if psf_convolve:
        psf = exp.getPsf().computeImage().getArray()
        psf_fn = os.path.join(prefix, 'psf-{}.fits'.format(label))
        fits.writeto(psf_fn, psf, clobber=True)
    else:
        psf_fn = None
    
    for obj in cat:
        num = str(obj['NUMBER'])
        coord_hsc = obj['x_hsc'], obj['y_hsc']
        cutout = imtools.get_cutout(coord_hsc, size, exp=exp)
        
        fn = os.path.join(prefix, 'cutout-{}-{}.fits'.format(label, num))
        cutout.writeFits(fn)

        X0 = coord_hsc[0] - cutout.getX0()
        Y0 = coord_hsc[1] - cutout.getY0()

        init_params = {
            'X0': X0, 
            'Y0': Y0,
            'PA': [obj['THETA_IMAGE'] + 90, 0, 180],
            'ell': [obj['ELLIPTICITY'], 0, 0.999],
        }

        fit_prefix = os.path.join(prefix, label)
        fit = sersic_fit(
            fn, init_params=init_params, prefix=fit_prefix, 
            clean='config', psf_fn=psf_fn)

        x0_hsc, y0_hsc = fit.X0 + cutout.getX0(), fit.Y0 + cutout.getY0()
        coord = wcs.pixelToSky(x0_hsc, y0_hsc)
        ra, dec = coord.getPosition(afwGeom.degrees)
        dr0 = np.sqrt((fit.X0 - X0)**2 + (fit.Y0 - Y0)**2)
        dsize = (fit.r_e - obj['FLUX_RADIUS'])*utils.pixscale
        dmu = fit.mu_0 - obj['mu_aper_0']

        data = [ra, dec, fit.n, fit.m_tot, fit.mu_0, fit.ell, fit.X0, fit.Y0,
                fit.r_e*utils.pixscale, fit.PA, x0_hsc, y0_hsc, dr0, dmu, dsize]

        names = ['ra_imfit', 'dec_imfit', 'n', 'm_tot('+band+')', 
                 'mu_0('+band+')', 'ell('+band+')', 'x_img_imfit', 
                 'y_img_imfit', 'r_e('+band+')', 'PA('+band+')', 
                 'x_hsc_imfit', 'y_hsc_imfit', 'dr0', 'dmu', 'dr_e']

        results = vstack([results, Table(rows=[data], names=names)])

        if viz_fit:
            img_mod_res(fn,
                        fit.params, 
                        fit_prefix+'_photo_mask.fits',
                        band='i', 
                        show=True, 
                        subplots=plt.subplots(1, 3, figsize=(15,5)))

        if clean:
            os.remove(fn)
            os.remove(fit_prefix+'_bestfit_params.txt')
            os.remove(fit_prefix+'_photo_mask.fits')

    results = hstack([cat, results])

    return results


def sex_measure(exp, label='run'):
    """
    Perform mesurements using SExtractor (because I'm feeling desperate). 

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
    import sexpy
        
    #########################################################
    # setup configuration and extra save params
    #########################################################

    apertures = [3., 4., 5., 6., 7., 8., 16., 32.]
    param_fn = label+'.params'
    config = {
        'DETECT_THRESH': 0.7, 
        'FILTER_NAME': 'gauss_5.0_21x21.conv',
        'BACK_SIZE': 128,
        'VERBOSE_TYPE': 'NORMAL', 
        'PHOT_APERTURES': ','.join([str(int(a)) for a in apertures]),
        'DETECT_MINAREA': 100,
        'PARAMETERS_NAME': param_fn,
    }
    params = ['MAG_APER('+str(i+1)+')' for i in range(len(apertures))]

    sw = sexpy.SexWrapper(config, params=params)
    sw.set_check_images('s', prefix=label+'-')

    #########################################################
    # write exposure for sextractor input and run
    #########################################################

    exp_fn = 'exp-'+label+'.fits'
    exp.writeFits(sw.get_indir(exp_fn))

    cat_label = 'sex-'+label
    sw.run(exp_fn+'[1]', cat=cat_label+'.cat')

    #########################################################
    # open seg map and add to bit mask
    #########################################################

    seg_fn = sw.get_outdir(label+'-'+'SEGMENTATION.fits')
    seg = fits.getdata(seg_fn)

    mask = exp.getMaskedImage().getMask()
    mask.addMaskPlane('SEX_SEG')
    mask.getArray()[seg>0] += mask.getPlaneBitMask('SEX_SEG')

    #########################################################
    # read catalog, add params, and write to csv file
    #########################################################

    cat = sexpy.read_cat(sw.get_outdir(cat_label+'.cat'))

    if len(cat)>0:
        x0, y0 = exp.getXY0()
        cat['x_img'] = cat['X_IMAGE'] 
        cat['y_img'] = cat['Y_IMAGE']
        cat['x_hsc'] = cat['X_IMAGE'] + x0
        cat['y_hsc'] = cat['Y_IMAGE'] + y0
        cat.rename_column('MAG_APER', 'MAG_APER_0')

        for i, diam in enumerate(apertures):
            r = utils.pixscale*diam/2 # arcsec
            sb = cat['MAG_APER_'+str(i)] + 2.5*np.log10(np.pi*r**2)
            cat['mu_aper_'+str(i)] = sb

        sex_plane = mask.getPlaneBitMask('SEX_SEG')
        edge_plane = mask.getPlaneBitMask('EDGE')
        bright_plane = mask.getPlaneBitMask('BRIGHT_OBJECT')

        fpset = afwDet.FootprintSet(
            mask, afwDet.Threshold(sex_plane, afwDet.Threshold.BITMASK))

        cat['num_edge_pix'] = -1 
        cat['num_bright_pix'] = -1
        cat['parent_fp_area'] = -1
        cat['fp_id'] = -1
        cat['nchild'] = -1
        
        for obj_i in range(len(cat)):
            x, y = cat['x_hsc', 'y_hsc'][obj_i]
            point = afwGeom.Point2I(int(x), int(y))
            for fp_id, fp in enumerate(fpset.getFootprints()):
                if fp.contains(point):
                    hfp = afwDet.HeavyFootprintF(fp, exp.getMaskedImage())
                    pix = hfp.getMaskArray()
                    num_edge = (pix & edge_plane != 0 ).sum()
                    num_bright = (pix & bright_plane != 0).sum()
                    cat['num_edge_pix'][obj_i] = num_edge
                    cat['num_bright_pix'][obj_i] = num_bright
                    cat['parent_fp_area'][obj_i] = hfp.getArea()
                    cat['fp_id'][obj_i] = fp_id
                    break

        for fp_id in np.unique(cat['fp_id']):
            obj_mask = cat['fp_id']==fp_id
            cat['nchild'][obj_mask] = obj_mask.sum()

    #########################################################
    # delete files created by and for sextractor
    #########################################################

    os.remove(sw.get_indir(exp_fn))
    os.remove(sw.get_outdir(cat_label+'.cat'))
    os.remove(sw.get_configdir(param_fn))
    os.remove(seg_fn)

    return cat
