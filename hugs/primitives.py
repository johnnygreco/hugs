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

__all__ = [
    'clean', 
    'image_threshold', 
    'run_imfit',
    'sex_measure',
    'detection',
    'photometry'
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


def run_imfit(exp, cat, label='run', master_band=None, bbox_grow=120, 
              clean=True, save_fit_fig=False, psf_convolve=False, quiet=False, 
              relpath=None):
    """
    Run imfit on postage stamps to make cuts on sample.

    Parameters
    ----------
    exp : lsst.afw.image.ExposureF
        Exposure object. 
    cat : astropy.table.Table
        Source catalog (from sextractor run)
    label : string, optional
        Label for this run. 
    master_band : str, optional 
        Master band for "forced" photometry.
    bbox_grow : int, optional
        Number of pixels to grow bbox in all directions.
    clean : bool, optional
        If True, delete files created by this function.
    psf_convolve : bool, optional
        If True, convolve image with psf.

    Returns
    -------
    results : astropy.table.Table
        Source catalog with columns added for fit params.
    """
    from hugs.tasks import sersic_fit
    if save_fit_fig:
        import matplotlib.pyplot as plt
        from hugs.imfit.viz import img_mod_res

    band = exp.getFilter().getName().lower()
    wcs = exp.getWcs()
    prefix = os.environ.get('TEMP_IO')
    if relpath:
        prefix = os.path.join(prefix, relpath)
    label += '-'+band

    results = Table()

    if psf_convolve:
        psf = exp.getPsf().computeImage().getArray()
        psf_fn = os.path.join(prefix, 'psf-{}.fits'.format(label))
        fits.writeto(psf_fn, psf, clobber=True)
    else:
        psf_fn = None
    
    for num, obj in enumerate(cat):
        coord_hsc = obj['x_hsc'], obj['y_hsc']
        cutout = imtools.get_cutout(coord_hsc, bbox_grow, exp=exp)
        dimension = cutout.getDimensions()
        
        fn = os.path.join(prefix, 'cutout-{}-{}.fits'.format(label, num))
        cutout.writeFits(fn)

        X0 = coord_hsc[0] - cutout.getX0()
        Y0 = coord_hsc[1] - cutout.getY0()
        img_shift = [cutout.getX0() - exp.getX0(), 
                     cutout.getY0() - exp.getY0()] 
        xmax = np.minimum(X0 + bbox_grow, dimension[0]-1)
        ymax = np.minimum(Y0 + bbox_grow, dimension[1]-1)

        if master_band is None:
            init_params = {
                'X0': [X0, 1, dimension[0]-1],
                'Y0': [Y0, 1, dimension[1]-1],
                'PA': [obj['THETA_IMAGE('+band+')'] + 90, 0, 180],
                'ell': [obj['ELLIPTICITY('+band+')'], 0, 0.999],
            }
        else:
            init_params = {
                'X0': [obj['x_img_imfit'] - img_shift[0], 'fixed'], 
                'Y0': [obj['y_img_imfit'] - img_shift[1], 'fixed'], 
                'n': [obj['n'], 'fixed'], 
                'ell': [obj['ell'], 'fixed'], 
                'PA': [obj['PA'], 'fixed'], 
                'r_e': obj['r_e('+master_band+')'], 
                'I_e': obj['I_e('+master_band+')']
            }

        fit_prefix = os.path.join(prefix, label)
        fit = sersic_fit(
            fn, init_params=init_params, prefix=fit_prefix, 
            clean='config', psf_fn=psf_fn, quiet=quiet)

        data = [fit.m_tot, fit.mu_0, fit.r_e*utils.pixscale, fit.I_e]
        names = ['m_tot('+band+')', 'mu_0('+band+')', 
                 'r_e('+band+')', 'I_e('+band+')']

        if 'FLUX_RADIUS('+band+')' in cat.colnames:
            dsize = (fit.r_e - obj['FLUX_RADIUS('+band+')'])*utils.pixscale
            dmu = fit.mu_0 - obj['mu_aper_0('+band+')']
            data.extend([dsize, dmu])
            names.extend(['dr_e('+band+')', 'dmu('+band+')'])

        if master_band is None:
            x0_hsc, y0_hsc = fit.X0 + cutout.getX0(), fit.Y0 + cutout.getY0()
            coord = wcs.pixelToSky(x0_hsc, y0_hsc)
            ra, dec = coord.getPosition(afwGeom.degrees)
            dR0 = np.sqrt((fit.X0 - X0)**2 + (fit.Y0 - Y0)**2)
            
            x_img = fit.X0 + img_shift[0]
            y_img = fit.Y0 + img_shift[1]

            master_data = [ra, dec, fit.n, fit.ell, dR0, x_img, y_img,
                           fit.PA, x0_hsc, y0_hsc]
            master_names = ['ra', 'dec', 'n', 'ell', 'dR0', 'x_img_imfit', 
                            'y_img_imfit', 'PA', 'x_hsc_imfit', 'y_hsc_imfit']

            names.extend(master_names)
            data.extend(master_data)

        results = vstack([results, Table(rows=[data], names=names)])

        if save_fit_fig:
            fig, ax = plt.subplots(1, 3, figsize=(16,5))
            fig.subplots_adjust(wspace=0.01)
            img_mod_res(fn,
                        fit.params, 
                        fit_prefix+'_photo_mask.fits',
                        band=band, 
                        show=False, 
                        subplots=(fig, ax),
                        save_fn='fit-'+fit_prefix+'-{}.png'.format(num))
            plt.close('all')

        if clean:
            os.remove(fn)
            os.remove(fit_prefix+'_bestfit_params.txt')
            os.remove(fit_prefix+'_photo_mask.fits')

    if clean and psf_convolve:
        os.remove(psf_fn)

    results = hstack([cat, results])

    return results


def sex_measure(exp, config, apertures, label, add_params, clean, 
                dual_exp=None, sf=None, relpath=''):
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

    param_fn = label+'.params'
    config['PHOT_APERTURES'] = ','.join([str(int(a)) for a in apertures])
    config['PARAMETERS_NAME'] = param_fn
    params = ['MAG_APER('+str(i+1)+')' for i in range(len(apertures))]

    sw = sexpy.SexWrapper(config, params=params, relpath=relpath)
    sw.set_check_images('s', prefix=label+'-')

    #########################################################
    # write exposure for sextractor input and run
    #########################################################

    exp_fn = 'exp-'+label+'.fits'
    exp.writeFits(sw.get_indir(exp_fn))
    cat_label = 'sex-'+label
    run_fn = exp_fn+'[1]'

    if dual_exp is not None:
        meas_band = dual_exp.getFilter().getName().lower()
        dual_label = label+'-'+meas_band
        dual_fn = 'exp-'+dual_label+'.fits'
        dual_exp.writeFits(sw.get_indir(dual_fn))
        cat_label = 'sex-'+dual_label
        run_fn = exp_fn+'[1],'+dual_fn+'[1]'

    sw.run(run_fn, cat=cat_label+'.cat')

    cat = sexpy.read_cat(sw.get_outdir(cat_label+'.cat'))

    if len(cat)>0:

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

        cat.rename_column('MAG_APER', 'MAG_APER_0')
        for i, diam in enumerate(apertures):
            r = utils.pixscale*diam/2 # arcsec
            sb = cat['MAG_APER_'+str(i)] + 2.5*np.log10(np.pi*r**2)
            cat['mu_aper_'+str(i)] = sb

    if add_params and (len(cat)>0):
        x0, y0 = exp.getXY0()
        cat['x_img'] = cat['X_IMAGE'] - 1
        cat['y_img'] = cat['Y_IMAGE'] - 1
        cat['x_hsc'] = cat['X_IMAGE'] + x0 - 1
        cat['y_hsc'] = cat['Y_IMAGE'] + y0 - 1

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

        find_synths = False
        if 'SYNTH' in mask.getMaskPlaneDict().keys() and sf is not None:
            synth_plane = mask.getPlaneBitMask('SYNTH')
            cat['synth_id'] = -1
            cat['synth_offset'] = -1
            find_synths = True
        
        for obj_i in range(len(cat)):
            x, y = cat['x_hsc', 'y_hsc'][obj_i]
            point = afwGeom.Point2I(int(x), int(y))
            for fp_id, fp in enumerate(fpset.getFootprints()):
                if fp.contains(point):
                    hfp = afwDet.HeavyFootprintF(fp, exp.getMaskedImage())
                    pix = hfp.getMaskArray()
                    num_edge = (pix & edge_plane != 0).sum()
                    num_bright = (pix & bright_plane != 0).sum()
                    cat['num_edge_pix'][obj_i] = num_edge
                    cat['num_bright_pix'][obj_i] = num_bright
                    cat['parent_fp_area'][obj_i] = hfp.getArea()
                    cat['fp_id'][obj_i] = fp_id
                    if find_synths:
                        if (pix & synth_plane !=0).sum() > 0:
                            syn_pos = sf.get_psets()[['X0', 'Y0']]
                            x_obj, y_obj = cat['x_img', 'y_img'][obj_i]
                            dist_sq = (syn_pos['X0']-x_obj)**2 +\
                                      (syn_pos['Y0']-y_obj)**2
                            cat['synth_id'][obj_i] = dist_sq.argmin()
                            cat['synth_offset'][obj_i] = np.sqrt(dist_sq.min())
                    break

        for fp_id in np.unique(cat['fp_id']):
            obj_mask = cat['fp_id']==fp_id
            cat['nchild'][obj_mask] = obj_mask.sum()

    #########################################################
    # delete files created by and for sextractor
    #########################################################

    if clean:
        if dual_exp is not None:
            os.remove(sw.get_indir(dual_fn))
        os.remove(sw.get_indir(exp_fn))
        os.remove(sw.get_outdir(cat_label+'.cat'))
        os.remove(sw.get_configdir(param_fn))
        os.remove(sw.get_outdir(label+'-'+'SEGMENTATION.fits'))

    return cat

def _byteswap(arr):
    """
    If array is in big-endian byte order (as astropy.io.fits
    always returns), swap to little-endian for SEP.
    """
    if arr.dtype.byteorder=='>':
        arr = arr.byteswap().newbyteorder()
    return arr


def detection(exp, thresh=0.7, kern_fwhm=1.0, wcs=None, bkg_kws={}, 
              ext_kws={}):
    """
    Source extraction using SEP.
    """

    img = imtools.get_image_ndarray(exp)
    img = _byteswap(img)
    bkg = sep.Background(img, **bkg_kws)
    bkg.subfrom(img)

    kern_sigma = gaussian_fwhm_to_sigma*(kern_fwhm/utils.pixscale)
    kernel = Gaussian2DKernel(kern_sigma)
    kernel.normalize()
    kernel = kernel.array

    sources = sep.extract(
        img, thresh, err=bkg.globalrms, filter_kernel=kernel, **ext_kws)
    sources = Table(sources)

    if wcs is not None:
        x0, y0 = exp.getXY0()
        sources['x_hsc'] = sources['x'] + x0
        sources['y_hsc'] = sources['y'] + y0
        coords = []
        for obj in sources:
            coord = wcs.pixelToSky(obj['x_hsc'], obj['y_hsc'])
            ra, dec = coord.getPosition(afwGeom.degrees)
            coords.append([ra, dec])
        coords = np.array(coords)
        sources['ra'] = coords[:, 0]
        sources['dec'] = coords[:, 1]

    return sources


def photometry(exp, sources, circ_ap_radii=[1.5, 3, 6, 9]):
    """
    Perform aperture photometry with SEP.
    """

    img = imtools.get_image_ndarray(exp)
    img = _byteswap(img)

    # calc kron radius
    x, y = sources['x'], sources['y']
    a, b, theta = sources['a'], sources['b'], sources['theta']
    kronrad, _ = sep.kron_radius(img, x, y, a, b, theta, 6.0)
    
    # flux within 2.5 kron radius
    flux, _, _ = sep.sum_ellipse(img, x, y, a, b, theta, 2.5*kronrad, subpix=1)

    # flux within large ellipse
    flux_ell, _, _ = sep.sum_ellipse(img, x, y, a, b, theta, 6, subpix=1)

    # calc flux radius and magnitude
    flux_radius, flag = sep.flux_radius(
        img,  x, y, 6.*a, 0.5, normflux=flux, subpix=5)

    sources['flux_radius'] = flux_radius
    sources['kron_radius'] = kronrad
    sources['mag_auto'] = 27 - 2.5*np.log10(flux)
    sources['mag_ell'] = 27 - 2.5*np.log10(flux_ell)

    # sum flux in circles of radius=3.0
    for i, r in enumerate(circ_ap_radii):
        flux_circ, _, _ = sep.sum_circle(img, x, y, r)
        mag = 27 - 2.5*np.log10(flux_circ)
        area = np.pi*(r*utils.pixscale)**2
        sources['mag_aper_'+str(i)] = mag
        sources['mu_aper_'+str(i)] = mag + 2.5*np.log10(area)

    return sources
