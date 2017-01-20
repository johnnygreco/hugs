from __future__ import division, print_function

import os
import numpy as np
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from astropy.table import Table

__all__ = [
    'io', 'pixscale', 'zpt', 'annuli', 'get_astropy_wcs',
    'get_psf_sigma', 'get_test_exp', 'add_cat_params', 
    'get_time_label', 'remove_mask_planes', 'get_fpset', 
    'combine_cats', 'check_random_state', 'embed_slices', 
    'get_group_patches', 'calc_mask_bit_fracs', 'get_exposure',
    'check_kwargs_defaults', 'make_noise_image', 'add_band_to_name'
]

io = os.environ.get('HUGS_PIPE_IO')
local_io = os.environ.get('LOCAL_IO')
pixscale = 0.168
zpt = 27.0


def annuli(row_c, col_c, r_in, r_out, shape):
    """
    Find the indices within an annulus embedded in 
    an array of the specified shape. 

    Parameters
    ----------
    row_c, col_c : float
        Central coordinates of the annulus.
    r_in : float
        Inner radius of annulus. If r_in=0, the 
        annulus becomes a circle. 
    r_out : float
        Outer radius of annulus. 
    shape : tuple
        Shape of the host array.

    Returns
    -------
    row_idx, col_idx : ndarray of int
        Indices of annulus within host array.
    """

    # find bounding box of annulus
    center = np.array([row_c, col_c])
    ul = np.ceil(center - r_out).astype(int)
    lr = np.floor(center + r_out).astype(int)
    ul = np.maximum(ul, np.array([0, 0]))
    lr = np.minimum(lr, np.array(shape[:2]) - 1)
    bb_shape = lr - ul + 1

    # generate the coords within annulus in bbox
    ctr_shift = center - ul
    bb_row, bb_col = np.ogrid[0:float(bb_shape[0]), 0:float(bb_shape[1])]
    radius = np.sqrt((bb_row-ctr_shift[0])**2 + (bb_col-ctr_shift[1])**2)
    row_idx, col_idx = np.nonzero((radius >= r_in) & (radius < r_out))

    # shift back to original coords
    row_idx.flags.writeable = True
    col_idx.flags.writeable = True
    row_idx += ul[0]
    col_idx += ul[1]

    return row_idx, col_idx


def embed_slices(center, arr_shape, img_shape):
    """
    Get slices to embed smaller array into larger image.

    Parameters
    ----------
    center : ndarray
        Center of array in the image coordinates.
    arr_shape : tuple 
        Shape of the array to embed (dimensions must be odd).
    img_shape : tuple
        Shape of the main image array.  

    Returns
    -------
    img_slice, arr_slice : tuples of slices
        Slicing indices. To embed array in image, 
        use the following:
        img[img_slice] = arr[arr_slice]
    """
    arr_shape = np.asarray(arr_shape)
    img_shape = np.asarray(img_shape)

    assert np.alltrue(arr_shape%2 != np.array([0,0]))

    imin = center - arr_shape//2
    imax = center + arr_shape//2 

    amin = (imin < np.array([0,0]))*(-imin)
    amax = arr_shape*(imax<=img_shape-1) +\
           (arr_shape-(imax-(img_shape-1)))*(imax>img_shape-1)

    imin = np.maximum(imin, np.array([0, 0]))
    imax = np.minimum(imax, np.array(img_shape)-1)
    imax += 1

    img_slice = np.s_[imin[0]:imax[0], imin[1]:imax[1]]
    arr_slice = np.s_[amin[0]:amax[0], amin[1]:amax[1]]

    return img_slice, arr_slice

def get_psf_sigma(exposure):
    """
    Get sigma of point-spread function.

    Parameters
    ----------
    exposure : lsst.afw.image.imageLib.ExposureF
        Exposure object. 

    Returns
    -------
    sigma : float
        Approximate standard deviation of the PSF
        of the image in this exposure.
    """
    psf = exposure.getPsf()
    sigma = psf.computeShape().getDeterminantRadius()
    return sigma


def get_astropy_wcs(fn):
    """
    Get an astropy WCS object from fits header.
    """
    from astropy import wcs
    from astropy.io import fits
    header = fits.open(fn)[1].header
    exp_wcs = wcs.WCS(header)
    return exp_wcs


def get_test_exp():
    """
    Return a small test exposure for testing.
    """
    import os
    import lsst.afw.image
    dataDIR = os.environ.get('TEST_DATA_DIR')
    fn = os.path.join(dataDIR, 'test_exposure.fits')
    exposure = lsst.afw.image.ExposureF(fn)
    return exposure


def get_exposure(data_id, butler=None, datadir=os.environ.get('HSC_DIR')):
    """
    Return HSC exposure. 

    Parameters
    ----------
    data_id : str, afwImage.ExposureF, or dict
        Exposure info. Can be file name, an exposure object, 
        or dict of HSC data ID.
    butler : lsst.daf.persistence.Butler
        The Butler.
    datadir : str, optional
        Location of exposure data.
    """
    if type(data_id)==str:
        exp = afwImage.ExposureF(data_id)
    elif type(data_id)==afwImage.ExposureF:
        exp = data_id
    elif type(data_id)==dict:
        for key in ['tract', 'patch', 'filter']:
            assert key in data_id.keys()
        if butler is None:
            import lsst.daf.persistence.Butler
            butler = lsst.daf.persistence.Butler(datadir)
        exp = butler.get('deepCoadd_calexp', data_id, immediate=True)
    return exp


def mkdir_if_needed(directory):
    """"
    Create directory if it does not exist.
    """
    import os
    if not os.path.isdir(directory):
        os.mkdir(directory)


def get_time_label():
    """
    Return time label for new files & directories.
    """
    import time
    label = time.strftime("%Y%m%d-%H%M%S")
    return label


def add_cat_params(sources, tract=None, patch=None):
    """
    Add useful params to final catalog.
    """
    if tract:
        sources['tract'] = [tract]*len(sources)
    if patch:
        sources['patch'] = [patch]*len(sources)
    
    for s in [2, 3]:
        scale = s*pixscale
        sources['a_'+str(s)+'_sig'] = scale*sources['semimajor_axis_sigma']
        sources['b_'+str(s)+'_sig'] = scale*sources['semiminor_axis_sigma']
    sources['r_circ_seg'] = pixscale*sources['equivalent_radius']
    q = 1-sources['ellipticity']
    sources['r_circ_ell'] = sources['a_2_sig']*np.sqrt(q)


def remove_mask_planes(mask, planes):
    """
    Remove mask planes if they are in mask.

    Parameters
    ----------
    mask : lsst.afw.MaskU
	Bit mask.
    planes : list
        Planes to clear
    """
    for plane in planes:
	if plane in list(mask.getMaskPlaneDict().keys()):
	    mask.removeAndClearMaskPlane(plane, True)


def get_fpset(mask, plane):
    """
    Get footprint set associated with mask bit plane.
    """
    import lsst.afw.detection as afwDet
    plane = mask.getPlaneBitMask(plane)
    fpset = afwDet.FootprintSet(
        mask, afwDet.Threshold(plane, afwDet.Threshold.BITMASK))
    return fpset


def combine_cats(catdir, label='master', outdir='same'):
    """
    Combine catalogs from a group of runs.
    """
    from astropy.table import vstack

    files = [fn for fn in os.listdir(catdir) if fn[-3:]=='csv']

    master = Table()
    for fn in files:
        tab = Table.read(os.path.join(catdir, fn), format='ascii')
        master = vstack([master, tab])

    outdir = catdir if outdir=='same' else outdir
    outfile = os.path.join(outdir, label+'.csv')
    master.write(outfile)


def check_random_state(seed):
    """
    Turn seed into a `numpy.random.RandomState` instance.

    Parameters
    ----------
    seed : `None`, int, list of ints, or `numpy.random.RandomState`
        If ``seed`` is `None`, return the `~numpy.random.RandomState`
        singleton used by ``numpy.random``.  If ``seed`` is an `int`,
        return a new `~numpy.random.RandomState` instance seeded with
        ``seed``.  If ``seed`` is already a `~numpy.random.RandomState`,
        return it.  Otherwise raise ``ValueError``.

    Returns
    -------
    random_state : `numpy.random.RandomState`
        RandomState object.

    Notes
    -----
    This routine is adapted from scikit-learn.  See
    http://scikit-learn.org/stable/developers/utilities.html#validation-tools.
    """
    import numbers

    if seed is None or seed is np.random:
        return np.random.mtrand._rand
    if isinstance(seed, (numbers.Integral, np.integer)):
        return np.random.RandomState(seed)
    if isinstance(seed, np.random.RandomState):
        return seed
    if type(seed)==list:
        if type(seed[0])==int:
            return np.random.RandomState(seed)

    raise ValueError('{0!r} cannot be used to seed a numpy.random.RandomState'
                     ' instance'.format(seed))


def get_group_patches(group_id, z_max=0.065, Mh_lims=[12.75, 14.0]):
    """
    Get HSC patches associated with galaxy groups.
    """
    prefix = 'cat_z{}_Mh{}-{}'.format(z_max, Mh_lims[0], Mh_lims[1])
    prefix = os.path.join(local_io, prefix)
    patches_fn = prefix+'_tracts_n_patches.npy'
    patches_dict = np.load(patches_fn).item()
    return Table(patches_dict[group_id]) if group_id else patches_dict


def calc_mask_bit_fracs(exp):
    mask = exp.getMaskedImage().getMask()
    msk_arr = mask.getArray()
    npix = float(msk_arr.size)
    getBitVal = mask.getPlaneBitMask
    fracs = {}
    planes = ['CLEANED', 'BRIGHT_OBJECT', 'BLEND']
    for p in planes:
        if p in mask.getMaskPlaneDict().keys():
            npix_p = (msk_arr & getBitVal(p) != 0).sum()
            fracs.update({p.lower()+'_frac': [npix_p/npix]})
    return fracs


def check_astropy_to_pandas(cat):
    if type(cat)==Table:
        cat = cat.to_pandas()
    return cat


def check_kwargs_defaults(kwargs, defaults):
    kw = defaults.copy()
    for k, v in kwargs.items():
        kw[k] = v
    return kw


def make_noise_image(masked_image, random_state=None):
    from .stats import get_clipped_sig_task
    rng = check_random_state(random_state)
    bad_mask_planes = ('EDGE', 'BRIGHT_OBJECT', 'DETECTED', 'SAT', 'CLIPPED')
    task = get_clipped_sig_task(bad_mask_planes=bad_mask_planes)
    stats = task.run(masked_image)
    back_rms = stats.stdev
    dims = masked_image.getDimensions()
    noise_array = back_rms*rng.randn(dims[1], dims[0])
    return noise_array


def add_band_to_name(cat, band, num_aps):
    params = [
        'MAG_AUTO', 'MAG_ISO', 'MU_THRESHOLD', 'MU_MAX', 'FLUX_RADIUS',
        'FWHM_IMAGE', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE', 'ELLIPTICITY',
        'ISO0', 'ISO3', 'ISO6', 'FLAGS', 'MAG_APER_', 'mu_aper_',
    ]

    for p in params:
        if p[-1]=='_':
            for n in range(num_aps):
                n = str(n)
                cat.rename_column(p+n, p+n+'('+band+')')
        else:
            cat.rename_column(p, p+'('+band+')')
