"""
Random utilities 
"""
from __future__ import division, print_function

import os
import numpy as np
import yaml
import lsst.afw.image as afwImage
from astropy.table import Table
from astropy.io import fits
from collections import namedtuple

project_dir = os.path.dirname(os.path.dirname(__file__))
default_config_fn = os.path.join(
    project_dir, 'pipe-configs/default_config.yml')
andy_mask_path = '/tigress/HSC/HSC/rerun/goulding/S18A_STARMASK/IMAGE_MASKS'

pixscale = 0.168
zpt = 27.0

#Extinction correction factor for HSC
#A_lambda = Coeff * E(B-V)
ExtCoeff = namedtuple('ExtCoeff', 'g r i z y')
ext_coeff = ExtCoeff(g=3.233, r=2.291, i=1.635, z=1.261, y=1.076)


# Patch metadata
PatchMeta = namedtuple(
    'PatchMeta', 'x0 y0 cleaned_frac small_frac bright_obj_frac good_data_frac')


def read_config(fn=default_config_fn):
    """
    Parse hugs pipeline configuration file.
    """
    with open(fn, 'r') as f:
        config_params = yaml.load(f)
    return config_params


def get_dust_map():
    """
    Use sdfmap to get E(B-V) values from the Schlegel, Finkbeiner & Davis 
    (1998) dust map.

    Notes 
    -----
    sfdmap can be downloaded here http://github.com/kbarbary/sfdmap.
    """
    import sfdmap
    dustmap = sfdmap.SFDMap()
    return dustmap


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
    Get sigma of point-spread function of HSC patch.

    Parameters
    ----------
    exposure : lsst.afw.image.ExposureF
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


def get_exposure(data_id, butler=None, datadir=os.environ.get('HSC_DIR')):
    """
    Return HSC exposure. 

    Parameters
    ----------
    data_id : str, lsst.afw.image.ExposureF, or dict
        Exposure info. Can be file name, an exposure object, 
        or dict of HSC data ID.
    butler : lsst.daf.persistence.Butler
        The Butler.
    datadir : str, optional
        Location of exposure data.

    Returns
    -------
    exp : lsst.afw.image.ExposureF  
        HSC exposure for given data ID.
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
        exp = butler.get('deepCoadd_calexp_hsc', data_id, immediate=True)
    return exp


def mkdir_if_needed(directory):
    """"
    Create directory if it does not exist.
    """
    if not os.path.isdir(directory):
        os.mkdir(directory)


def get_time_label():
    """
    Return time label for new files & directories.
    """
    import time
    label = time.strftime("%Y%m%d-%H%M%S")
    return label


def remove_mask_planes(mask, planes):
    """
    Remove mask planes if they are in mask.

    Parameters
    ----------
    mask : lsst.afw.image.MaskU
	Bit mask.
    planes : list
        Planes to clear
    """
    for plane in planes:
        if plane in list(mask.getMaskPlaneDict().keys()):
            mask.removeAndClearMaskPlane(plane, True)


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


def calc_mask_bit_fracs(exp, planes=['SMALL', 'CLEANED', 'BRIGHT_OBJECT']):
    """
    Calculate the fraction of image flagged as 'cleaned' and 'bright_object'. 

    Parameters
    ----------
    exp : lsst.afw.image.ExposureF

    Returns
    -------
    fracs : dict
        Fraction of pixels with keys cleaned_frac and bright_object_frac.
    """
    mask = exp.getMaskedImage().getMask()
    msk_arr = mask.getArray()
    npix = float(msk_arr.size)
    getBitVal = mask.getPlaneBitMask
    fracs = {}
    for p in planes:
        if p in mask.getMaskPlaneDict().keys():
            npix_p = (msk_arr & getBitVal(p) != 0).sum()
            fracs.update({p.lower()+'_frac': npix_p/npix})
    return fracs


def check_astropy_to_pandas(cat):
    """
    Change astropy table to pandas dataframe if necessary.
    """
    if type(cat)==Table:
        cat = cat.to_pandas()
    return cat


def check_kwargs_defaults(kwargs, defaults):
    """
    Build keyword argument by changing a default set of parameters.

    Parameters
    ----------
    kwargs : dict
        Keyword arguments that are different for default values.
    defaults : dict
        The default parameter values.

    Returns
    -------
    kw : dict
        A new keyword argument.
    """
    kw = defaults.copy()
    for k, v in kwargs.items():
        kw[k] = v
    return kw


def make_noise_image(masked_image, random_state=None):
    """
    Generate Gaussian noise image with zero mean and std equal to the rms level
    of the background in the given image. 

    Parameters
    ----------
    masked_image : lsst.afw.image.MaskedImageF
        Masked image for calculating noise image shape and the scale of the 
        noise fluctuations.
    random_state : 
        Input for check_random_state above.

    Returns 
    -------
    noise_array : ndarray
        The noise image. 
    """
    from .stats import get_clipped_sig_task
    rng = check_random_state(random_state)
    bad_mask_planes = ('EDGE', 'BRIGHT_OBJECT', 'DETECTED', 'SAT')
    task = get_clipped_sig_task(bad_mask_planes=bad_mask_planes)
    stats = task.run(masked_image)
    back_rms = stats.stdev
    dims = masked_image.getDimensions()
    noise_array = back_rms*rng.randn(dims[1], dims[0])
    return noise_array


def solid_angle(ra_lim, dec_lim):
    """
    Calculate solid angle with given ra & dec range.
    All angles in degrees.

    Parameters
    ----------
    ra_lim : list-like, optional
        ra limits.
    dec_lim : list-like, optional
        dec limits.
    Returns
    -------
    area : float
        Solid angle in square degrees.
    """
    ra_lim = np.deg2rad(np.asarray(ra_lim))
    dec_lim = np.deg2rad(np.asarray(dec_lim))
    dsin_dec = np.sin(dec_lim[1]) - np.sin(dec_lim[0])
    area = ra_lim.ptp() * dsin_dec * (180.0/np.pi)**2
    return area


def angsep(ra1, dec1, ra2, dec2, sepunits='arcsec'):
    """
    Return angular separation btwn points

    Inputs
    ------
    ra1 :  float or ndarray
        Right ascension(s) of first point(s) (in degrees). The shape
        must be the same as dec1.
    dec1 : float or ndarray
        Declination(s) of first point(s) (in degrees). The shape must
        be the same as ra1.
    ra2 :  float or ndarray
        Right ascension(s) of second(s) point (in degrees). The shape
        must be the same as dec2.
    dec2 : float or ndarray
        Declination(s) of second point(s) (in degrees). The shape must
        be the same as ra2.
    sepunits : string, optional, default = 'arcsec'
        Angular unit of the returned separation
        (radian, degree, arcsec, or arcmin)

    Returns
    -------
    sep : same type as input angles, angular separation between points

    Note
    ----
    This is function uses the Vincenty formula:
    https://en.wikipedia.org/wiki/Great-circle_distance.
    """

    deg2rad = np.pi/180.0

    ra1 = ra1*deg2rad
    dec1 = dec1*deg2rad
    ra2 = ra2*deg2rad
    dec2 = dec2*deg2rad

    sin_dRA = np.sin(ra2 - ra1)
    cos_dRA = np.cos(ra2 - ra1)
    sin_dec1 = np.sin(dec1)
    sin_dec2 = np.sin(dec2)
    cos_dec1 = np.cos(dec1)
    cos_dec2 = np.cos(dec2)

    num1 = cos_dec2*sin_dRA
    num2 = cos_dec1*sin_dec2 - sin_dec1*cos_dec2*cos_dRA
    denom = sin_dec1*sin_dec2 + cos_dec1*cos_dec2*cos_dRA
    sep = np.arctan2(np.sqrt(num1*num1 + num2*num2), denom)

    conversion = {'radian':1.0, 
                  'arcsec':206264.806, 
                  'arcmin':206264.806/60.0, 
                  'degree':180./np.pi}

    sep *= conversion[sepunits]

    return sep


def fetch_andy_mask(tract, patch, band):
    """
    Get Andy's new bright object mask.
    """
    path = os.path.join(andy_mask_path, 'HSC-' + band.upper())
    fn = os.path.join(path, '{}/{}/mask.fits'.format(tract, patch))
    mask = fits.getdata(fn)
    return mask


# the following functinos were modified from astroML
# https://github.com/astroML/astroML

def ra_dec_to_xyz(ra, dec):
    """
    Convert ra & dec to Euclidean points
    Parameters
    ----------
    ra, dec : ndarrays
        RA and DEC in degrees.
    Returns
    -------
    x, y, z : ndarrays
    """
    sin_ra = np.sin(ra*np.pi/180.)
    cos_ra = np.cos(ra*np.pi/180.)

    sin_dec = np.sin(np.pi/2-dec*np.pi/180.)
    cos_dec = np.cos(np.pi/2-dec*np.pi/180.)

    return (cos_ra*sin_dec, sin_ra*sin_dec, cos_dec)


def angular_dist_to_euclidean_dist(theta, r=1):
    """
    Convert angular distances to Euclidean distances.
    Parameters
    ----------
    theta : float or ndarray
        Angular distance.
    r : float, optional
        Radius. Unit sphere assumed as default.
    Returns
    -------
    d : float or ndarray
        Euclidean distance
    """
    d = 2*r*np.sin(0.5*theta*np.pi/180.)
    return d


def euclidean_dist_to_angular_dist(d, r=1):
    """
    Convert Euclidean distances to angular distances.
    Parameters
    ----------
    d : float or ndarray
        Euclidean distance
    r : float, optional
        Radius. Unit sphere assumed as default.
    Returns
    -------
    theta : float or ndarray
        Angular distance.
    """
    theta = 2*np.arcsin(d/r/2.)*180./np.pi
    return theta
