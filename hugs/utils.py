from __future__ import division, print_function

import os
import numpy as np
import yaml
import lsst.afw.image as afwImage
from astropy.table import Table
from collections import namedtuple

project_dir = os.path.dirname(os.path.dirname(__file__))
default_config_fn = os.path.join(
    project_dir, 'pipe-configs/default_config.yml')

pixscale = 0.168
zpt = 27.0

#Extinction correction factor for HSC
#A_lambda = Coeff * E(B-V)
ExtCoeff = namedtuple('ExtCoeff', 'g r i z y')
ext_coeff = ExtCoeff(g=3.233, r=2.291, i=1.635, z=1.261, y=1.076)


def read_config(fn):
    with open(fn, 'r') as f:
        config_params = yaml.load(f)
    return config_params


def get_dust_map():
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
        exp = butler.get('deepCoadd_calexp_hsc', data_id, immediate=True)
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


def calc_mask_bit_fracs(exp):
    mask = exp.getMaskedImage().getMask()
    msk_arr = mask.getArray()
    npix = float(msk_arr.size)
    getBitVal = mask.getPlaneBitMask
    fracs = {}
    planes = ['CLEANED', 'BRIGHT_OBJECT']
    for p in planes:
        if p in mask.getMaskPlaneDict().keys():
            npix_p = (msk_arr & getBitVal(p) != 0).sum()
            fracs.update({p.lower()+'_frac': npix_p/npix})
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
