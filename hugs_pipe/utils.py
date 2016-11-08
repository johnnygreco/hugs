from __future__ import division, print_function

import os
import numpy as np
import lsst.afw.math as afwMath

__all__ = ['io', 'pixscale', 'annuli', 'get_astropy_wcs',
           'get_psf_sigma', 'get_test_exp', 'add_cat_params', 
           'get_time_label', 'remove_mask_planes', 'get_fpset']

io = os.environ.get('HUGS_PIPE_IO')
pixscale = 0.168

def annuli(row_c, col_c, r_in, r_out, shape):
    """
    Find the indices within an annulus embedded in 
    an array of the specified shape. 

    Parameters
    ----------
    row_c, col_c : float
        Ceentral coordinates of the annulus.
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
    sources['a_3_sig'] = 3.0*pixscale*sources['semimajor_axis_sigma']
    sources['b_3_sig'] = 3.0*pixscale*sources['semiminor_axis_sigma']
    sources['r_circ_seg'] = pixscale*sources['equivalent_radius']
    q = 1-sources['ellipticity']
    sources['r_circ_ell'] = sources['a_3_sig']*np.sqrt(q)


def remove_mask_planes(mask, planes):
    """
    Reomve mask planes if they are in mask.

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
