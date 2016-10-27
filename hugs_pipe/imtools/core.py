from __future__ import division, print_function

import os
import lsst.daf.persistence
import lsst.afw.image as afwImage
import lsst.afw.display as afwDisp
import lsst.afw.geom as afwGeom
import lsst.afw.math as afwMath
hscdir = os.environ.get('HSC_DIR')

__all__ = ['get_cutout', 'smooth_gauss']


def get_cutout(center, size, exp=None, data_id=None, butler=None):
    """
    Generate a cutout of exposure. Must give exposure object or 
    data_id and, optionally, a butler.
     
    Parameters
    ----------
    center : tuple
        Center of desired cutout in the tract
        and patch system.
    size : float
        Size to grow bbox in all directions.
    exp : lsst.afw.image.ExposureF, optional
        Exposure from which to get cutout.    
    data_id : dict, optional
        HSC data ID.
    butler : lsst.daf.persistence.Butler, optional
        the butler.

    Returns
    -------
    cutout : lsst.afw.image.ExposureF
        Desired exposure cutout.
    """

    # get exposure
    if exp is None:
        if butler is None:
            butler = lsst.daf.persistence.Butler(hscdir)
        exp = butler.get('deepCoadd_calexp', data_id, immediate=True)

    # generate bbox and get cutout
    center = afwGeom.Point2I(int(center[0]), int(center[1]))
    bbox = afwGeom.Box2I(center, center)
    bbox.grow(size)
    bbox.clip(exp.getBBox())
    cutout = exp.Factory(exp, bbox, afwImage.PARENT)

    return cutout


def smooth_gauss(masked_image, sigma, nsigma=7.0):
    """
    Smooth image with a Gaussian kernel. 

    Parameters
    ----------
    masked_image : lsst.afw.image.imageLib.MaskedImageF
        Masked image object to be smoothed
    sigma : float
        Standard deviation of Gaussian
    nsigma : float, optional
        Number of sigma for kernel width

    Returns
    -------
    convolved_image : lsst.afw.image.imageLib.MaskedImageF
        The convolved masked image
    """
    width = (int(sigma*nsigma + 0.5) // 2)*2 + 1 # make sure it is odd
    gauss_func = afwMath.GaussianFunction1D(sigma)
    gauss_kern = afwMath.SeparableKernel(width, width, gauss_func, gauss_func)
    convolved_image = masked_image.Factory(masked_image.getBBox())
    afwMath.convolve(convolved_image, masked_image, gauss_kern,
                     afwMath.ConvolutionControl())
    return convolved_image
