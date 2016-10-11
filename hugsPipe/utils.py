from __future__ import division, print_function

import numpy as np
import lsst.afw.math as afwMath

__all__ = ['pixscale', 'annuli', 'smooth_gauss', 'get_psf_sigma', 'viz']

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


def viz(image, transparency=75, frame=1, 
        colors=[('THRESH_HIGH', 'magenta'), ('THRESH_LOW', 'yellow')]):
    """
    Visualize retsults with ds9.
    """
    import lsst.afw.display as afwDisplay
    disp = afwDisplay.Display(frame)
    for name, color in colors:
        disp.setMaskPlaneColor(name, color)
    disp.setMaskTransparency(transparency)
    disp.mtv(image)
    return disp
