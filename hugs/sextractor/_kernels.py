from __future__ import division, print_function

import os
import numpy as np

__all__ = ['exp_kern', 'gauss_kern', 'make_conv_file']

def make_conv_file(kern_name, width_param, size, outdir='default', 
                   norm=True, fmt='%.8f', **kwargs):
    """
    Make a convolution file for given kernel.

    Parameters
    ----------
    kern_name : string
        The kernel name (currently only gauss or exp).
    width_param : float
        The param that describes the width of the kernel in
        pixels (fwhm for gauss and scale length for exp).
    size : odd int
        Number of pixel in x & y directions.
    outdir : string, optional
        The directory to write file in. If 'default', 
        will use the SEX_CONFIG_DIR env variable.
    norm : bool, optional
        If True, sextractor will normalize the kernel 
    fmt : string, optional
        Numpy savetxt output format.
    **kwargs : dict, optional
        Optional parameters for kernel functions.
    """
    assert (type(size) is int) and (size%2!=0)

    if outdir=='default':
        from ._settings import SEX_CONFIG_DIR
        outdir = SEX_CONFIG_DIR

    comment = '# {0}x{0} convolution mask '.format(size)

    if kern_name=='exp':
        kern = exp_kern(width_param, size, **kwargs)
        comment += 'of an exponential with alpha = {} pixels'.format(width_param)
    elif kern_name=='gauss':
        kern = gauss_kern(width_param, size, **kwargs)
        comment += 'of a Gaussian with fwhm = {} pixels'.format(width_param)
    else:
        raise NameError(kern_name+' is not an available kernel')

    convfile = kern_name+'_{0:.1f}_{1}x{1}.conv'.format(width_param, size)
    header = 'CONV NORM' if norm else 'CONV NONORM'
    header += '\n'+comment
    outfile = os.path.join(outdir, convfile)
    np.savetxt(outfile, kern, header=header, comments='', fmt=fmt)

def exp_kern(alpha, size, norm_array=False, mode='center', factor=10):
    """
    Generate 2D, radially symmetric exponential kernal for sextractor. 
    The kernels are discretized using `astropy.convolution.discretize_model`.

    Parameters
    ----------
    alpha : float
        The scale length of the exponential.
    size : odd int
        Number of pixel in x & y directions.
    norm_array : bool, optional
        If True, normalize the kern array. 
    mode : str, optional
        One of the following discretization modes: 
        'center', 'oversample', 'linear_interp', 
        or 'integrate'. See astropy docs for details. 
    factor : float or int
        Factor of oversampling. 

    Returns
    -------
    kern : 2D ndarray, if (return_fn=False)
        The convolution kernel. 
    kern, fn, comment : ndarray, string, string (if return_fn=True)
        The kernel, file name, and the comment
        for the file.

    Notes
    -----
    The kernel will be normalized by sextractor by 
    default, and the non-normalized files are a bit 
    cleaner due to less leading zeros. 
    """
    assert size%2!=0, 'ERROR: size must be odd'
    from astropy.convolution import discretize_model

    x_range = (-(int(size) - 1) // 2, (int(size) - 1) // 2 + 1)
    model = lambda x, y: np.exp(-np.sqrt(x**2 + y**2)/alpha)
    kern = discretize_model(model, x_range, x_range, mode=mode, factor=factor)

    if norm_array:
        kern /= kern.sum()
    
    return kern

def gauss_kern(width, size, width_type='fwhm', norm_array=False, mode='center', 
               factor=10):
    """
    Generate a 2D radially symmetric Gaussian kernel for sextractor,
    The kernels are discretized using `astropy.convolution.discretize_model`.

    Parameters
    ----------
    width : float
        The width of the Gaussian. 
    size : odd int
        Number of pixel in x & y directions.
    width_type : string, optional
        The width type given ('fwhm' or 'sigma').
    norm_array : bool, optional
        If True, normalize the kern array. 
    mode : str, optional
        One of the following discretization modes: 
        'center', 'oversample', 'linear_interp', 
        or 'integrate'. See astropy docs for details. 
    factor : float or int
        Factor of oversampling. 

    Returns
    -------
    kern : 2D ndarray, if (return_fn=False)
        The convolution kernel. 
    kern, fn, comment : ndarray, string, string (if return_fn=True)
        The kernel, file name, and the comment
        for the file.

    Notes
    -----
     i) The kernel will be normalized by sextractor by 
        default, and the non-normalized files are a bit 
        cleaner due to less leading zeros. 
    """
    assert size%2!=0, 'ERROR: size must be odd'
    from astropy.convolution import discretize_model

    convert = 2.0*np.sqrt(2*np.log(2))
    if width_type=='fwhm':
        sigma = width/convert
        fwhm = width
    elif width_type=='sigma':
        sigma = width
        fwhm = width*convert

    x_range = (-(int(size) - 1) // 2, (int(size) - 1) // 2 + 1)
    model = lambda x, y: np.exp(-(x**2 + y**2)/(2*sigma**2))
    kern = discretize_model(model, x_range, x_range, mode=mode, factor=factor)

    if norm_array:
        kern /= kern.sum()

    return kern
