from __future__ import division, print_function

import numpy as np
from astropy.convolution import discretize_model

__all__ = ['exp_kern']

def exp_kern(alpha, size, norm=True, mode='center', factor=10):
    """
    Generate 2D, radially symmetric exponential kernel The kernels are 
    discretized using astropy.convolution.discretize_model.

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
    kern : 2D ndarray
        The convolution kernel. 
    """
    assert size%2!=0, 'ERROR: size must be odd'
    x_range = (-(int(size) - 1) // 2, (int(size) - 1) // 2 + 1)
    model = lambda x, y: np.exp(-np.sqrt(x**2 + y**2)/alpha)
    kern = discretize_model(model, x_range, x_range, mode=mode, factor=factor)
    if norm:
        kern /= kern.sum()
    return kern
