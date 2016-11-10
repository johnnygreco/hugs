from __future__ import division, print_function

import numpy as np
from scipy.special import gammaincinv, gamma
import scipy.ndimage as ndi
from .sersic import Sersic
from ..utils import pixscale, zpt, check_random_state

__all__ = ['SynthFactory']

PARAM_LIMS = {
    'mu_0': [24, 27],  # mag/arcsec^2
    'r_e': [1.5, 4.0], # arcsec
    'PA': [0, 180],    # degree
    'n': [0.5, 1.1],
    'ell': [0.1, 0.4],
    'X0': [0, 4100],  
    'Y0': [0, 4200]
}

class SynthFactory(object):
    """
    Synthetic (Sersic) galaxy factory.
    """

    def __init__(self, param_lims={}, seed=None):
        """
        Initialize.
        
        Parameters
        ----------
        param_lims : dict, optional
            Parameter limits. 
        seed : int, optional
            Random number generator seed.
        """

        self.param_lims = PARAM_LIMS.copy()
        for k, v in param_lims.items():
            self.param_lims[k] = v
        self.rng = check_random_state(seed)

    def random_params(self, num):
        """
        Generate list of random parameters.

        Parameters
        ----------
        num : int
            Number are random parameters. 
        """

        params_list = []
        for i in range(num):
            params = {}
            for p, lim in self.param_lims.items():
                val = self.rng.uniform(lim[0], lim[1])
                params.update({p:val})
            params_list.append(params)
        return params_list

    def make_galaxy(self, params, img_shape):
        """
        Make synthetic Sersic galaxy.

        Parameters
        ----------
        params : dict
            Sersic parameters. Uses imfit's convention:
            mu_0, r_e, n, X0, Y0, ell, and PA (except for mu_0)
        img_shape : tuple
            Shape of image to embed galaxy in.

        Returns
        -------
        img : ndarray
            Image with synthetic galaxy. 
        """

        # convert mu_0 to I_e 
        mu_0 = params.pop('mu_0')
        b_n = gammaincinv(2.*params['n'], 0.5)
        mu_e = mu_0 + 2.5*b_n/np.log(10)
        I_e = (pixscale**2)*10**((zpt-mu_e)/2.5)
        params['I_e'] = I_e
        params['r_e'] /= pixscale

        # generate image with synth
        model = Sersic(params)
        img = model.array(img_shape)

        return img

    def create_image(self, num_synths, img_shape):
        """
        Inject synthetic galaxies into image.

        Parameters
        ----------
        num_synths : int
            Number of synthetics to inject.
        img_shape : tuple
            Shape of image.

        Returns
        -------
        image : ndarray
            Image of synthetics.
        """

        image = np.zeros(img_shape)
        params_list = self.random_params(num_synths)
        
        for params in params_list:
            galaxy = self.make_galaxy(params, img_shape)
            image += galaxy

        return image
