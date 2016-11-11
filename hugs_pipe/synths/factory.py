from __future__ import division, print_function

import numpy as np
from scipy.special import gammaincinv, gamma
import scipy.ndimage as ndi
from .sersic import Sersic
from ..utils import pixscale, zpt, check_random_state
from ..utils import embed_slices

__all__ = ['SynthFactory']

PSET_LIMS = {
    'mu0_i': [23.5, 27.5],   # i-band mag/arcsec^2
    'r_e': [1.5, 4.0],   # arcsec
    'PA': [0, 180],      # degree
    'n': [0.5, 1.1],
    'ell': [0.1, 0.4],
    'X0': [0, 4100],  
    'Y0': [0, 4200],
    'g-i': [0.8, 1.2]
}

class SynthFactory(object):
    """
    Synthetic (Sersic) galaxy factory.
    """

    def __init__(self, num_synths=None, pset_list=None,
                 pset_lims={}, seed=None):
        """
        Initialize.
        
        Parameters
        ----------
        num_synths : int, optional
            Number of synthetics to inject.
        pset_list : list of dicts, optional
            Set of galaxy paramters.
        pset_lims : dict, optional
            Parameter limits. 
        seed : int, optional
            Random number generator seed.
        """

        self.param_lims = PSET_LIMS.copy()
        for k, v in pset_lims.items():
            self.pset_lims[k] = v
        self.rng = check_random_state(seed)
        self._pset_list = None
        if pset_list or num_synths:
            self.set_pset_list(pset_list=pset_list, 
                                 num_synths=num_synths)

    def random_pset_list(self, num):
        """
        Generate list of random parameters.

        Parameters
        ----------
        num : int
            Number are random parameter sets. 

        Returns
        -------
        pset_list : list of dicts
            Set of random paramters
        """

        pset_list = []
        for i in range(num):
            pset = {}
            for p, lim in self.param_lims.items():
                val = self.rng.uniform(lim[0], lim[1])
                pset.update({p:val})
            pset['mu0_g'] = pset['mu0_i'] + pset['g-i']
            pset_list.append(pset)
        return pset_list

    def set_pset_list(self, pset_list=None, num_synths=None):
        """
        Set the parameter list. 

        Parameters
        ----------
        pset_list : list of dicts, optional
            List of galaxy parameters. If None, random set will be 
            generated according to param_lims (must give num). 
        num_synths : int, optional
            Number of random param sets to generate. 
        """
        if pset_list:
            self._pset_list = pset_list
        else:
            self._pset_list = self.random_pset_list(num_synths)

    def get_pset_list(self):
        """
        Return current parameter set list.
        """
        return self._pset_list

    def make_galaxy(self, pset, bbox_num_reff=10, band='i'):
        """
        Make synthetic Sersic galaxy.

        Parameters
        ----------
        pset : dict
            Sersic parameters. Uses imfit's convention:
            mu_0, r_e, n, X0, Y0, ell, and PA (except for mu_0)
        bbox_num_reff : int, optional
            Number of r_eff to extend the bounding box.
        band : string, optional
            Photometric band (need for central surface brightness).

        Returns
        -------
        galaxy : ndarray
            Image with synthetic galaxy. 
        """

        # convert mu_0 to I_e 
        p = pset.copy()
        mu_0 = p['mu0_'+band.lower()]
        b_n = gammaincinv(2.*p['n'], 0.5)
        mu_e = mu_0 + 2.5*b_n/np.log(10)
        I_e = (pixscale**2)*10**((zpt-mu_e)/2.5)
        p['I_e'] = I_e
        p['r_e'] /= pixscale
        
        # calculate image shape
        side = 2*int(bbox_num_reff*p['r_e']) + 1
        img_shape = (side, side)
        p['X0'] = img_shape[1]//2
        p['Y0'] = img_shape[0]//2

        # generate image with synth
        model = Sersic(p)
        galaxy = model.array(img_shape)

        return galaxy 

    def create_synth_image(self, img_shape, band='i', **kwargs):
        """
        Inject synthetic galaxies into image.

        Parameters
        ----------
        img_shape : tuple
            Shape of image.
        bands : string, optional
            Photometric band.
        **kwargs : dictm optional
            make_galaxy args. 

        Returns
        -------
        image : ndarrays
            Image of synthetics.

        Notes
        -----
        Must set parameters list before running. 
        """
        assert self._pset_list is not None, 'must set galaxy param set list'

        image = np.zeros(img_shape)

        for pset in self._pset_list:
            galaxy = self.make_galaxy(pset, band=band.lower(), **kwargs)
            gal_pos = np.array([int(pset['Y0']), int(pset['X0'])])
            img_slice, gal_slice = embed_slices(gal_pos, galaxy, image)
            image[img_slice] += galaxy[gal_slice]
        return image

    def inject(self, exp, set_mask=True, **kwargs):
        """
        Inject synths into exposure. 

        Parameters
        ----------
        exp : lsst.afw.image.ExposureF or ndarray
            Exposure object.
        set_mask : bool, optional
            If True, set SYNTH bit mask plane.
        **kwargs : dict
            Arguments for create_image.
        """
        import lsst.afw.image
        import lsst.afw.geom

        assert self._pset_list is not None, 'must set galaxy param set list'

        if type(exp)==lsst.afw.image.imageLib.ExposureF:
            mi = exp.getMaskedImage()
            img = mi.getImage()
            mask = mi.getMask()
            img_arr = img.getArray()
            img_shape = img_arr.shape
            if set_mask:
                mask.addMaskPlane('SYNTH')
                for p in self._pset_list:
                    center = lsst.afw.geom.Point2I(int(p['X0']), int(p['Y0']))
                    bbox = lsst.afw.geom.Box2I(center, center)
                    bbox.grow(20)
                    bbox.clip(exp.getBBox(lsst.afw.image.LOCAL))
                    cutout = mask.Factory(mask, bbox, lsst.afw.image.LOCAL)
                    cutout.getArray()[:] += mask.getPlaneBitMask('SYNTH')
        else:
            img_arr = exp

        # embed synthetics
        synths =  self.create_synth_image(img_shape, **kwargs)
        img_arr[:] += synths
