from __future__ import division, print_function

import numpy as np
import pandas as pd
from scipy.signal import fftconvolve
from scipy.special import gammaincinv
from .sersic import Sersic
from ..utils import pixscale, zpt, check_random_state
from ..utils import embed_slices

__all__ = ['SynthFactory']

PSET_LIMS = {
    'mu_0(i)': [23.5, 27.5], # i-band SB [mag/arcsec^2]
    'r_e': [1.5, 4.0],     # effective radius [arcsec]
    'PA': [0, 180],        # position angle [degree]
    'n': [0.5, 1.1],       # sersic n
    'ell': [0.1, 0.4],     # ellipticity 
    'X0': [0, 4100],       # x position [pixel] 
    'Y0': [0, 4200],       # y position [pixel]
    'g_i': [0.8, 1.2],     # g-i color
    'g_r': [0.2, 1.0]      # g-r color
}

class SynthFactory(object):
    """
    Synthetic (Sersic) galaxy factory.
    """

    def __init__(self, num_synths=None, psets=None,
                 pset_lims={}, random_state=None):
        """
        Initialize.
        
        Parameters
        ----------
        num_synths : int, optional
            Number of synthetics to inject.
        psets : list of dicts, optional
            Set of galaxy parameters.
        pset_lims : dict, optional
            Parameter limits. 
        random_state : int, list of ints, RandomState instance, or None, optional 
            If int or list of ints, random_state is the rng seed.
            If RandomState instance, random_state is the rng.
            If None, the rng is the RandomState instance used by np.random.
	
        random_state : int, optional
            Random number generator random_state.
        """

        self.pset_lims = PSET_LIMS.copy()
        for k, v in pset_lims.items():
            self.pset_lims[k] = v
        self.rng = check_random_state(random_state)
        self._psets = None
        if psets or num_synths:
            self.set_psets(psets=psets, num_synths=num_synths)

    @property
    def num_synths(self):
        return len(self._psets)

    def random_psets(self, num):
        """
        Generate list of random parameters.

        Parameters
        ----------
        num : int
            Number are random parameter sets. 

        Returns
        -------
        psets : pandas.DataFrame
            Set of random parameters
        """

        psets = []
        for i in range(num):
            pset = {}
            for p, lim in self.pset_lims.items():
                val = self.rng.uniform(lim[0], lim[1])
                pset.update({p:val})
            pset['mu_0(g)'] = pset['mu_0(i)'] + pset['g_i']
            pset['mu_0(r)'] = pset['mu_0(i)'] + pset['g_r']
            psets.append(pset)
        psets = pd.DataFrame(psets)
        return psets

    def set_psets(self, psets=None, num_synths=None):
        """
        Set the parameter list. 

        Parameters
        ----------
        psets : list of dicts or pandas.DataFrame, optional
            List of galaxy parameters. If None, random set will be 
            generated according to pset_lims (must give num). 
        num_synths : int, optional
            Number of random param sets to generate. 
        """

        if psets:
            self._psets = pd.DataFrame(psets)
        else:
            self._psets = self.random_psets(num_synths)

    def get_psets(self):
        """
        Return current parameter set list.
        """

        assert self._psets is not None, 'must set galaxy param sets'
        return self._psets

    def write_cat(self, fn, **kwargs):
        """
        Save parameter set list to csv file. 
        """

        assert self._psets is not None, 'must set galaxy param sets'
        self._psets.to_csv(fn, **kwargs)

    def make_galaxy(self, pset, bbox_num_reff=10, band='i', index=None):
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

        # convert mu_0 to I_e and r_e to pixels
        p = pset.copy()
        mu_0 = p['mu_0('+band.lower()+')']
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
        model = Sersic(p, calc_params=True)
        galaxy = model.array(img_shape)

        # add total mag to dataframe
        if index is not None:
            self._psets.set_value(index, 
                                  'm_tot('+band.lower()+')', 
                                  model.m_tot)

        return galaxy 

    def create_image(self, img_shape, band='i', cat_fn=None, **kwargs):
        """
        Generate image of synthetic galaxies.

        Parameters
        ----------
        img_shape : tuple
            Shape of image.
        bands : string, optional
            Photometric band.
        cat_fn : string, optional
            Synth catalog file name. 
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

        assert self._psets is not None, 'must set galaxy param sets'

        image = np.zeros(img_shape)

        for index, pset in self._psets.iterrows():
            galaxy = self.make_galaxy(
                pset, band=band.lower(), index=index, **kwargs)
            gal_pos = np.array([int(pset.Y0), int(pset.X0)])
            img_slice, gal_slice = embed_slices(gal_pos, 
                                                galaxy.shape, 
                                                image.shape)
            image[img_slice] += galaxy[gal_slice]
        if cat_fn:
            self.write_cat(cat_fn)
        return image

    def inject(self, exp, set_mask=True, psf_convolve=True, **kwargs):
        """
        Inject synths into exposure. 

        Parameters
        ----------
        exp : lsst.afw.image.ExposureF 
            Exposure object.
        set_mask : bool, optional
            If True, set SYNTH bit mask plane.
        **kwargs : dict
            Arguments for create_image.
        """

        import lsst.afw.image
        import lsst.afw.geom

        assert self._psets is not None, 'must set galaxy param sets'
        assert type(exp)==lsst.afw.image.ExposureF

        mi = exp.getMaskedImage()
        img = mi.getImage()
        mask = mi.getMask()
        img_arr = img.getArray()
        if set_mask:
            mask.addMaskPlane('SYNTH')
            for index, pset in self._psets.iterrows():
                center = lsst.afw.geom.Point2I(int(pset['X0']), 
                                               int(pset['Y0']))
                bbox = lsst.afw.geom.Box2I(center, center)
                bbox.grow(20)
                bbox.clip(exp.getBBox(lsst.afw.image.LOCAL))
                cutout = mask.Factory(mask, bbox, lsst.afw.image.LOCAL)
                cutout.getArray()[:] += mask.getPlaneBitMask('SYNTH')

        # embed synthetics
        synths =  self.create_image(img_arr.shape, **kwargs)

        if psf_convolve:
            psf = exp.getPsf().computeKernelImage().getArray()
            synths = fftconvolve(synths, psf, 'same')

        img_arr[:] += synths
