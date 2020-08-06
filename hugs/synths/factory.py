from __future__ import division, print_function

import numpy as np
import pandas as pd
from scipy.signal import fftconvolve
from scipy.special import gammaincinv
import lsst.afw.image
import lsst.geom
from .sersic import Sersic
from ..utils import pixscale, zpt, check_random_state
from ..utils import embed_slices

try:
    import galsim
    HAS_GALSIM = True
except:
    print('galsim not installed --> no inclined disk synths')
    HAS_GALSIM = False

__all__ = ['inject_synths']


def _make_galaxy(pset, bbox_num_reff=10, band='i'):
    """
    Make synthetic Sersic galaxy.

    Parameters
    ----------
    pset : dict, astropy.table.Row
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
    mu_0 = pset['mu_0_' + band.lower()]
    b_n = gammaincinv(2.*pset['n'], 0.5)
    mu_e = mu_0 + 2.5*b_n/np.log(10)
    I_e = (pixscale**2)*10**((zpt-mu_e)/2.5)
    r_e = pset['r_e'] / pixscale
    
    # calculate image shape
    side = 2*int(bbox_num_reff*r_e) + 1
    img_shape = (side, side)

    params = dict(X0=img_shape[1]//2,
                  Y0=img_shape[0]//2,
                  I_e=I_e, 
                  r_e=r_e, 
                  n=pset['n'], 
                  ell=pset['ell'],
                  PA=pset['PA'])

    # generate image with synth
    model = Sersic(params)
    galaxy = model.array(img_shape)

    return galaxy 


def _make_inclined_disk(src, bbox_num_reff=10, band='i'):
    assert HAS_GALSIM, 'you need galsim to make inclined disks'
    r_pix = src['r_e'] / pixscale
    side = 2*int(bbox_num_reff * r_pix) + 1
    nx, ny = (side, side)
    incl = src['incl'] * galsim.degrees
    model = galsim.InclinedExponential(incl, half_light_radius=src['r_e'], 
                                       scale_h_over_r=src['q0'])
    model = model.rotate(src['PA'] * galsim.degrees)
    flux = 10**(0.4 * (zpt - src[f'm_{band}']))
    galaxy = model.drawImage(nx=nx, ny=ny, scale=pixscale).array * flux
    return galaxy


def inject_synths(cat, exp, bbox_num_reff=10, band='i', psf_convolve=True, 
                  set_mask=True, return_synths=False, synth_model='sersic'):


    image_shape = exp.getDimensions().getY(), exp.getDimensions().getX()
    synth_image = np.zeros(image_shape)

    # make synthetic image
    for src in cat:
        if synth_model == 'sersic':
            galaxy = _make_galaxy(
                src, band=band.lower(), bbox_num_reff=bbox_num_reff)
        elif synth_model == 'disk' or synth_model == 'inclined disk':
            galaxy = _make_inclined_disk(
                src, band=band.lower(), bbox_num_reff=bbox_num_reff)
        gal_pos = np.array([int(src['y']), int(src['x'])])
        img_slice, gal_slice = embed_slices(gal_pos, 
                                            galaxy.shape, 
                                            synth_image.shape)
        synth_image[img_slice] += galaxy[gal_slice]

    if psf_convolve:
        psf = exp.getPsf().computeKernelImage().getArray()
        synth_image = fftconvolve(synth_image, psf, 'same')

    if set_mask:
        mask = exp.getMask()
        mask.addMaskPlane('SYNTH')
        for src in cat:
            center = lsst.geom.Point2I(int(src['x']), 
                                           int(src['y']))
            bbox = lsst.geom.Box2I(center, center)
            bbox.grow(20)
            bbox.clip(exp.getBBox(lsst.afw.image.LOCAL))
            cutout = mask.Factory(mask, bbox, lsst.afw.image.LOCAL)
            cutout.getArray()[:] += mask.getPlaneBitMask('SYNTH')

    exp.getImage().getArray()[:] += synth_image

    if return_synths:
        return synth_image
