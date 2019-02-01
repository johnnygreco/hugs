"""
Class to help creating a source detection pipeline with sep
by creating 'steps' with different configurations. 
"""
import os
import numpy as np
from scipy import ndimage
from astropy.table import Table
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma, gaussian_sigma_to_fwhm
import sep

from . import utils
from .utils import check_kwargs_defaults, check_random_state
from .log import logger

__all__ = ['SepStepperBase', 'SepLsstStepper', 'sep_ellipse_mask']
          

def _byteswap(arr):
    """
    If array is in big-endian byte order (as astropy.io.fits
    always returns), swap to little-endian for SEP.
    """
    if arr.dtype.byteorder=='>':
        arr = arr.byteswap().newbyteorder()
    return arr


def sep_ellipse_mask(sources, image_shape, scale=5.0):

    logger.info('building ellipse mask')

    mask = np.zeros(image_shape, dtype=bool)
    sep.mask_ellipse(mask, sources['x'], sources['y'], sources['a'], 
                     sources['b'], sources['theta'], scale)

    logger.debug('{:.2f}% of patch masked'.format(100 * mask.sum()/mask.size))

    return mask


class SepStepperBase(object):

    SEP_EXTRACT_DEFAULTS = dict(
        thresh=5.0, minarea=10, deblend_nthresh=32,
        deblend_cont=0.005, clean=True, clean_param=1.0, 
        filter_type='conv', filter_num_fwhm=1.0
    )

    SEP_BACK_DEFAULTS = dict(bw=64, bh=64, fw=3, fh=3, fthresh=0.0)

    def __init__(self, config_fn=None, config={}):

        if config_fn is not None:
            config = utils.read_config(config_fn)

        self.step_kws = config
        self.sources = {}
         
    def setup_image(self, image_object, **kwargs):
        """
        Must make the followiing attributes:
            - image
            - noise_image (ndarray or number)
            - psf_fwhm (pixels)
        """
        raise NotImplementedError

    def _measure(self, img, sources, mask=None):

        logger.info('measuring source parameters')

        # calculate "AUTO" parameters
        kronrad, krflag = sep.kron_radius(
            img, sources['x'], sources['y'], sources['a'], sources['b'],
            sources['theta'], 6.0, mask=mask)
        flux, fluxerr, flag = sep.sum_ellipse(
            img, sources['x'], sources['y'], sources['a'], sources['b'],
            sources['theta'], 2.5*kronrad, subpix=1, mask=mask)
        flag |= krflag  # combine flags into 'flag'
        flux[flux<=0] = np.nan
        mag_auto = utils.zpt - 2.5*np.log10(flux)
        r, flag = sep.flux_radius(
            img, sources['x'], sources['y'], 6.*sources['a'], 0.5,
            normflux=flux, subpix=5, mask=mask)

        sources['mag_auto'] = mag_auto
        sources['flux_auto'] = flux
        sources['flux_radius'] = r * utils.pixscale

        # approximate fwhm 
        r_squared = sources['a']**2 + sources['b']**2
        sources['fwhm'] = 2 * np.sqrt(np.log(2) * r_squared) *  utils.pixscale
        
        q = sources['b'] / sources['a']
        area = np.pi * q * sources['flux_radius']**2
        sources['mu_ave_auto'] = sources['mag_auto']  + 2.5 * np.log10(2*area)
        
        area_arcsec = np.pi * (self.psf_fwhm/2)**2 * utils.pixscale**2
        flux, fluxerr, flag = sep.sum_circle(
            img, sources['x'], sources['y'], self.psf_fwhm/2, 
            subpix=5, mask=mask)
        flux[flux<=0] = np.nan
        mu_0 = utils.zpt - 2.5*np.log10(flux / area_arcsec)

        sources['mu_0_aper'] = mu_0

        return sources

    def apply_mask(self, image, mask, inplace=True):
        img = image if inplace else image.copy()
        if type(self.noise_image) == np.ndarray:
            img[mask] = self.noise_image[mask]
        else:
            img[mask] = self.noise_image
        return img 

    def run(self, step_name, mask=None):

        try:
            img = self.image.copy()
        except AttributeError:
            logger.error('You must setup image before running pipeline!')
            exit(1)

        logger.info('running ' + step_name)

        kws = self.step_kws[step_name].copy()

        sep_extract_kws = kws.pop('sep_extract_kws', {})
        sep_extract_kws = check_kwargs_defaults(sep_extract_kws,
                                                self.SEP_EXTRACT_DEFAULTS)

        sep_back_kws = kws.pop('sep_back_kws', {})
        sep_back_kws = check_kwargs_defaults(sep_back_kws,
                                             self.SEP_BACK_DEFAULTS)

        if mask is not None:
            logger.info('applying mask') 
            sep_extract_kws['mask'] = mask
            sep_back_kws['mask'] = mask

        # define gaussian kernel for detection
        # TODO: make filter shape optional
        num_fwhm = sep_extract_kws.pop('filter_num_fwhm', 1)
        kernel_fwhm = self.psf_fwhm * num_fwhm
        logger.info('smoothing with kernel with fwhm = {:.2f} arcsec'.\
                    format(kernel_fwhm * utils.pixscale))
        kern = Gaussian2DKernel(kernel_fwhm * gaussian_fwhm_to_sigma, 
                                mode='oversample')
        kern.normalize()
        sep_extract_kws['filter_kernel'] = kern.array

        # estimate background and subtract from image
        bkg = sep.Background(img, **sep_back_kws)
        img_sub = img - bkg

        # extract sources
        logger.info('detecting with a threshold of {} x background'.\
                    format(sep_extract_kws['thresh']))
        sources, segmap =  sep.extract(
            img_sub, err=bkg.rms(), segmentation_map=True, **sep_extract_kws)

        sources = Table(sources)

        logger.info('found {} sources'.format(len(sources)))

        if kws['do_measure']:
            sources = self._measure(img, sources, mask)

        sources['seg_label'] = np.arange(1, len(sources) + 1)
        self.sources[step_name] = sources

        return sources, segmap


class SepLsstStepper(SepStepperBase):

    def setup_image(self, exposure, random_state=None):
        self.image = exposure.getImage().getArray()
        self.psf_fwhm = utils.get_psf_sigma(exposure) * gaussian_sigma_to_fwhm
        self.noise_image =  utils.make_noise_image(exposure.getMaskedImage(), 
                                                    random_state)
