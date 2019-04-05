import numpy as np
import skimage
import photutils
from scipy import signal
from astropy.utils import lazyproperty
from photutils import EllipticalAperture
from photutils import CircularAnnulus, CircularAperture, aperture_photometry


__all__ = ['Morphology']


class Morphology(object):
    """
    Class for calculating morphological properties of a single source.
    """
    
    def __init__(self, image, centroid, a, b, ellip, theta, mask=None):
        
        self.image = image.astype(np.float64)
        if mask is None:
            mask = np.zeros_like(image, dtype=bool)
        self.mask = mask
        self.masked_image = np.ma.masked_array(self.image, mask)

        self.centroid = centroid
        self.semimajor_axis = a
        self.semiminor_axis = b
        self.ellipticity = ellip
        self.theta = theta

    def elliptical_mask(self, scale=2):
        ell = EllipticalAperture(
            self.centroid, scale * self.semimajor_axis,
            scale * self.semiminor_axis,  np.deg2rad(self.theta)
        )
        ell_mask = ell.to_mask()[0].to_image(self.image.shape).astype(bool)
        return ell_mask

    def circular_mask(self, radius=4):
        circ =  CircularAperture(self.centroid, radius)
        circ_mask = circ.to_mask()[0].to_image(self.image.shape).astype(bool)
        return circ_mask

    def gini(self, ell_mask_scale=0, circ_mask_radius=0):
        """
        The Gini coefficient.
        References
        ----------
        Abraham, R. G. et al. 2003, ApJ, 588, 218
        Lotz, J. M. et al. 2004, AJ, 128, 163 
        """

        try:
            masked_image = self.masked_image.copy() 
            if ell_mask_scale > 0:
                masked_image.mask |=  ~self.elliptical_mask(ell_mask_scale)
            if circ_mask_radius > 0:
                masked_image.mask |= self.circular_mask(circ_mask_radius)

            data = masked_image.filled(np.nan).flatten()
            data = data[~np.isnan(data)]
            sorted_pix = np.sort(data)

            N = len(sorted_pix)
            norm = 1. / (np.abs(np.mean(sorted_pix)) * N * (N - 1))
            g_coeff = (2 * np.arange(1, N + 1) - N - 1) * np.abs(sorted_pix)
            g_coeff = norm * np.sum(g_coeff)
        except:
            g_coeff = np.nan

        return g_coeff

    def autocorr(self, ell_mask_scale=2, aperture_radius=5, annulus_width=4):
        # Compute 2D autocorrelation function

        try:
            masked_image = self.masked_image.copy() 
            if ell_mask_scale > 0:
                masked_image.mask |=  ~self.elliptical_mask(ell_mask_scale)
            masked_image = masked_image.filled(0.0)

            fft_imgae = np.fft.fft2(masked_image)
            acorr_image = np.fft.ifft2(fft_imgae * np.conjugate(fft_imgae)).real
            acorr_image = np.fft.ifftshift(acorr_image)

            ny, nx = masked_image.shape
            yy, xx = np.mgrid[:ny, :nx]

            circ = CircularAperture([nx // 2, ny // 2], r=aperture_radius)
            ann = CircularAnnulus([nx // 2, ny // 2],
                                  r_in=aperture_radius,
                                  r_out=aperture_radius + annulus_width)

            ann_mean = aperture_photometry(
                acorr_image, ann)['aperture_sum'][0] / ann.area()
            circ_mean = aperture_photometry(
                acorr_image, circ)['aperture_sum'][0] / circ.area()
        except:
            acorr_image = np.nan
            circ_mean = np.nan
            ann_mean = np.nan

        return acorr_image, circ_mean, ann_mean
