from __future__ import division, print_function

import os
import numpy as np
from . import utils
from . import primitives as prim


__all__ = ['run']

############################################################
# Default parameters for thresholding and association
############################################################

HSC_DIR = os.environ.get('HSC_DIR')
THRESH = {'high': 15.0, 'low':3.0, 'det':3.5}
NPIX = {'high': 1, 'low':50, 'det':100}
ASSOC = {'r_in': 5, 'r_out': 15, 'max_on_bit': 10, 
        'min_pix': '6 psf sigma', 'dilate': None, 
        'plane_name': 'THRESH_HIGH'}
DEFAULTS = [THRESH, NPIX, ASSOC]


def _get_params(dataID, thresh, npix, assoc, butler, data_dir):
    """
    Get function parameters.
    """
    if type(dataID)==str:
        import lsst.afw.image
        exposure = lsst.afw.image.ExposureF(dataID)
    else:
        if butler is None:
            import lsst.daf.persistence
            butler = lsst.daf.persistence.Butler(data_dir)
        exposure = butler.get('deepCoadd_calexp', dataID, immediate=True)
    params = [thresh, npix, assoc]
    for i in range(len(params)):
        for k,v in list(params[i].items()):
            assert k in list(DEFAULTS[i].keys()), 'invalid parameter: '+k
            DEFAULTS[i][k] = v
        params[i] = DEFAULTS[i].copy()
    params.append(exposure)
    return params


def run(dataID, thresh={}, npix={}, assoc={}, butler=None, 
        kern_fwhm=3.5, data_dir=HSC_DIR, visualize=True):
    """
    Run pipeline.

    Parameters
    ----------

    Returns
    -------
    """

    ############################################################
    # Get parameters and setup the mask planes. We don't need
    # the negative detection mask.
    ############################################################

    thresh, npix, assoc, exposure = _get_params(dataID, thresh, npix, 
                                                assoc, butler, data_dir)
    mi = exposure.getMaskedImage()
    mask = mi.getMask()
    mask.clearMaskPlane(mask.getMaskPlane('DETECTED'))
    mask.removeAndClearMaskPlane('DETECTED_NEGATIVE', True) 

    ############################################################
    # Smooth image at psf scale.
    ############################################################

    psf_sigma = utils.get_psf_sigma(exposure)
    mi_smooth = utils.smooth_gauss(mi, psf_sigma)
    rgrow = int(2.4*psf_sigma + 0.5)

    ############################################################
    # Image thesholding at low and high thresholds.
    ############################################################

    fp_low = prim.image_threshold(
        mi_smooth, thresh['low'], mask=mask,
        plane_name='THRESH_LOW', rgrow=rgrow, npix=npix['low'])
    fp_high = prim.image_threshold(
        mi_smooth, thresh['high'], mask=mask,
        plane_name='THRESH_HIGH', rgrow=rgrow, npix=npix['high'])

    ############################################################
    # Generate noise array, which we will use to replace 
    # unwanted sources with noise. 
    ############################################################

    shape = mask.getArray().shape
    back_rms = mi.getImage().getArray()[mask.getArray()==0].std()
    noise_array = back_rms*np.random.randn(shape[0], shape[1])

    ############################################################
    # Get "association" segmentation image; will be non-zero for 
    # low-thresh footprints that are associated with high-thresh
    # footprints. Then, replace these sources with noise.
    ############################################################

    if 'psf sigma' in str(assoc['min_pix']):
        nsig = int(assoc['min_pix'].split()[0])
        assoc['min_pix'] = np.pi*(nsig*psf_sigma)**2
    assoc = prim.associate(mask, fp_low, **assoc)
        
    exp_clone = exposure.clone()
    mi_clone = exp_clone.getMaskedImage()
    mi_clone.getImage().getArray()[assoc!=0] = noise_array[assoc!=0]
    if visualize:
        displays = []
        displays.append(utils.viz(exp_clone, 75, 1)) 

    ############################################################
    # Smooth with large kernel for detection.
    ############################################################

    fwhm = kern_fwhm/utils.pixscale # pixels
    sigma = fwhm/(2*np.sqrt(2*np.log(2)))
    mi_clone_smooth = utils.smooth_gauss(mi_clone, sigma)

    ############################################################
    # Image thresholding at medium threshold for detection.
    ############################################################

    fp_det = prim.image_threshold(
        mi_clone_smooth, thresh['det'], mask=mask,
        plane_name='DETECTED', npix=npix['det'])

    ############################################################
    # Deblend sources in footprints
    ############################################################
    


    return_objects = [fp_low, fp_high, fp_det, exposure, 
                      exp_clone, mi_clone_smooth]
    if visualize:
       displays.append(utils.viz(exposure, 60, 2))
       return_objects.append(displays)

    return return_objects
