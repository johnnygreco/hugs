
from __future__ import division, print_function

import os
import yaml
import numpy as np
from . import utils

class Config(object):
    """
    
    """

    def __init__(self, param_fn):
        with open(param_fn, 'r') as f:
            params = yaml.load(f)
        if params['data_dir']=='local':
            dir = os.environ.get('LOCAL_DATA')
            self.data_dir = os.path.join(dir, 'hsc')
        elif params['data_dir']=='test':
            self.data_dir = os.environ.get('TEST_DATA_DIR')
        else:
            self.data_dir = params['data_dir']
        self._butler = None

        self.thresh_low = params['thresh_low']
        self.thresh_high = params['thresh_high']
        self.thresh_det = params['thresh_det']
        self.assoc = params['assoc']
        self.deblend_stamps = params['deblend_stamps']

    @property
    def bulter(self):
        """
        """
        if self._butlter is None:
            import lsst.daf.persistence
            self._butler = lsst.daf.persistence.Butler(self.data_dir)
        return self._butler

    def get_exposure(self, data_id):
        """
        """
        if type(data_id)==str:
            import lsst.afw.image
            fn = data_id
            fn = os.path.join(self.data_dir, fn)
            exposure = lsst.afw.image.ExposureF(fn)
        else:
            exposure = self.butler.get('deepCoadd_calexp', 
                                        data_id, immediate=True)
            fn = self.butler.get('deepCoadd_calexp_filename', 
                                 data_id, immediate=True)[0]
        return exposure, fn

    def set_data_id(self, data_id):
        """
        """
        self.exp, self.fn = self.get_exposure(data_id)
        self.wcs = utils.get_astropy_wcs(self.fn)
        self.mi = self.exp.getMaskedImage()
        self.mask = self.mi.getMask()

        self.mask.clearMaskPlane(self.mask.getMaskPlane('DETECTED'))
        if 'DETECTED_NEGATIVE' in list(self.mask.getMaskPlaneDict().keys()):
            self.mask.removeAndClearMaskPlane('DETECTED_NEGATIVE', True)
        self.psf_sigma = utils.get_psf_sigma(self.exp)

        ngrow = self.thresh_low.pop('n_sig_grow')
        self.thresh_low['rgrow'] = int(ngrow*self.psf_sigma + 0.5)

        ngrow = self.thresh_high.pop('n_sig_grow')
        self.thresh_high['rgrow'] = int(ngrow*self.psf_sigma + 0.5)

        ngrow = self.thresh_det.pop('n_sig_grow')
        self.thresh_det['rgrow'] = int(ngrow*self.psf_sigma + 0.5)

        if 'psf sigma' in str(self.assoc['min_pix']):
            nsig = int(self.assoc['min_pix'].split()[0])
            self.assoc['min_pix'] = np.pi*(nsig*self.psf_sigma)**2
