from __future__ import division, print_function

import os, sys
import logging
import yaml
import numpy as np
from . import utils
try:
    import coloredlogs
except ImportError:
    pass

class Config(object):
    """
    Class for parsing the hugs_pipe configuration.
    """

    def __init__(self, param_fn=None, data_id=None, log_level='info',
                 log_fn=None):
        """
        Initialization

        Parameters
        ----------
        param_fn : string, optional
            The parameter file name (must be a yaml file).
            If None, will use default config.
        data_id : dict or string, optional
            The HSC calibrated exposure data id (dict) or 
            filename (string). 
        log_level : string, optional
            Level of python logger.
        log_fn : string, optional
            Log file name.
        """

        # read parameter file & setup param dicts
        if param_fn is None:
            dir = os.path.dirname(os.path.realpath(__file__))
            param_fn = os.path.join(dir, 'default_config.yaml')
        with open(param_fn, 'r') as f:
            params = yaml.load(f)
        if params['data_dir']=='hsc':
            self.data_dir = os.environ.get('HSC_DIR')
        else:
            self.data_dir = params['data_dir']
        self.thresh_low = params['thresh_low']
        self.thresh_high = params['thresh_high']
        self.thresh_det = params['thresh_det']
        self.assoc = params['assoc']
        self.deblend_stamps = params['deblend_stamps']

        # set data id if given
        if data_id:
            self.set_data_id(data_id)
        self._butler = None

        # setup logger
        self.logger = logging.getLogger('hugs-pipe logger')
        self.logger.setLevel(getattr(logging, log_level.upper()))
        fmt = '%(name)s: %(asctime)s %(levelname)s: %(message)s'
        try: 
            formatter = coloredlogs.ColoredFormatter(fmt=fmt, 
                                                     datefmt='%m/%d %H:%M:%S')
        except: 
            formatter = logging.Formatter(fmt=fmt, datefmt='%m/%d %H:%M:%S')

        sh = logging.StreamHandler(sys.stdout)
        sh.setFormatter(formatter)
        self.logger.addHandler(sh)

        if log_fn:
            fh = logging.FileHandler(log_fn)
            fh.setFormatter(formatter)
            self.logger.addHandler(fh)

    @property
    def butler(self):
        """
        Let's only load the butler once.
        """
        if self._butler is None:
            import lsst.daf.persistence
            self._butler = lsst.daf.persistence.Butler(self.data_dir)
        return self._butler

    def get_exposure(self, data_id):
        """
        Get the HSC calibrated exposure.

        Parameters
        ----------
        data_id : dict or string
            The HSC calibrated exposure data id (dict) or 
            filename (string).

        Returns
        -------
        exposure : lsst.afw.image.ExposureF
            HSC calibrated exposure.
        fn : string
            The exposure filename.
        """
        if type(data_id)==str:
            import lsst.afw.image
            fn = data_id
            exposure = lsst.afw.image.ExposureF(fn)
        else:
            exposure = self.butler.get('deepCoadd_calexp', 
                                        data_id, immediate=True)
            fn = self.butler.get('deepCoadd_calexp_filename', 
                                 data_id, immediate=True)[0]
        return exposure, fn

    def set_data_id(self, data_id):
        """
        Setup the data id. This must be done before passing a 
        config object to hugs_pipe.run.

        Parameters
        ----------
        data_id : dict or string, optional
            The HSC calibrated exposure data id (dict) or 
            filename (string).
        """

        # get exposure and wcs
        self.exp, self.fn = self.get_exposure(data_id)
        self.wcs = utils.get_astropy_wcs(self.fn)
        self.mi = self.exp.getMaskedImage()
        self.mask = self.mi.getMask()

        # clear detected mask and remove unnecessary plane
        self.mask.clearMaskPlane(self.mask.getMaskPlane('DETECTED'))
        if 'DETECTED_NEGATIVE' in list(self.mask.getMaskPlaneDict().keys()):
            self.mask.removeAndClearMaskPlane('DETECTED_NEGATIVE', True)
        self.psf_sigma = utils.get_psf_sigma(self.exp)

        # setup thresh low/high/det: n_sig_grow --> rgrow
        ngrow = self.thresh_low.pop('n_sig_grow')
        self.thresh_low['rgrow'] = int(ngrow*self.psf_sigma + 0.5)

        ngrow = self.thresh_high.pop('n_sig_grow')
        self.thresh_high['rgrow'] = int(ngrow*self.psf_sigma + 0.5)

        ngrow = self.thresh_det.pop('n_sig_grow')
        self.thresh_det['rgrow'] = int(ngrow*self.psf_sigma + 0.5)

        # convert assoc min_pix to psf units 
        if 'psf sigma' in str(self.assoc['min_pix']):
            nsig = int(self.assoc['min_pix'].split()[0])
            self.assoc['min_pix'] = np.pi*(nsig*self.psf_sigma)**2

        # convert deblend_stamps kernel sigma to psf units
        if 'psf sigma' in self.deblend_stamps['kern_sig_pix']:
            nsig = float(self.deblend_stamps['kern_sig_pix'].split()[0])
            self.deblend_stamps['kern_sig_pix'] = nsig*self.psf_sigma
        
        # set wcs to astropy.wcs.WCS object
        self.deblend_stamps['wcs'] = self.wcs
