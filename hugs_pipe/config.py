from __future__ import division, print_function

import os, sys
import logging
import yaml
import numpy as np
import time
import lsst.afw.image
from . import utils
try:
    import coloredlogs
except ImportError:
    pass

class Config(object):
    """
    Class for parsing the hugs_pipe configuration.
    """

    def __init__(self, config_fn=None, data_id=None, log_level='info',
                 log_fn=None, random_state=None):
        """
        Initialization

        Parameters
        ----------
        config_fn : string, optional
            The parameter file name (must be a yaml file).
            If None, will use default config.
        data_id : dict or string, optional
            The HSC calibrated exposure data id (dict) or 
            filename (string). 
        log_level : string, optional
            Level of python logger.
        log_fn : string, optional
            Log file name.
        random_state : int, list of ints, RandomState instance, or None, optional 
            If int or list of ints, random_state is the rng seed.
            If RandomState instance, random_state is the rng.
            If None, the rng is the RandomState instance used by np.random.
        """

        # read parameter file & setup param dicts
        if config_fn is None:
            filedir = os.path.dirname(os.path.realpath(__file__))
            self.config_fn = os.path.join(filedir, 'default_config.yml')
        else:
            self.config_fn = config_fn
        with open(self.config_fn, 'r') as f:
            params = yaml.load(f)
        if params['data_dir']=='hsc':
            self.data_dir = os.environ.get('HSC_DIR')
        else:
            self.data_dir = params['data_dir']
        self._thresh_low = params['thresh_low']
        self._thresh_high = params['thresh_high']
        self._thresh_det = params['thresh_det']
        self._find_blends = params['find_blends']
        self._measure_sources = params['measure_sources']
        self._clean = params['clean']
        self._photometry = params['photometry']
        self._butler = None
        self._timer = None
        self.log_fn = log_fn
        self.log_level = log_level
        self.rng = utils.check_random_state(random_state)
        self._clean['random_state'] = self.rng

        self.phot_colors = params['phot_colors']

        # set data id if given
        if data_id:
            self.set_data_id(data_id)
        else:
            self.data_id = None

    def setup_logger(self, data_id):
        """
        Setup the python logger.
        """
        if type(data_id)==str:
            name = data_id
        elif type(data_id)==dict:
            t, p, b = data_id['tract'], data_id['patch'], data_id['filter']
            name = 'hugs-pipe: {} | {} | {}'.format(t, p, b)
        else:
            name = 'hugs-pipe'
        self.logger = logging.getLogger(name)
        self.logger.setLevel(getattr(logging, self.log_level.upper()))
        fmt = '%(name)s: %(asctime)s %(levelname)s: %(message)s'
        try: 
            formatter = coloredlogs.ColoredFormatter(fmt=fmt, 
                                                     datefmt='%m/%d %H:%M:%S')
        except: 
            formatter = logging.Formatter(fmt=fmt, datefmt='%m/%d %H:%M:%S')

        if not self.logger.handlers:
            sh = logging.StreamHandler(sys.stdout)
            sh.setFormatter(formatter)
            self.logger.addHandler(sh)

            if self.log_fn:
                fh = logging.FileHandler(self.log_fn)
                fh.setFormatter(formatter)
                self.logger.addHandler(fh)

            msg = 'starting hugs-pipe with config file '+self.config_fn
            self.logger.info(msg)

    @property
    def butler(self):
        """
        Let's only load the butler once.
        """
        if self._butler is None:
            import lsst.daf.persistence
            self._butler = lsst.daf.persistence.Butler(self.data_dir)
        return self._butler

    @property
    def timer(self):
        """
        Timer for pipeline. 
        """
        if self._timer is None:
            self._timer = time.time()
        else:
            delta_time = (time.time() - self._timer)/60.0
            self._timer = time.time()
            return delta_time

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
            fn = data_id
            exposure = lsst.afw.image.ExposureF(fn)
        else:
            exposure = self.butler.get('deepCoadd_calexp', 
                                        data_id, immediate=True)
            fn = self.butler.get('deepCoadd_calexp_filename', 
                                 data_id, immediate=True)[0]
        return exposure, fn

    def reset_mask_planes(self):
        """
        Remove mask planes created by hugs-pipe.
        """
        utils.remove_mask_planes(self.mask, ['BLEND', 
                                             'CLEANED', 
                                             'SYNTH',
                                             'THRESH_HIGH', 
                                             'THRESH_LOW'])

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
        
        self.data_id = data_id
        self.setup_logger(data_id)

        # careful not to modify parameters
        self.thresh_low = self._thresh_low.copy()
        self.thresh_high = self._thresh_high.copy()
        self.thresh_det = self._thresh_det.copy()
        self.clean = self._clean.copy()
        self.find_blends = self._find_blends.copy()
        self.measure_sources = self._measure_sources.copy()
        self.photometry = self._photometry.copy()

        # get exposure 
        if type(data_id)==lsst.afw.image.imageLib.ExposureF:
            self.exp = data_id
            self.fn = None
        else:
            self.exp, self.fn = self.get_exposure(data_id)
        self.mi = self.exp.getMaskedImage()
        self.mask = self.mi.getMask()

        # get color data for forced photometry
        if self.phot_colors:
            self.color_data = {}
            for color in self.phot_colors:
                if color.upper() != 'I':
                    if type(data_id)==str:
                        _id = self.fn.replace('HSC-I', 'HSC-'+color.upper())
                    else:
                        _id = data_id.copy()
                        _id['filter'] = 'HSC-'+color.upper()
                    _exp, _ = self.get_exposure(_id)
                    _img = _exp.getMaskedImage().getImage().getArray() 
                    self.color_data.update({color.upper(): _img})

        # clear detected mask and remove unnecessary plane
        self.mask.clearMaskPlane(self.mask.getMaskPlane('DETECTED'))
        utils.remove_mask_planes(self.mask, ['CR', 
                                             'CROSSTALK',
                                             'DETECTED_NEGATIVE', 
                                             'NOT_DEBLENDED']) 

        try:
            self.psf_sigma = utils.get_psf_sigma(self.exp)
        except AttributeError:
            self.logger.warning('no psf with exposure... using mean seeing')
            self.psf_sigma = self.mean_seeing_sigma

        # setup thresh low/high/det: n_sig_grow --> rgrow
        ngrow = self.thresh_low.pop('n_sig_grow')
        self.thresh_low['rgrow'] = int(ngrow*self.psf_sigma + 0.5)

        ngrow = self.thresh_high.pop('n_sig_grow')
        self.thresh_high['rgrow'] = int(ngrow*self.psf_sigma + 0.5)

        ngrow = self.thresh_det.pop('n_sig_grow')
        self.thresh_det['rgrow'] = int(ngrow*self.psf_sigma + 0.5)

        # convert clean min_pix to psf units 
        if 'psf sigma' in str(self.clean['min_pix_low_thresh']):
            nsig = int(self.clean['min_pix_low_thresh'].split()[0])
            self.clean['min_pix_low_thresh'] = np.pi*(nsig*self.psf_sigma)**2

        # clean n_sig_grow --> rgrow
        ngrow = self.clean.pop('n_sig_grow')
        if ngrow is None:
            self.clean['rgrow'] = None
        else:
            self.clean['rgrow'] = int(ngrow*self.psf_sigma + 0.5)

        # convert measure_sources kernel sigma to psf units
        if 'psf sigma' in self.measure_sources['kern_sig_pix']:
            nsig = float(self.measure_sources['kern_sig_pix'].split()[0])
            self.measure_sources['kern_sig_pix'] = nsig*self.psf_sigma
