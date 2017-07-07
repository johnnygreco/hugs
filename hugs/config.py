from __future__ import division, print_function

import os, sys
import logging
import yaml
import numpy as np
import time
import lsst.afw.image
from . import utils
from .exposure import HugsExposure
try:
    import coloredlogs
except ImportError:
    pass

class PipeConfig(object):
    """
    Class for parsing the hugs configuration.
    """

    def __init__(self, config_fn=None, tract=None, patch=None, 
                 log_level='info', log_fn=None, random_state=None, 
                 run_name='hugs'):
        """
        Initialization

        Parameters
        ----------
        config_fn : string, optional
            The parameter file name (must be a yaml file).
            If None, will use default config.
        tract : int, optional
            HSC tract.
        patch : str, optional
            HSC patch
        log_level : string, optional
            Level of python logger.
        log_fn : string, optional
            Log file name.
        random_state : int, list of ints, RandomState instance, or None 
            If int or list of ints, random_state is the rng seed.
            If RandomState instance, random_state is the rng.
            If None, the rng is the RandomState instance used by np.random.
        """

        # read parameter file & setup param dicts
        self.config_fn = config_fn if config_fn else utils.default_config_fn
        with open(self.config_fn, 'r') as f:
            params = yaml.load(f)

        self.data_dir = params['data_dir']
        self.hugs_io = params['hugs_io']

        self._thresh_low = params['thresh_low']
        self._thresh_high = params['thresh_high']
        self._clean = params['clean']
        self.rng = utils.check_random_state(random_state)
        self._clean['random_state'] = self.rng
        self.log_fn = log_fn
        self.log_level = log_level
        self.run_name = run_name
        self._butler = None
        self._timer = None

        self.band_detect = params['band_detect']
        self.band_verify = params['band_verify']
        self.band_meas = params['band_meas']

        # setup for sextractor
        sex_setup = params['sextractor']
        self.sex_config = sex_setup['config']
        self.delete_created_files = sex_setup['delete_created_files']
        self.sex_io_dir = sex_setup['sex_io_dir']
        self.verify_max_sep = sex_setup['verify_max_sep']

        self.circ_aper_radii = self.sex_config['PHOT_APERTURES'].split(',')
        self.circ_aper_radii = [
            utils.pixscale*float(a)/2 for a in self.circ_aper_radii] 

        # set patch id if given
        if tract is not None:
            assert patch is not None
            self.set_patch_id(tract, patch)
        else:
            self.tract = None
            self.patch = None

    def setup_logger(self, tract, patch):
        """
        Setup the python logger.
        """
        name = 'hugs-pipe: {} | {}'.format(tract, patch)
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

    def reset_mask_planes(self):
        """
        Remove mask planes created by hugs-pipe.
        """
        for band in self.bands:
            mask = self.exp[band].getMaskedImage().getMask()
            utils.remove_mask_planes(mask, ['CLEANED', 
                                            'THRESH_HIGH', 
                                            'THRESH_LOW'])

    def set_patch_id(self, tract, patch):
        """
        Setup the tract/patch. This must be done before passing a 
        config object to hugs_pipe.run.

        Parameters
        ----------
        tract : int
            HSC tract.
        patch : str
            HSC patch.
        """
        
        self.tract = tract
        self.patch = patch
        self.setup_logger(tract, patch)

        # careful not to modify parameters
        self.thresh_low = self._thresh_low.copy()
        self.thresh_high = self._thresh_high.copy()
        self.clean = self._clean.copy()

        # get exposure 
        bands = self.band_detect + self.band_verify + self.band_meas
        self.bands = ''.join(set(bands))
        self.exp = HugsExposure(tract, patch, self.bands, self.butler)

        # clear detected mask and remove unnecessary plane
        for band in self.bands:
            mask = self.exp[band].getMaskedImage().getMask()
            utils.remove_mask_planes(mask, ['CR', 
                                            'CROSSTALK',
                                            'DETECTED_NEGATIVE', 
                                            'NOT_DEBLENDED', 
                                            'SUSPECT', 
                                            'UNMASKEDNAN']) 

        try:
            self.psf_sigma = utils.get_psf_sigma(self.exp[self.band_detect])
        except AttributeError:
            self.logger.warning('no psf with exposure... using mean seeing')
            self.psf_sigma = self.mean_seeing_sigma

        # setup thresh low/high/det: n_sig_grow --> rgrow
        ngrow = self.thresh_low.pop('n_sig_grow')
        self.thresh_low['rgrow'] = int(ngrow*self.psf_sigma + 0.5)

        ngrow = self.thresh_high.pop('n_sig_grow')
        self.thresh_high['rgrow'] = int(ngrow*self.psf_sigma + 0.5)

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

        # sextractor parameter file
        fn = '{}-{}-{}-{}.params'.format(
            tract, patch[0], patch[-1], self.run_name)
        self.sex_config['PARAMETERS_NAME'] = os.path.join(self.sex_io_dir, fn)
