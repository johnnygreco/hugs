from __future__ import division, print_function

import os
import numpy as np
import lsst.daf.persistence
from lsst.pipe.base import Struct
from astropy.table import Table
from .stats import get_clipped_sig_task
from .synths import inject_synths
from .utils import fetch_andy_mask
from .log import logger

class HugsExposure(object):
    """
    Class to fetch and contain multi-band HSC exposures and some metadata
    associated with the hugs pipeline. 

    Parameters
    ----------
    tract : int 
        HSC tract
    patch : str
        HSC patch
    bands : str, optional
        Photometric bands 
    butler : lsst.daf.persistence.Butler
        HSC data butler
    """

    def __init__(self, tract, patch, bands='gri', butler=None, 
                 coadd_label='deepCoadd_calexp', band_detect='i', 
                 rerun='/tigress/HSC/DR/s18a_wide', use_andy_mask=True):
        self.tract = tract
        self.patch = patch
        self._butler = butler
        self.bands = bands
        self.rerun = rerun
        self.synths = None
        self.fn = {}
        self.stat = {}

        for band in bands:
            data_id = {'tract': tract, 
                       'patch': patch, 
                       'filter': 'HSC-' + band.upper()}
            exp = self.butler.get(coadd_label, data_id, immediate=True)
            fn = self.butler.get(
                coadd_label+'_filename', data_id, immediate=True)[0]
            self.fn[band] = fn
            stat_task = get_clipped_sig_task()
            self.stat[band] = stat_task.run(exp.getMaskedImage())

            if use_andy_mask:
                if band == bands[0]:
                    logger.warn("using andy's bright object masks")
                mask = exp.getMask()
                mask.clearMaskPlane(mask.getMaskPlane('BRIGHT_OBJECT'))
                andy_mask = fetch_andy_mask(tract, patch, band)
                val = mask.getPlaneBitMask('BRIGHT_OBJECT')
                mask.getArray()[andy_mask.astype(bool)] += val

            setattr(self, band.lower(), exp)

        self.x0, self.y0 = self.i.getXY0()
        self.patch_meta = Struct(
            x0 = float(self.x0),
            y0 = float(self.y0),
            good_data_frac = self.good_data_fraction(band_detect),
            small_frac = None,
            cleaned_frac = None,
            bright_obj_frac = None
        )

    def __getitem__(self, attr):
        return self.__getattribute__(attr)

    @property
    def butler(self):
        """
        Let's only load the butler once.
        """
        if self._butler is None:
            self._butler = lsst.daf.persistence.Butler(self.rerun)
        return self._butler

    def get_mask_array(self, band='i', planes=['CLEANED', 'BRIGHT_OBJECT']):
        """
        Get mask from given planes.

        Parameters
        ----------
        band : str
            Photometric band
        planes : list
            The mask planes to include

        Returns
        -------
        mask : ndarray
           Mask with masked pixels = 1 and non-masked pixels = 0 
        """
        mask = self[band].getMaskedImage().getMask()
        arr = np.zeros(mask.getArray().shape, dtype=bool)
        for p in planes: 
            if p in mask.getMaskPlaneDict().keys():
                arr |= mask.getArray() & mask.getPlaneBitMask(p) != 0
        return arr

    def good_data_fraction(self, band='i', 
            bad_masks=['NO_DATA', 'SUSPECT', 'SAT']):
        """
        Find the fraction of pixels that contain data.
        """
        mask = self[band].getMask()
        nodata = np.zeros_like(mask.getArray(), dtype=bool)
        for m in bad_masks:
            nodata |= mask.getArray() & mask.getPlaneBitMask(m) != 0
        nodata = nodata.astype(float)
        no_data_frac = nodata.sum()/nodata.size
        return 1.0 - no_data_frac

    def make_rgb(self, rgb_bands='irg', stretch=0.4, Q=8):
        from astropy.visualization import make_lupton_rgb
        rgb = make_lupton_rgb(self[rgb_bands[0]].getImage().getArray(), 
                              self[rgb_bands[1]].getImage().getArray(), 
                              self[rgb_bands[2]].getImage().getArray(), 
                              stretch=stretch, Q=Q)
        return rgb


class SynthHugsExposure(HugsExposure):

    def __init__(self, synth_cat, tract, patch, bands='gri', butler=None, 
                 coadd_label='deepCoadd_calexp', band_detect='i', 
                 rerun='/tigress/HSC/DR/s18a_wide', use_andy_mask=True, 
                 synth_model='sersic'):

        super(SynthHugsExposure, self).__init__(
            tract, patch, bands, butler, coadd_label, band_detect, rerun, 
            use_andy_mask)

        if  type(synth_cat)==Table:
            self.synths = synth_cat
        else:
            self.synths = synth_cat.get_exp_synths(self[bands[0]])

        msg = 'injecting {} synthetic galaxies into tract={} patch={} exposure'
        logger.warn(msg.format(len(self.synths), tract, patch))
        if len(self.synths) > 0:
            for b in bands:
                inject_synths(self.synths, exp=self[b], band=b, 
                              synth_model=synth_model)
