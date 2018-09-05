from __future__ import division, print_function

import os
import numpy as np
import lsst.daf.persistence
from lsst.pipe.base import Struct

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
                 coadd_label='deepCoadd_calexp'):
        self.tract = tract
        self.patch = patch
        self._butler = butler
        self.bands = bands
        self.fn = {}

        for band in bands:
            data_id = {'tract': tract, 
                       'patch': patch, 
                       'filter': 'HSC-'+band.upper()}
            exp = self.butler.get(coadd_label, data_id, immediate=True)
            setattr(self, band.lower(), exp)
            fn = self.butler.get(
                coadd_label+'_filename', data_id, immediate=True)[0]
            self.fn[band] = fn

        self.x0, self.y0 = self.i.getXY0()
        self.patch_meta = Struct(
            x0 = float(self.x0),
            y0 = float(self.y0),
            good_data_frac = self.good_data_fraction(),
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
            hsc_dir = os.environ.get('HSC_DIR')
            self._butler = lsst.daf.persistence.Butler(hsc_dir)
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

    def good_data_fraction(self, band='i'):
        """
        Find the fraction of pixels that contain data.
        """
        mask = self[band].getMask()
        nodata = mask.getArray() & mask.getPlaneBitMask('NO_DATA') != 0
        nodata = nodata.astype(float)
        no_data_frac = nodata.sum()/nodata.size
        return 1.0 - no_data_frac
