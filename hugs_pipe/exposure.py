
from __future__ import division, print_function

import os
import numpy as np
import lsst.daf.persistence

class HugsPipeExposure(object):
    
    def __init__(self, tract, patch, bands='gri', butler=None):
        self.tract = tract
        self.patch = patch
        self._butler = butler
        self.bands = bands
        for band in bands:
            data_id = {'tract': tract, 
                       'patch': patch, 
                       'filter': 'HSC-'+band.upper()}
            exp = self.butler.get('deepCoadd_calexp', data_id, immediate=True)
            setattr(self, band.lower(), exp)

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

    def get_mask_array(self, band='i', 
                       names=['CLEANED', 'BRIGHT_OBJECT', 'BLEND']):
        mask = self[band].getMaskedImage().getMask()
        arr = np.zeros(mask.getArray().shape, dtype=bool)
        for n in names:
            arr |= mask.getArray() & mask.getPlaneBitMask(n) != 0 
        return arr
