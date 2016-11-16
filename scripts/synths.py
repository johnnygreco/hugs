"""
Run hugs-pipe on data with synthetic UGDs.
"""
from __future__ import division, print_function

import os
import numpy as np
import lsst.afw.image as afwImage
import lsst.afw.display as afwDisp
from lsst.pipe.base import Struct
import hugs_pipe as hp

def calc_mask_bit_fracs(exp):
    mask = exp.getMaskedImage().getMask()
    msk_arr = mask.getArray()
    npix = float(msk_arr.size)
    getBitVal = mask.getPlaneBitMask

    npix_clean = (msk_arr & getBitVal('CLEANED') != 0).sum()
    npix_blend = (msk_arr & getBitVal('BLEND') != 0).sum()
    npix_bright = (msk_arr & getBitVal('BRIGHT_OBJECT') != 0).sum()

    return Struct(clean_frac=npix_clean/npix,
                  bright_frac=npix_bright/npix,
                  blend_frac=npix_blend/npix)

def main(config, pset_lims, num_synths=10, seed=None):

    synths_kwargs = {'num_synths': num_synths,
                     'seed': seed,
                     'pset_lims': pset_lims}

    results = hp.run(config,
                     debug_return=True,
                     inject_synths=True,
                     synths_kwargs=synths_kwargs)

    sources = results.sources
    sf = results.sf
    synths = sf.get_psets()

    recovered = sources[sources['synth_index']>0]
    indices = recovered['synth_index']
    truth = synths.iloc[indices]
    mask_fracs = calc_mask_bit_fracs(results.exp_clean)


if __name__=='__main__':

    seed = None
    tract = 9348
    patch = '8,6'
    #tract = 9617 
    #patch = '7,5'

    pset_lims = {
        'n': [0.7, 0.7],
        'mu0_i': [24, 26],
        'ell': [0.2, 0.2],
        'r_e': [2, 2]
    }

    config = hp.Config()
    data_id = {'tract': tract, 'patch': patch, 'filter': 'HSC-I'}
    config.set_data_id(data_id)

    main(config, pset_lims, seed=seed)
    config.reset_mask_planes()
