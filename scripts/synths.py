"""
Run hugs-pipe on data with synthetic UGDs.
"""
from __future__ import division, print_function

import os
import numpy as np
import hugs_pipe as hp
import lsst.afw.display as afwDisp

tract = 9617
patch = '7,5'
data_id = {'tract': tract, 'patch': patch, 'filter': 'HSC-I'}


synths_kwargs = {
    'num_synths': 10,
    'seed': None, 
    'pset_lims': {'mu0_i': [24., 24.],
                  'r_e': [2.0, 2.0],
                  'n': [1., 1.]}
}

cfg = hp.Config(data_id=data_id)

results = hp.run(cfg, 
                 debug_return=True, 
                 inject_synths=True, 
                 synths_kwargs=synths_kwargs)

hp.viz.display_candies(results.sources, results.exposure)
