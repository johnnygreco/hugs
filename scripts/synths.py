"""
Run hugs-pipe on data with synthetic UGDs.
"""
from __future__ import division, print_function

import os
import numpy as np
import hugs_pipe as hp
import lsst.afw.display as afwDisp
synthdir = os.path.join(hp.io, 'synthetic_udgs')
outdir = os.path.join(hp.io, 'synth_results')

# select data at random to run 
num_run = 1
files = np.array([f for f in os.listdir(synthdir)\
                  if f[-4:]=='fits' if f.split('-')[2]=='I'])
files = files[np.random.randint(files.size, size=num_run)]
#####

log_fn = os.path.join(outdir, 'synths.log')
config_fn = os.path.join(outdir, 'config.yaml')
config = hp.Config(config_fn=config_fn, log_fn=log_fn)

files = ['calexp-HSC-I-8766-4,3.fits']


for fn in files:
    tract, patch = int(fn.split('-')[3]), fn.split('-')[4][:3]
    fn = os.path.join(synthdir, fn)
    prefix = os.path.join(outdir, 'synths-{}-{}'.format(tract, patch))
    config.set_data_id(fn)
    results = hp.run(config, debug_return=True)
    sources = results.sources
    sources.write(prefix+'.csv')

    exp = results.exposure
    exp_clean = results.exp_clean
    mask = exp.getMaskedImage().getMask()	
    mask.removeAndClearMaskPlane('BAD', True)

    disp = hp.viz.display_candies(sources, exp)
    disp2 = afwDisp.Display(2)
    disp2.setMaskTransparency(75)
    disp2.setMaskPlaneColor('THRESH_HIGH', 'magenta')
    disp2.setMaskPlaneColor('THRESH_LOW', 'yellow')
    disp2.setMaskPlaneColor('FAKE', 'green')
    disp2.mtv(exp_clean)

    mask_clean = exp_clean.getMaskedImage().getMask()	
    mask_clean.removeAndClearMaskPlane('THRESH_HIGH', True)
    exp_clean.writeFits('/Users/protostar/Desktop/wvtest.fits')
