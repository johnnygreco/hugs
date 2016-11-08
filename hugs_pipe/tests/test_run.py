"""
Unit test hugs_pipe.run
"""
from __future__ import division, print_function

import os
import yaml
import lsst.afw.image
from ..run import run
from ..config import Config
from .. import utils

exposure = utils.get_test_exp()
dataDIR = os.environ.get('TEST_DATA_DIR')
fn = os.path.join(dataDIR, 'test_exposure.fits')
dir = os.path.dirname(os.path.realpath(__file__))
dir = os.path.dirname(dir)

def test_run():
    cfg_fn = os.path.join(dir, 'default_config.yml')
    cfg = Config(cfg_fn)
    cfg.set_data_id(fn)
    results = run(cfg, debug_return=True)
    mask = results.exposure.getMaskedImage().getMask()
    planes = list(mask.getMaskPlaneDict().keys())
    assert 'DETECTED' in planes
    assert 'THRESH_LOW' in planes
    assert 'THRESH_HIGH' in planes
    nfp_low = len(results.fpset_low.getFootprints())
    nfp_high = len(results.fpset_high.getFootprints())
    assert nfp_low > nfp_high
