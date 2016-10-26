#!/usr/env/bin python 

import os
import numpy as np
from astropy.table import Table
import hugs_pipe

pipe_io = os.environ.get('HUGS_PIPE_IO')
fn = os.path.join(pipe_io, 'testing/master_testing.cat')

# selection cuts
size_3sig = 2.5
r_circ = 1.95
mag_ell = 23
mu_0 = 24

cat = Table.read(fn, format='ascii')
cuts = 0.168*3*cat['semimajor_axis_sigma'] > size_3sig
cuts &= 0.168*cat['equivalent_radius'] > r_circ
cuts &= cat['mu_3'] > mu_0
cuts &= cat['mag_ell'] < mag_ell

cat = cat[cuts]
hugs_pipe.viz.view_candy(cat)
