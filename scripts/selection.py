#!/usr/env/bin python 
""" 
Make selection cuts on source catalog from hugs_pipe.run, 
and generate coordinate lists for a DAS quarry.
"""

from __future__ import division, print_function

import os
import numpy as np
from astropy.table import Table
from hugs_pipe.utils import pixscale
pipe_io = os.environ.get('HUGS_PIPE_IO')
fn = os.path.join(pipe_io, 'testing/master_testing.cat')

###################################
# selection cut params
###################################
size_3sig = 2.5
r_circ = 1.95
mag_ell = 23
mu_0 = 24

###################################
# quarry params
###################################
size = 30
rerun = 's16_wide'

###################################
# make cuts
###################################

cat = Table.read(fn, format='ascii')
cuts = 3*pixscale*cat['semimajor_axis_sigma'] > size_3sig
cuts &= pixscale*cat['equivalent_radius'] > r_circ
cuts &= cat['mu_3'] > mu_0
cuts &= cat['mag_ell'] < mag_ell
cat = cat[cuts]

###################################
# generate coord lists for quarry
###################################

print('#?   ra       dec        sw     sh   filter '  
      'image mask variance type  rerun')
for source in cat:
    for band in ['G', 'R', 'I']:
        ra, dec = source['ra'], source['dec']
        line = ' {:.7f} {:.7f} {}asec {}asec HSC-{}  true '\
               ' true  true    coadd {}'
        print(line.format(ra, dec, size, size, band, rerun))
