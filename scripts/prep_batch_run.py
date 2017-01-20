import os 
from argparse import ArgumentParser
import numpy as np
from astropy.table import Table, vstack
import hugs_pipe as hp

parser = ArgumentParser()
parser.add_argument('groups_fn', type=str)
parser.add_argument('--run_label', type=str, default=None)
args = parser.parse_args()

if args.run_label is None:
    run_label = hp.utils.get_time_label()
else:
    run_label = args.run_label
    
parent_dir = os.path.join(hp.io, 'batch-run-'+run_label)
hp.utils.mkdir_if_needed(parent_dir)

patches = Table()
groups = np.loadtxt(args.groups_fn, dtype=int)

for g in groups:
    group_dir = os.path.join(parent_dir, 'group-'+str(g))
    hp.utils.mkdir_if_needed(group_dir)
    group_patches = hp.get_group_patches(g)
    group_patches['group_id'] = g
    group_patches['outdir'] = group_dir
    patches = vstack([patches, group_patches])

patches.write(os.path.join(parent_dir, 'patches.csv'))
