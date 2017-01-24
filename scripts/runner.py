"""
Run hugs-pipe.
"""
from __future__ import division, print_function

import os
from time import time
import numpy as np
import pandas as pd
import mpi4py.MPI as MPI
import schwimmbad
import hugs_pipe as hp
from hugs_pipe.utils import calc_mask_bit_fracs

def worker(p):

    rank = MPI.COMM_WORLD.Get_rank()
    if p['seed'] is None:
        tract, p1, p2 = p['tract'], int(p['patch'][0]), int(p['patch'][-1])
        seed = [int(time()), tract, p1, p2, rank]
    else:
        seed = p['seed']

    prefix = os.path.join(p['outdir'], 'hugs-pipe')
    log_fn = prefix+'.log'
    config = hp.Config(config_fn=p['config_fn'], 
                       log_fn=log_fn, 
                       random_state=seed)
    config.set_patch_id(p['tract'], p['patch'])
    config.logger.info('random seed set to {}'.format(seed))
    
    results = hp.run(config)

    # write source catalog with all objects
    all_det = results.all_detections.to_pandas()
    ext = 'csv' if len(all_det)>0 else 'empty'
    fn = prefix+'-cat-all.'+ext
    all_det.to_csv(fn, index=False)

    # write source catalog with all objects
    verified = results.sources.to_pandas()
    ext = 'csv' if len(verified)>0 else 'empty'
    fn = prefix+'-cat-verified.'+ext
    verified.to_csv(fn, index=False)

    # write final source catalog
    candy = results.candy.to_pandas()
    ext = 'csv' if len(candy)>0 else 'empty'
    fn = prefix+'-cat-candy.'+ext
    candy.to_csv(fn, index=False)

    # write mask fractions 
    fn = prefix+'-mask-fracs.csv'
    mask_fracs = pd.DataFrame(results.mask_fracs)
    mask_fracs.to_csv(fn, index=False)


if __name__=='__main__':
    from astropy.table import Table
    args = hp.parse_args()

    if args.tract is not None:
        assert args.patch is not None
        tract, patch = args.tract, args.patch
        patches = Table([[tract], [patch]], names=['tract', 'patch'])
        outdir = os.path.join(args.outdir, 
                              'solo-run-{}-{}'.format(tract, patch))
        outdir = outdir+'_'+args.label if args.label else outdir
        hp.utils.mkdir_if_needed(outdir)
        patches['outdir'] = outdir
    else:
        assert args.patches_fn is not None
        rank = MPI.COMM_WORLD.Get_rank()
        patches = Table.read(args.patches_fn)
        if rank==0:
            outdirs = []
            for tract, patch in patches['tract', 'patch']:
                path = os.path.join(args.outdir, str(tract))
                hp.utils.mkdir_if_needed(path)
                path = os.path.join(path, patch)
                hp.utils.mkdir_if_needed(path)
                outdirs.append(path)
            patches['outdir'] = outdirs

    patches['seed'] = args.seed
    patches['config_fn'] = args.config_fn

    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)

    list(pool.map(worker, patches))
    pool.close()
