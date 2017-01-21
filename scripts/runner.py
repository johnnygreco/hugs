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
        seed = [int(time()), tract, p1, p2, p['group_id'], rank]
    else:
        seed = p['seed']

    prefix = os.path.join(p['outdir'], 'hugs-{}-{}'.format(p['tract'],
                                                           p['patch']))
    log_fn = prefix+'.log'
    config = hp.Config(config_fn=p['config_fn'], 
                       log_fn=log_fn, 
                       random_state=seed)
    config.set_patch_id(p['tract'], p['patch'])
    config.logger.info('random seed set to {}'.format(seed))
    
    results = hp.run(config)

    # write source catalog with all objects
    all_det = results.all_detections.to_pandas()
    if len(all_det)>0:
        all_det['tract'] = p['tract']
        all_det['patch'] = p['patch']
    fn = prefix+'-cat-all.csv'
    all_det.to_csv(fn, index=False)

    # write source catalog with all objects
    verified = results.sources.to_pandas()
    if len(verified)>0:
        verified['tract'] = p['tract']
        verified['patch'] = p['patch']
    fn = prefix+'-cat-verified.csv'
    verified.to_csv(fn, index=False)

    # write final source catalog
    candy = results.candy.to_pandas()
    if len(candy)>0:
        candy['tract'] = p['tract']
        candy['patch'] = p['patch']
    fn = prefix+'-cat-candy.csv'
    candy.to_csv(fn, index=False)

    # write mask fractions 
    fn = prefix+'-mask-fracs.csv'
    mask_fracs = pd.DataFrame(results.mask_fracs)
    mask_fracs['tract'] = p['tract']
    mask_fracs['patch'] = p['patch']
    mask_fracs.to_csv(fn, index=False)


def combine_results(outdir):

    all_files = os.listdir(outdir)
    join = os.path.join
    parse = lambda key: [join(outdir, f) for f in all_files if key in f]

    prefix = join(outdir, 'hugs-pipe')
    suffixes = ['-cat-candy.csv',
                '-cat-verified.csv',
                '-cat-all.csv', 
                '-mask-fracs.csv']

    files = [parse(s) for s in suffixes]
     
    for fnames, suffix in zip(files, suffixes):
        df = []
        for fn in fnames:
            try: 
                df_patch = pd.read_csv(fn)
                if len(df_patch>0):
                    df.append(df_patch)
                    os.remove(fn)
                else:
                    os.rename(fn, fn[:-3]+'empty')
            except pd.io.common.EmptyDataError:
                os.rename(fn, fn[:-3]+'empty')
        df = pd.concat(df, ignore_index=True)
        df.to_csv(prefix+suffix, index=False)

    log_fn = prefix+'.log'
    with open(log_fn, 'w') as outfile:
        for fn in parse('log'):
            with open(fn) as infile:
                outfile.write(infile.read())
            os.remove(fn)


def main(pool, patches, config_fn, seed=None):

    patches['seed'] = seed
    patches['config_fn'] = config_fn

    list(pool.map(worker, patches))

    pool.close()

    if len(patches)>1:
        for outdir in np.unique(patches['outdir']):
            combine_results(outdir)
            

if __name__=='__main__':
    from astropy.table import Table
    args = hp.parse_args()

    if args.patches_fn is None:
        assert (args.tract is not None) and (args.patch is not None)
        tract, patch = args.tract, args.patch
        patches = Table([[tract], [patch]], names=['tract', 'patch'])
        outdir = os.path.join(args.outdir, 
                              'solo-run-{}-{}'.format(tract, patch))
        outdir = outdir+'_'+args.label if args.label else outdir
        hp.utils.mkdir_if_needed(outdir)
        patches['outdir'] = outdir
        patches['group_id'] = 0
    else:
        patches = Table.read(args.patches_fn)

    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)
    main(pool, patches, config_fn=args.config_fn, seed=args.seed)
