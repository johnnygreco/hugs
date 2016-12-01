"""
Run hugs-pipe.
"""
from __future__ import division, print_function

import os
from time import time
import multiprocessing
import numpy as np
import pandas as pd
import schwimmbad
import hugs_pipe as hp
from hugs_pipe.utils import calc_mask_bit_fracs

def worker(p):

    data_id = {'tract': p['tract'], 'patch': p['patch'], 'filter': 'HSC-I'}
    prefix = os.path.join(p['outdir'], 'hugs-{}-{}'.format(p['tract'],
                                                           p['patch']))
    log_fn = prefix+'.log'

    if p['seed'] is None:
        tract, p1, p2 = p['tract'], int(p['patch'][0]), int(p['patch'][-1])
        pid = multiprocessing.current_process().pid
        seed = [int(time())+pid, tract, p1, p2, pid]
    else:
        seed = p['seed']

    config = hp.Config(config_fn=p['config_fn'], 
                       log_fn=log_fn, 
                       random_state=seed)
    config.set_data_id(data_id)
    config.logger.info('random seed set to {}'.format(seed))
    
    results = hp.run(config)

    # write source catalog
    sources = results.sources.to_pandas()
    fn = prefix+'-cat.csv'
    sources.to_csv(fn, index=False)

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

    files = [parse('cat'), parse('mask-fracs')]

    prefix = join(outdir, 'hugs-pipe')
    suffixes = ['-cat.csv', '-mask-fracs.csv']
     
    for fnames, suffix in zip(files, suffixes):
        df = []
        for fn in fnames:
            try: 
                df.append(pd.read_csv(fn))
                os.remove(fn)
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


def main(pool, patches, outdir, config_fn, seed=None):

    patches['outdir'] = outdir
    patches['seed'] = seed
    patches['config_fn'] = config_fn

    pool.map(worker, patches)
    pool.close()

    if len(patches)>1:
        combine_results(outdir)


if __name__=='__main__':
    args = hp.parse_args()

    if args.group_id is None:
        from astropy.table import Table
        assert (args.tract is not None) and (args.patch is not None)
        tract, patch = args.tract, args.patch
        patches = Table([[tract], [patch]], names=['tract', 'patch'])
        outdir = os.path.join(args.outdir, 'run-{}-{}'.format(tract, patch))
        outdir = outdir+'_'+args.label if args.label else outdir
        hp.utils.mkdir_if_needed(outdir)
    else:
        patches = hp.get_group_patches(group_id=args.group_id) 
        outdir = args.group_dir
        hp.utils.mkdir_if_needed(outdir)
        print('searching in', len(patches), 'patches')

    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)
    main(pool, patches, outdir, config_fn=args.config_fn, seed=args.seed)
