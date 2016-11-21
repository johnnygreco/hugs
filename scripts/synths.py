"""
Run hugs-pipe on data with synthetic UGDs.
"""
from __future__ import division, print_function

import os
import numpy as np
import pandas as pd
import schwimmbad
from multiprocessing import Pool
from astropy.table import Table
import hugs_pipe as hp
from hugs_pipe.utils import calc_mask_bit_fracs

pset_lims = {'n': [0.6, 1.2],
             'mu0_i': [23, 26],
             'ell': [0.2, 0.2],
             'r_e': [2, 5]}

def worker(p):
    data_id = {'tract': p['tract'], 'patch': p['patch'], 'filter': 'HSC-I'}
    prefix = os.path.join(p['outdir'], 'hugs-{}-{}'.format(p['tract'],
                                                           p['patch']))
    log_fn = prefix+'.log'
    config = hp.Config(config_fn=p['config_fn'], log_fn=log_fn)
    config.set_data_id(data_id)
    
    synths_kwargs = {'num_synths': p['num_synths'],
                     'seed': p['seed'],
                     'pset_lims': pset_lims}

    sf = hp.SynthFactory(**synths_kwargs)

    results = hp.run(config, debug_return=True, synth_factory=sf)

    # write source catalog
    sources = results.sources.to_pandas()
    fn = prefix+'-cat.csv'
    sources.to_csv(fn, index=False)

    # write synth catalog
    fn = prefix+'-synths.csv'
    sf_df = results.sf.get_psets()
    sf_df['tract'] = tract
    sf_df['patch'] = patch
    sf_df['patch_id'] = np.arange(len(sf_df))
    sf_df.to_csv(fn, index=False)

    # write mask fractions 
    fn = prefix+'-mask-fracs.csv'
    mask_fracs = calc_mask_bit_fracs(results.exp_clean)
    mask_fracs = pd.DataFrame(mask_fracs)
    mask_fracs['tract'] = tract
    mask_fracs['patch'] = patch 
    mask_fracs.to_csv(fn, index=False)

def main(pool, patches, outdir, config_fn, num_synths=10, seed=None):
    patches['outdir'] = outdir
    patches['num_synths'] = num_synths
    patches['seed'] = seed
    patches['config_fn'] = config_fn

    pool.map(worker, patches)
    pool.close()

    # write synth param lims
    global pset_lims 
    pset_lims = pd.DataFrame(pset_lims, index=['min', 'max'])
    fn = os.path.join(outdir, 'synth_param_lims.csv')
    pset_lims.to_csv(fn, index=True, index_label='limit')

if __name__=='__main__':
    args = hp.parse_args(os.path.join(hp.io, 'synth-results'))

    if args.group_id is None:
        assert (args.tract is not None) and (args.patch is not None)
        tract, patch = args.tract, args.patch
        patches = Table([[tract], [patch]], names=['tract', 'patch'])
        outdir = os.path.join(args.outdir, 'synths-{}-{}'.format(tract, patch))
        hp.utils.mkdir_if_needed(outdir)
        pool = schwimmbad.choose_pool(processes=args.nthreads)
    else:
        patches = hp.get_group_patches(group_id=args.group_id) 
        outdir = args.group_dir
        hp.utils.mkdir_if_needed(outdir)
        print('searching in', len(patches), 'patches')
        pool = schwimmbad.choose_pool(processes=args.nthreads)

    main(pool, patches, outdir, config_fn=args.config_fn, 
         num_synths=args.num_synths, seed=args.seed)
