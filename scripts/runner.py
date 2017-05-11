"""
Run hugs-pipe on data with synthetic UGDs.
"""
from __future__ import division, print_function

import os
from time import time
import multiprocessing
import numpy as np
import pandas as pd
import mpi4py.MPI as MPI
import schwimmbad
import hugs_pipe as hp
from hugs_pipe.utils import calc_mask_bit_fracs


def callback(randoms_results):
    if randoms_results is None:
        pass
    else:
        df = randoms_results.df
        randoms_results.db.set_detected(df.loc[df['detected']==1, 'id'])


def worker(p):

    rank = MPI.COMM_WORLD.Get_rank()
    if p['seed'] is None:
        tract, p1, p2 = p['tract'], int(p['patch'][0]), int(p['patch'][-1])
        seed = [int(time()), tract, p1, p2, rank]
    else:
        seed = p['seed']

    prefix = 'hugs-pipe' if p['num_synths']==0 else 'synths'
    prefix = os.path.join(p['outdir'], prefix)
    log_fn = prefix+'.log'
    config = hp.Config(config_fn=p['config_fn'],
                       log_fn=log_fn,
                       random_state=seed)
    config.set_patch_id(p['tract'], p['patch'])
    config.logger.info('random seed set to {}'.format(seed))
    config.sex_measure['relpath'] = p['relpath']
    config.run_imfit['relpath'] = p['relpath']
    
    if p['num_synths']>0:
        synths_kwargs = {'num_synths': p['num_synths'],
                         'random_state': config.rng,
                         'pset_lims': config.synth_pset_lims}
        sf = hp.SynthFactory(**synths_kwargs)
    else:
        sf = None

    results = hp.run(config, synth_factory=sf, randoms_only=p['randoms_only'])

    if not p['randoms_only']:
        if p['num_synths']>0:
            # write synth catalog
            fn = prefix+'.csv'
            sf_df = sf.get_psets()
            sf_df['synth_id'] = np.arange(len(sf_df))
            sf_df.to_csv(fn, index=False)

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

    return results.randoms_results


if __name__=='__main__':
    import shutil
    from astropy.table import Table
    args = hp.parse_args()

    rank = MPI.COMM_WORLD.Get_rank()
    temp_io = os.environ.get('TEMP_IO')
    sex_io = os.environ.get('SEX_IO_DIR')
    runlabel = 'run' if args.num_synths==0 else 'synths'

    if args.tract is not None:
        assert args.patch is not None
        tract, patch = args.tract, args.patch
        patches = Table([[tract], [patch]], names=['tract', 'patch'])
        outdir = os.path.join(args.outdir,
                              'solo-{}-{}-{}'.format(runlabel, tract, patch))
        outdir = outdir+'_'+args.label if args.label else outdir
        hp.utils.mkdir_if_needed(outdir)
        patches['outdir'] = outdir
        maindir = outdir
    else:
        assert args.patches_fn is not None
        patches = Table.read(args.patches_fn)
        if rank==0:
            if args.outdir==hp.io:
                time_label = hp.utils.get_time_label()
                rundir = 'batch-{}-{}'.format(runlabel, time_label)
                maindir = os.path.join(hp.io, rundir)
                hp.utils.mkdir_if_needed(maindir)
            else:
                maindir = args.outdir
            outdirs = []
            for tract, patch in patches['tract', 'patch']:
                path = os.path.join(maindir, str(tract))
                hp.utils.mkdir_if_needed(path)
                path = os.path.join(path, patch)
                hp.utils.mkdir_if_needed(path)
                outdirs.append(path)
            patches['outdir'] = outdirs

    if rank==0:
        tempdir_relpath = maindir.split('/')[-1].replace(',', '-')
        tempdir = os.path.join(temp_io, tempdir_relpath)
        hp.utils.mkdir_if_needed(tempdir)
        sexin = os.path.join(sex_io, 'sexin/'+tempdir_relpath)
        hp.utils.mkdir_if_needed(sexin)
        sexout = os.path.join(sex_io, 'sexout/'+tempdir_relpath)
        hp.utils.mkdir_if_needed(sexout)
        patches['relpath'] = tempdir_relpath

    patches['num_synths'] = args.num_synths
    patches['seed'] = args.seed
    patches['config_fn'] = args.config_fn
    patches['randoms_only'] = args.randoms_only

    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.n_cores)
    list(pool.map(worker, patches, callback=callback))
    pool.close()

    if rank==0:
        shutil.rmtree(tempdir)
        shutil.rmtree(sexin)
        shutil.rmtree(sexout)
