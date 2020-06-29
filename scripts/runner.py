"""
Run hugs pipeline. 
"""
from __future__ import division, print_function

import os, shutil
from time import time
import mpi4py.MPI as MPI
import schwimmbad
from hugs.pipeline import next_gen_search, find_lsbgs
from hugs.utils import PatchMeta
import hugs


def ingest_data(args):
    """
    Write data to database with the master process.
    """
    timer = time()
    success, sources, meta_data, synth_ids = args
    run_name, tract, patch, patch_meta = meta_data
    db_ingest = hugs.database.HugsIngest(session, run_name)
    if success:
        db_ingest.add_all(tract, patch, patch_meta, sources)
        if synth_ids is not None:
            db_ingest.add_injected_synths(synth_ids)
    else:
        db_ingest.add_tract(tract)
        db_ingest.add_patch(patch, patch_meta)
    delta_time = time() - timer
    hugs.log.logger.info('time to ingest = {:.2f} seconds'.format(delta_time))


def worker(p):
    """
    Workers initialize pipe configuration and run pipeline.
    """
    rank = MPI.COMM_WORLD.Get_rank()
    if p['seed'] is None:
        tract, p1, p2 = p['tract'], int(p['patch'][0]), int(p['patch'][-1])
        seed = [int(time()), tract, p1, p2, rank]
    else:
        seed = p['seed']

    config = hugs.PipeConfig(run_name=p['run_name'], 
                             config_fn=p['config_fn'],
                             random_state=seed, 
                             log_fn=p['log_fn'],
                             rerun_path=p['rerun_path'])
    config.set_patch_id(p['tract'], p['patch'])
    config.logger.info('random seed set to {}'.format(seed))
    
    if p['use_old_pipeline']:
        results = find_lsbgs.run(config)
    else:
        results = next_gen_search.run(config, False)

    pm = results.hugs_exp.patch_meta

    if (results.synths is not None) and results.success:
        if len(results.synths) > 0:
            synth_ids = results.synths.to_pandas().loc[:, ['synth_id']]
            for plane in config.synth_check_masks:
                masked = hugs.synths.find_masked_synths(results.synths, 
                                                        results.exp_clean, 
                                                        planes=plane)
                synth_ids['mask_' + plane.lower()] = masked
        else:
            synth_ids = None
    else:
        synth_ids = None

    patch_meta = PatchMeta(
        x0 = pm.x0,
        y0 = pm.y0,
        small_frac = pm.small_frac,
        cleaned_frac = pm.cleaned_frac,
        bright_obj_frac = pm.bright_obj_frac,
        good_data_frac = pm.good_data_frac
    )

    meta_data = [
        config.run_name,
        config.tract,
        config.patch,
        patch_meta,
    ]

    if results.success:
        df = results.sources.to_pandas()
        df['flags'] = df['flags'].astype(int)
    else:
        df = None

    config.reset_mask_planes()
    config.logger.info('writing results to database')
    return results.success, df, meta_data, synth_ids


if __name__=='__main__':
    from argparse import ArgumentParser
    from astropy.table import Table
    rank = MPI.COMM_WORLD.Get_rank()
    
    # parse command-line arguments
    parser = ArgumentParser('Run hugs pipeline')
    parser.add_argument('-t', '--tract', type=int, help='HSC tract')
    parser.add_argument('-p', '--patch', type=str, help='HSC patch')
    parser.add_argument('-c', '--config_fn', help='hugs config file',
                        default=hugs.utils.default_config_fn)
    parser.add_argument('--patches_fn', help='patches file')
    parser.add_argument('--use-old-pipeline', action="store_true")
    parser.add_argument('-r', '--run_name', type=str, default='hugs-pipe-run')
    parser.add_argument('--seed', help='rng seed', default=None)
    parser.add_argument('--rerun_path', help='full rerun path', default=None)
    parser.add_argument('--overwrite', type=bool, 
                        help='overwrite database', default=True)
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--ncores', default=1, type=int, 
                          help='Number of processes (uses multiprocessing).')
    group.add_argument('--mpi', default=False, action="store_true", 
                          help="Run with MPI.")
    args = parser.parse_args()
    config_params = hugs.utils.read_config(args.config_fn)
    outdir = config_params['hugs_io']

    #######################################################################
    # run on a single patch
    #######################################################################

    if args.tract is not None:
        assert args.patch is not None
        tract, patch = args.tract, args.patch
        patches = Table([[tract], [patch]], names=['tract', 'patch'])
        run_dir_name = '{}-{}-{}'.format(args.run_name, tract, patch)
        outdir = os.path.join(outdir, run_dir_name)
        hugs.utils.mkdir_if_needed(outdir)
        log_fn = os.path.join(outdir, 'hugs-pipe.log')
        patches['outdir'] = outdir
        patches['log_fn'] = log_fn


    #######################################################################
    # OR run on all patches in file
    #######################################################################

    elif args.patches_fn is not None:
        patches = Table.read(args.patches_fn)
        if rank==0:
            time_label = hugs.utils.get_time_label()
            outdir = os.path.join(
                outdir, '{}-{}'.format(args.run_name, time_label))
            hugs.utils.mkdir_if_needed(outdir)
            log_dir = os.path.join(outdir, 'log')
            hugs.utils.mkdir_if_needed(log_dir)
            log_fn = []
            for tract, patch in patches['tract', 'patch']:
                fn = os.path.join(log_dir, '{}-{}.log'.format(tract, patch))
                log_fn.append(fn)
            patches['outdir'] = outdir
            patches['log_fn'] = log_fn

    else:
        print('\n**** must give tract and patch --or-- a patch file ****\n')
        parser.print_help()
        exit()

    patches['rerun_path'] = args.rerun_path
    patches['seed'] = args.seed
    patches['config_fn'] = args.config_fn
    patches['run_name'] = args.run_name
    patches['use_old_pipeline'] = args.use_old_pipeline

    if rank==0:
        # open database session with master process
        db_fn = os.path.join(outdir, args.run_name+'.db')
        engine = hugs.database.connect(db_fn, args.overwrite)
        session = hugs.database.Session()
        shutil.copyfile(args.config_fn, os.path.join(outdir, 'config.yml'))

    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.ncores)
    list(pool.map(worker, patches, callback=ingest_data))
    pool.close()
