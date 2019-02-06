"""
Run hugs pipeline. 
"""
from __future__ import division, print_function

import os
from time import time
import mpi4py.MPI as MPI
import schwimmbad
from hugs.pipeline import next_gen_search
from hugs.utils import PatchMeta, project_dir
from astropy.table import vstack
import hugs


def ingest_data(args):
    """
    Write data to database with the master process.
    """
    timer = time()
    success, cats, meta_data = args
    sources, recovered, injected, synth_cat = cats
    run_name, tract, patch, patch_meta = meta_data
    db_ingest = hugs.database.HugsIngest(session, run_name)
    if success and (len(sources) > 0):
        db_ingest.add_all(tract, patch, patch_meta, sources.to_pandas())
        all_recovered.append(recovered)
        all_injected.append(injected)
        all_synth_cat.append(synth_cat)
    else:
        failed_patches['tract'].append(tract)
        failed_patches['patch'].append(patch)
        failed_patches['good_data_frac'].append(patch_meta.good_data_frac)
        failed_patches['success'].append(success)
    delta_time = time() - timer
    print('time to ingest =', delta_time)


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
                             rerun_path=p['rerun_path'])
    config.set_patch_id(p['tract'], p['patch'])
    config.logger.info('random seed set to {}'.format(seed))
    
    results = next_gen_search.run(config)

    pm = results.hugs_exp.patch_meta
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
        sources = results.sources
        sources['flags'] = sources['flags'].astype(int)

        synth_cat = config.synth_cat
        synth_cat['tract'] = config.tract
        synth_cat['patch'] = config.patch

        synth_cat.rename_column('x', 'x_image')
        synth_cat.rename_column('y', 'y_image')
        (match, match_synth), _  = hugs.cattools.xmatch(
            sources, synth_cat, max_sep=config.synth_max_match_sep)

        recovered = results.sources[match]
        injected = synth_cat[match_synth]
        injected['tract'] = config.tract
        injected['patch'] = config.patch

        txt = '{} injected, {} recovered'.format(len(synth_cat), 
                                                 len(injected))
        config.logger.info(txt)


    else:
        sources = None
        recovered = None
        injected = None
        synth_cat = None

    config.logger.info('passing results to master process')
    cats = [sources, recovered, injected, synth_cat]
    return results.success, cats, meta_data


if __name__=='__main__':
    from argparse import ArgumentParser
    from astropy.table import Table
    rank = MPI.COMM_WORLD.Get_rank()
    config_dir = os.path.join(project_dir, 'pipe-configs')
    
    # parse command-line arguments
    parser = ArgumentParser('Run hugs pipeline')
    parser.add_argument('-t', '--tract', type=int, help='HSC tract')
    parser.add_argument('-p', '--patch', type=str, help='HSC patch')
    parser.add_argument('-c', '--config_fn', help='hugs config file',
                        default=os.path.join(config_dir, 'hugs-run-dev.yml'))
    parser.add_argument('--patches_fn', help='patches file')
    parser.add_argument('-r', '--run_name', type=str, default='synth-run')
    parser.add_argument('--seed', help='rng seed', default=None)
    parser.add_argument('--rerun_path', help='full rerun path', default=None)
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

    if rank==0:
        # master process lists for results
        db_fn = os.path.join(outdir, args.run_name+'.db')
        engine = hugs.database.connect(db_fn, True)
        session = hugs.database.Session()
        all_recovered = []
        all_injected = []
        all_synth_cat = []
        failed_patches = {'tract': [], 
                          'patch': [], 
                          'good_data_frac': [], 
                          'success': []}

    pool = schwimmbad.choose_pool(mpi=args.mpi, processes=args.ncores)
    list(pool.map(worker, patches, callback=ingest_data))
    pool.close()

    if rank==0:
        fn = lambda lab: os.path.join(outdir, args.run_name + lab + '.csv')
        if len(all_recovered) > 0:
            all_recovered = vstack(all_recovered)
            all_injected = vstack(all_injected)
            all_synth_cat = vstack(all_synth_cat)

            all_recovered.write(fn('-recovered'), overwrite=True)
            all_injected.write(fn('-injected'), overwrite=True)
            all_synth_cat.write(fn('-synth-cat'), overwrite=True)

        failed_patches = Table(failed_patches)
        failed_patches.write(fn('-failed-patches'), overwrite=True)
