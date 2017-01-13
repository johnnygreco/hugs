from __future__ import division, print_function

import os, sys
from argparse import ArgumentParser
from .utils import io

def parse_args(which='all', default_outdir=io, parser=None):

    parser = parser if parser else ArgumentParser('python '+sys.argv[0]) 

    if 'patch' in which or which=='all':
        parser.add_argument('-t', 
                            '--tract', 
                            type=int, 
                            help='HSC tract')
        parser.add_argument('-p', 
                            '--patch', 
                            type=str, 
                            help='HSC patch')

    if 'config_fn' in which or which=='all':
        parser.add_argument('-c', 
                            '--config_fn', 
                            type=str, 
                            help='config file name', 
                            default=None)

    if 'outdir' in which or which=='all':
        parser.add_argument('-o', 
                            '--outdir', 
                            type=str, 
                            help='output directory', 
                            default=default_outdir)

    if 'seed' in which or which=='all':
        parser.add_argument('-s', 
                            '--seed', 
                            type=int, 
                            help='rng seed', 
                            default= None)

    if 'synths' in which or which=='all':
        parser.add_argument('--num_synths', 
                            type=int, 
                            help='number of synths to inject', 
                            default=10)

    if 'label' in which or which=='all':
        parser.add_argument('--label',
                            type=str,
                            help='extra label for run',
                            default=None)

    if 'parallel' in which or which=='all':
        group = parser.add_mutually_exclusive_group()
        group.add_argument("--ncores", 
                           dest="n_cores", 
                           default=1,
                           type=int, 
                           help="Number of processes (uses multiprocessing).")
        group.add_argument("--mpi", 
                           dest="mpi",
                           default=False,
                           action="store_true", 
                           help="Run with MPI.")

    if 'group id' in which or which=='all':
        parser.add_argument('-g', 
                            '--group_id', 
                            type=int, 
                            help='group id', 
                            default=None)

    if 'sex' in which or which=='all':
        parser.add_argument('--use_sex',
                            help='use sextractor for final detection',
                            action='store_true')

    args = parser.parse_args()
    if hasattr(args, 'group_id'):
        args.group_dir = os.path.join(args.outdir, 'group-'+str(args.group_id))

    if len(sys.argv)==1:
        parser.print_help()
        exit(0)

    return args
