from __future__ import division, print_function

import os, sys
from argparse import ArgumentParser
from .utils import io

def parse_args(default_outdir=io):

    parser = ArgumentParser('python '+sys.argv[0]) 

    parser.add_argument('-t', 
                        '--tract', 
                        type=int, 
                        help='HSC tract')
    parser.add_argument('-p', 
                        '--patch', 
                        type=str, 
                        help='HSC patch')
    parser.add_argument('-g', 
                        '--group_id', 
                        type=int, 
                        help='group id', 
                        default=None)
    parser.add_argument('-c', 
                        '--config_fn', 
                        type=str, 
                        help='config file name', 
                        default=None)
    parser.add_argument('-o', 
                        '--outdir', 
                        type=str, 
                        help='output directory', 
                        default=default_outdir)
    parser.add_argument('-s', 
                        '--seed', 
                        type=int, 
                        help='rng seed', 
                        default= None)
    parser.add_argument('--num_synths', 
                        type=int, 
                        help='number of synths to inject', 
                        default=10)
    args = parser.parse_args()

    if args.group_id:
        args.group_dir = os.path.join(args.outdir, 'group_'+str(args.group_id))
    elif (args.tract is None) or (args.patch is None):
        parser.print_help()
        exit(0)

    return args
