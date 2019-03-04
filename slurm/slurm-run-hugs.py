##############################################
# Submit slurm jobs for running hugs
##############################################
import os
import subprocess
from argparse import ArgumentParser
from hugs.utils import project_dir, mkdir_if_needed
from hugs import slurm
from hugs.log import logger
script_dir = '/tigress/jgreco/hsc-s18a/slurm-scripts'
mkdir_if_needed(script_dir)

parser = ArgumentParser()
parser.add_argument('-r', '--run-name', dest='run_name', required=True)
parser.add_argument('-n', '--ntasks-per-node', dest='ncores', 
                    type=int, required=True)
parser.add_argument('-p', '--patch-fn', dest='patch_fn', required=True)
parser.add_argument('-c', '--config-fn', dest='config_fn', required=True)
parser.add_argument('--nodes', type=int, default=1)
parser.add_argument('--time', default='05:00:00')
parser.add_argument('--test', action='store_true')
args = parser.parse_args()

assert args.ncores <= 40, 'there are 40 cores per node!'
total_cores = int(args.ncores * args.nodes)

logger.info('submitting slurm job with {} cores across {} nodes'.\
            format(total_cores, args.nodes))

task_options = [total_cores, args.run_name, args.patch_fn, args.config_fn]

slurm_configs = [
    {'job-name': args.run_name,
     'nodes': args.nodes,
     'ntasks-per-node': args.ncores,
     'time': args.time, 
     'filename': os.path.join(script_dir, args.run_name + '-hugs'), 
     'task-options': task_options}
]

template_dir = os.path.join(project_dir, 'slurm/templates') 
templates = [os.path.join(template_dir, 'run-hugs')]

slurm.drive(templates, slurm_configs, test=args.test, wait=False)
