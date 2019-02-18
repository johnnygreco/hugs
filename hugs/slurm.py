import time
import os, re
import subprocess
import numpy as np
from collections import OrderedDict
from . import LSST_INSTALLED
scratch_dir = '/scratch/gpfs/jgreco/hugs-io/'

slurm_options = OrderedDict([
    ('job-name', 'slurm-driver'),
    ('output', scratch_dir + 'slurm-%j.out'),
    ('error', scratch_dir + 'slurm-%j.err'),
    ('nodes', '1'),
    ('ntasks-per-node', '16'),
    ('time', '01:00:00'),
    ('mail-type', 'end'),
    ('mail-user', 'jgreco@princeton.edu')
])


def _print_cmd(cmd):
    print('----- running -----')
    print(cmd)
    print('-------------------')


def _write_slurm_script(cmd, filename=None, test=False, **batch_options):
    if not test:
        assert filename, 'must give file name!'
        with open(filename, 'w') as fout:
            fout.write('#!/bin/bash\n')
            for opts in batch_options.items():
                fout.write('#SBATCH --{0}={1}\n'.format(*opts))
            
            fout.write('\n')
            fout.write('{0}\n'.format(cmd))
    else:
        print('#!/bin/bash')
        for opts in batch_options.items():
            print('#SBATCH --{0}={1}'.format(*opts))
        print()
        print('{0}\n'.format(cmd))


def _get_job_status(jobid, wait=30):
    """Returns status of slurm job <jobid>
    Currently parses output of `sacct`.  Perhaps would
    be a good idea to move this to pyslurm (though this would 
    add a dependency.)
    ref: https://github.com/timothydmorton/lsst-utils
    """
    cmd = 'scontrol show job {0}'.format(jobid)
    try:
        output = subprocess.check_output(cmd, shell=True).decode('utf-8')
        m = re.search('JobState=(\w+)', output)
    except subprocess.CalledProcessError:
        m = False

    status = None
    if m:
        status = m.group(1)
    else:
        repeat = 0
        while not m and repeat < wait:    
            cmd = 'sacct -b -j {0}'.format(jobid)
            output = subprocess.check_output(cmd, shell=True).decode('utf-8')
            m = re.search('{0}\s+([A-Z]+)'.format(jobid), output)
            time.sleep(1)
            repeat += 1 
        if m:
            status = m.group(1)   

    if status is None:
        raise ValueError('Job not found: {0}'.format(jobid))
    else:
        return status


def _get_jobid(output):
    jid = re.search('batch job (\d+)', output)
    if not jid:
        print('Cannot find job id: {0}'.format(output))
        raise RuntimeError('Cannot find job id.')
    return jobid 


def _wait_for_jobs(job_ids, test=False):
    while len(job_ids) > 0:
        completed = []
        for jid in job_ids:
            if test:
                if np.random.random() < 0.3:
                    status = 'COMPLETED'
                else:
                    status = 'RUNNING'
            else:
                status = _get_job_status(jid)
            if status in ('RUNNING', 'PENDING', 'COMPLETING'):
                continue
            elif status == 'COMPLETED':
                print('job id {0} completed. '.format(jid))
                completed.append(jid)
            elif status in ('FAILED', 'CANCELLED', 'TIMEOUT'):
                time.sleep(5)
                status = _get_job_status(jid)
                if status in ('FAILED', 'CANCELLED', 'TIMEOUT'):
                    raise RuntimeError('Unexpected status: Job {0} is {1}.'.\
                                       format(jid, status))
                else:
                    continue
            else:
                raise Exception('Unknown status for {0}: {1}.'.\
                                format(jid, status))
        for jid in completed:
            job_ids.remove(jid)
        time.sleep(5)
    return len(job_ids)==0


def drive(templates, configs, test=False, wait=True):
    """
    Submit slurm jobs using template files, each with an associated config 
    dictionary. The jobs will run in order, waiting for the previous job 
    to finish before starting. This is useful if you have jobs that depend
    on the output of a previous job. 

    Parameters
    ----------
    templates : list of str
        Slurm template file names.
    configs: list of dict
        Config dictionaries for each slurm job.
    test : bool
        If true, will perform dry run for testing.
    wait : bool
        Wait for job to finish before starting the next one. 
    """

    assert LSST_INSTALLED, 'load lsst before running!'

    if test:
        print('\n************ TESTING ************\n')

    for template, config in zip(templates, configs):

        with open(template, 'r') as slurm_file:
            slurm_txt = slurm_file.read()

        if 'task-options' in config.keys():
            slurm_txt = slurm_txt.format(*config.pop('task-options'))

        slurm_fn = config.pop('filename', scratch_dir + 'slurm-driver')

        for key, val in config.items():
            slurm_options[key] = val

        print('writing slurm script:', slurm_fn)
        _write_slurm_script(slurm_txt, slurm_fn, test, **slurm_options)

        cmd = 'sbatch {}'.format(slurm_fn)
        _print_cmd(cmd)

        if test:
            job_ids = np.random.randint(10000, size=1).tolist()
        else:
            output = subprocess.check_output(cmd, shell=True)
            job_ids = re.findall('batch job (\d+)', output.decode('utf-8'))    

        if wait:
            finished = _wait_for_jobs(job_ids, test)
            assert finished, 'job did not finish!'
