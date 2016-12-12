#!/bin/bash
for group_id in $(cat /home/jgreco/projects/hugs-pipe/slurm/thegroups.txt); do
    python runner.py --ncores 1 -g $group_id -c $LOCAL_IO/config-11-17-2016.yml \
        -o $HUGS_PIPE_IO/run-results
done
