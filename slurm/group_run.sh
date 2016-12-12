#! /bin/bash 
# 
#SBATCH -J hugs-pipe-run      # job name
#SBATCH -o /scratch/network/jgreco/run-%j.out
#SBATCH -e /scratch/network/jgreco/run-%j.err             
#SBATCH -N 2
#SBATCH --ntasks-per-node=16
#SBATCH -t 3:00:00 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end 
#SBATCH --mail-user=jgreco@princeton.edu 

cd /home/jgreco/projects/hugs-pipe/scripts

for group_id in $(cat /home/jgreco/projects/hugs-pipe/slurm/thegroups.txt); do
    python runner.py --ncores 32 -g $group_id -c $LOCAL_IO/config-11-17-2016.yml \
        -o $HUGS_PIPE_IO/run-results
done
wait
