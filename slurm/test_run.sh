#! /bin/bash 
# 
#SBATCH -J hugs-pipe-run      # job name
#SBATCH -o /scratch/network/jgreco/run-%j.out
#SBATCH -e /scratch/network/jgreco/run-%j.err             
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -t 0:30:00 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end 
#SBATCH --mail-user=jgreco@princeton.edu 

cd /home/jgreco/projects/hugs-pipe/scripts

group_id=9552

python runner.py --ncores 10 -g $group_id -c $LOCAL_IO/config-11-17-2016.yml \
    -o $HUGS_PIPE_IO/sex-results
