#! /bin/bash 
# 
#SBATCH -J hugs-pipe-run      # job name
#SBATCH -o /scratch/network/jgreco/run-%j.out
#SBATCH -e /scratch/network/jgreco/run-%j.err             
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -t 0:20:00 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end 
#SBATCH --mail-user=jgreco@princeton.edu 

cd /home/jgreco/projects/hugs-pipe/scripts

group_id=1246

python runner.py --ncores 16 -g $group_id -c $LOCAL_IO/config-01-16-2017.yml \
    -o $HUGS_PIPE_IO/sex-results 
