#! /bin/bash 
# 
#SBATCH -J synths      # job name
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH -t 1:00:00 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end 
#SBATCH --mail-user=jgreco@princeton.edu 

cd /home/jgreco/projects/hugs-pipe/scripts

group_id=4736
num_synths=15

python synths.py --ncores 10 -g $group_id --num_synths $num_synths\
    -c $HUGS_PIPE_IO/synth-results/config-11-17-2016.yml
