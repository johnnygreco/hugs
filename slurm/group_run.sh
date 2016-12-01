#! /bin/bash 
# 
#SBATCH -J hugs-pipe-run      # job name
#SBATCH -o /scratch/network/jgreco/hugs-pipe-output/run-%j.out
#SBATCH -e /scratch/network/jgreco/hugs-pipe-output/run-%j.err             
#SBATCH -N 2
#SBATCH --ntasks-per-node=16
#SBATCH -t 2:00:00 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end 
#SBATCH --mail-user=jgreco@princeton.edu 

cd /home/jgreco/projects/hugs-pipe/scripts

groups=(3102 2713 4736 6165 7572 8528)

for group_id in ${groups[*]}; do
    python runner.py --ncores 32 -g $group_id -c $HUGS_PIPE_IO/config-11-17-2016.yml \
        -o /scratch/network/jgreco/hugs-pipe-output/run-results
done
wait
