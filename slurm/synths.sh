#! /bin/bash 
# 
#SBATCH -J synths      # job name
#SBATCH -n 20          # total number of mpi tasks requested
#SBATCH -t 4:00:00 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end 
#SBATCH --mail-user=jgreco@princeton.edu 

cd /home/jgreco/projects/hugs-pipe/scripts
module load openmpi/gcc/1.10.2/64

outdir=/home/jgreco/hugs-pipe-out/group_52983
group_id=4736
num_synths=15

python synths.py --mpi \
-c /Volumes/tiger/hugs-pipke-io/synth-results/config-11-17-2016.yml \
-g $group_id \
--num_synths $num_synths
