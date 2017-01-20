#! /bin/bash 
# 
#SBATCH -J hugs-pipe-run      # job name
#SBATCH -o /scratch/network/jgreco/run-%j.out
#SBATCH -e /scratch/network/jgreco/run-%j.err             
#SBATCH -N 4
#SBATCH --ntasks-per-node=16
#SBATCH -t 8:00:00 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end 
#SBATCH --mail-user=jgreco@princeton.edu 

cd /home/jgreco/projects/hugs-pipe/scripts

RUN_LABEL=`date +%Y%m%d-%H%M%S`
PATCHES_FN=$HUGS_PIPE_IO/batch-run-$RUN_LABEL/patches.csv
GROUP_FN=$LOCAL_IO/thegroups.txt

python prep_batch_run.py $GROUP_FN --run_label $RUN_LABEL

mpiexec -n 64 python runner.py --mpi \
    --patches_fn $PATCHES_FN \
    -c $LOCAL_IO/hugs-pipe-config-01-19-2017.yml
