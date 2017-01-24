#! /bin/bash 
# 
#SBATCH -J hugs-pipe-run      # job name
#SBATCH -o /scratch/network/jgreco/run-%j.out
#SBATCH -e /scratch/network/jgreco/run-%j.err             
#SBATCH -N 4
#SBATCH --ntasks-per-node=16
#SBATCH -t 5:30:00 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end 
#SBATCH --mail-user=jgreco@princeton.edu 

cd /home/jgreco/projects/hugs-pipe/scripts

RUN_LABEL=`date +%Y%m%d-%H%M%S`
OUTDIR=$HUGS_PIPE_IO/batch-run-$RUN_LABEL
PATCHES_FN=$LOCAL_IO/patches_z0.05_Mh12.5-15.0.csv

mkdir $OUTDIR
cp $PATCHES_FN $OUTDIR/

mpiexec -n 64 python runner.py --mpi \
    --patches_fn $PATCHES_FN \
    -o $OUTDIR \
    -c $LOCAL_IO/hugs-pipe-config-01-19-2017.yml
