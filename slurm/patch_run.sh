#! /bin/bash 
# 
#SBATCH -J hugs-pipe-run-patches     # job name
#SBATCH -o /scratch/network/jgreco/run-%j.out
#SBATCH -e /scratch/network/jgreco/run-%j.err             
#SBATCH -N 8
#SBATCH --ntasks-per-node=16
#SBATCH -t 23:00:00 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end 
#SBATCH --mail-user=jgreco@princeton.edu 

cd /home/jgreco/projects/hugs-pipe/scripts

RUN_LABEL=`date +%Y%m%d-%H%M%S`
OUTDIR=$HUGS_PIPE_IO/patches-run-$RUN_LABEL
PATCHES_FN=$LOCAL_IO/patch-files/hsc-wide-patches-full.csv

mkdir $OUTDIR
cp $PATCHES_FN $OUTDIR/

mpiexec -n 128 python runner.py --mpi \
    --patches_fn $PATCHES_FN \
    -o $OUTDIR \
    -c $LOCAL_IO/pipe-configs/01-25-2017.yml
