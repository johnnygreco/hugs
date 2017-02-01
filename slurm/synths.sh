#! /bin/bash 
# 
#SBATCH -J hugs-pipe-run-synths      # job name
#SBATCH -o /scratch/network/jgreco/run-%j.out
#SBATCH -e /scratch/network/jgreco/run-%j.err             
#SBATCH -N 5
#SBATCH --ntasks-per-node=16
#SBATCH -t 6:00:00 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end 
#SBATCH --mail-user=jgreco@princeton.edu 

cd /home/jgreco/projects/hugs-pipe/scripts

RUN_LABEL=`date +%Y%m%d-%H%M%S`
OUTDIR=$HUGS_PIPE_IO/batch-synths-$RUN_LABEL
ZLABEL='z0.05'
MLABEL='Mh12.5-15.0'
INDIR=$LOCAL_IO/group-patches
PATCHES_FN=$INDIR/patches_$ZLABEL\_$MLABEL.csv
NUM_SYNTHS=15

mkdir $OUTDIR
cp $PATCHES_FN $OUTDIR/
cp $INDIR/cat_$ZLABEL\_$MLABEL\_group_info.txt $OUTDIR/
cp $INDIR/cat_$ZLABEL\_$MLABEL\_tracts_n_patches.npy $OUTDIR/

mpiexec -n 80 python runner.py --mpi \
    --patches_fn $PATCHES_FN \
    -o $OUTDIR \
    -c $LOCAL_IO/pipe-configs/01-25-2017.yml \
    --num_synths $NUM_SYNTHS
