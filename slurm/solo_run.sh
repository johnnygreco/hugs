#! /bin/bash 
# 
#SBATCH -J hugs-pipe-run      # job name
#SBATCH -o /scratch/network/jgreco/run-%j.out
#SBATCH -e /scratch/network/jgreco/run-%j.err             
#SBATCH -N 3
#SBATCH --ntasks-per-node=16
#SBATCH -t 0:50:00 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end 
#SBATCH --mail-user=jgreco@princeton.edu 

cd /home/jgreco/projects/hugs-pipe/scripts

group_id=4737
outdir=$HUGS_PIPE_IO/run-results/group-$group_id

mkdir -p $outdir

mpiexec -n 48 python runner.py --mpi -g $group_id \
    -c $LOCAL_IO/hugs-pipe-config-01-19-2017.yml \
    -o $outdir
