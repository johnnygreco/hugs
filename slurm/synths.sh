#! /bin/bash 
# 
#SBATCH -J synths      # job name
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH -t 2:00:00 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end 
#SBATCH --mail-user=jgreco@princeton.edu 

cd /home/jgreco/projects/hugs-pipe/scripts

if [[ $# -ne 1 ]]; then
    echo "must give group id"
    return 1
fi

group_id=$1
num_synths=15
nrepeat=10

for label in $(seq 1 $nrepeat); do
    python synths.py --ncores 16 -g $group_id --num_synths $num_synths --label $label \
        -c $HUGS_PIPE_IO/synth-results/config-11-17-2016.yml \
        -o /scratch/network/jgreco/hugs-pipe-output/synth-results
done
wait
