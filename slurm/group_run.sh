#! /bin/bash 
# 
#SBATCH -N 2  # node count 
#SBATCH --ntasks-per-node=15
#SBATCH -t 4:00:00 
#SBATCH --mail-type=begin
#SBATCH --mail-type=end 
#SBATCH --mail-user=jgreco@princeton.edu 

cd /home/jgreco/projects/hugs-pipe/scripts

outdir=/home/jgreco/hugs-pipe-out/group_52983
region_file=/home/jgreco/data/mycats/patches_group_52983.txt

if [ ! -d $outdir ]; then
    mkdir $outdir
fi

IFS=$'\n'
for region in $(cat $region_file); do 
    IFS=$' '
    arr=($region)
    python runner.py ${arr[0]} ${arr[1]} -o $outdir &
done
wait
