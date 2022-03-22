#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6G
#SBATCH --time=24:00:00
#SBATCH --job-name=simulations
#SBATCH --mail-user=matthieu.gendre@unifr.ch
#SBATCH --mail-type=begin,end
#SBATCH --output=/home/mgendre/Cluster/output_%j.o
#SBATCH --error=/home/mgendre/Cluster/error_%j.e
#SBATCH --array=1-150%100


for r in 0.05 0.333333 0.5;
do 
Rscript --vanilla /home/mgendre/Cluster/scripts/main.R $r script_$SLURM_ARRAY_TASK_ID; 
done;


#Rscript /home/mgendre/Cluster/scripts/overallMean.R
