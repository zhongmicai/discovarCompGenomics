#!/bin/bash

#slurm script to run coveragePipeline.sh

#SBATCH -p serial_requeue
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=10000
#SBATCH -t 0-10:00
#SBATCH -e coveragePipeline.err
#SBATCH -o coveragePipeline.out

coveragePipeline.sh $1 $2 $3 $4 $5
