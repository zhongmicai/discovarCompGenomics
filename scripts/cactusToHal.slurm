#!/bin/bash

#cactusToHal.slurm

#SBATCH -n 16
#SBATCH -t 2-00:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o intermediateCactusToHal.out
#SBATCH -e intermediateCactusToHal.err

source /n/sw/progressiveCactus-latest/progressiveCactus/environment

project=$1
halFile=$2

cactus2hal.py --append $project $halFile
