#!/bin/bash

#findDeletions.slurm

#SBATCH -n 4
#SBATCH -t 2-00:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o findDeletions.out
#SBATCH -e findDeletions.err

halFile=$1
refGenome=$2
ingroupGenomes=$3
outputFile=$4

findRegionsExclusivelyInGroup $halFile $refGenome $ingroupGenomes > $outputFile
