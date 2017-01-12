#!/bin/bash

#getAlignmentDepth.slurm

#SBATCH -n 16
#SBATCH -t 2-00:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o getAlignmentDepth.out
#SBATCH -e getAlignmentDepth.err

module load hal

halPath=$1
targetGenome=$2
referenceGenome=$3
outputWig=$4
outputBed=$5

halAlignmentDepth --outWiggle $outputWig --targetGenomes $targetGenome $halPath $referenceGenome
wigToBed.py $outputWig $outputBed 1 0 1
