#!/bin/bash

#cactusToHal.slurm

#SBATCH -n 16
#SBATCH -t 2-00:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o intermediateCactusToHal.out
#SBATCH -e intermediateCactusToHal.err

module load hal

halPath=$1
tagetGenome=$2
referenceGenome=$3
outputWig=$4
outputBed=$5

halAlignmentDepth --targetGenomes $targetGenome --outWiggle $outputFile $halPath $referenceGenome
wigToBed.py $outputWig $outputBed 1 0 1
