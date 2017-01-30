#!/bin/bash

#phyloptrain.slurm


#SBATCH -n 16
#SBATCH -t 1-00:00
#SBATCH --mem=50000
#SBATCH -p serial_requeue
#SBATCH -o phylopTrain.out
#SBATCH -e phylopTrain.err

halFile=$1
refGenome=$2
bedFile=$3
outFile=$4

halPhyloPTrain.py --numProc 16 --noAncestors $halFile $refGenome $bedFile $outFile
