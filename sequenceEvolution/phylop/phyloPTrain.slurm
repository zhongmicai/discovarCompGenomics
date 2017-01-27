#!/bin/bash

#phyloptrain.slurm


#SBATCH -n 4
#SBATCH -t 1-00:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o phylopTrain.out
#SBATCH -e phylopTrain.err

halFile=$1
refGenome=$2
bedFile=$3
outFile=$4

halPhyloPTrain.py --numProc 4 --noAncestors $halFile $refGenome $bedFile $outFile
