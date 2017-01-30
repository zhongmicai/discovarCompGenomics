#!/bin/bash

#phyloP.slurm


#SBATCH -n 4
#SBATCH -t 1-00:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o phyloP.out
#SBATCH -e phyloP.err

halPath=$1
refGenome=$2
mod=$3
outFile=$4

halPhyloP $halPath $refGenome $mod $outFile
