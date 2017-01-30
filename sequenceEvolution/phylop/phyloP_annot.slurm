#!/bin/bash

#phyloP_annot.slurm


#SBATCH -n 16
#SBATCH -t 1-00:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o phyloP_annot.out
#SBATCH -e phyloP_annot.err

halPath=$1
refGenome=$2
mod=$3
outFile=$4
bedFile=$5

halPhyloP --refBed $bedFile $halPath $refGenome $mod $outFile
