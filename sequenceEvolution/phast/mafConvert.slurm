#!/bin/bash

#mafConvert.slurm
#Nate Edelman, 2/6/17
#Takes a maf file and output format
#returns a file of the alignment in the specified format

#SBATCH -n 16
#SBATCH -t 2-00:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o mafConvert.out
#SBATCH -e mafConvert.err

mafFile=$1
format=$2
outFile=$3

msa_view --in-format MAF $mafFile --out-format $format > $outFile
