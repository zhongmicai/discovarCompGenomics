#!/bin/bash

#subsetMaf.slurm
#Nate Edelman, 2/6/17
#Takes a maf file and certain species
#returns a maf file of the alignment of only those species

#SBATCH -n 16
#SBATCH -t 2-00:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o subsetMaf.out
#SBATCH -e subsetMaf.err

mafFile=$1
seqs=$2
outFile=$3

msa_view --in-format MAF --out-format PHYLIP $mafFile --seqs $seqs > $outFile
