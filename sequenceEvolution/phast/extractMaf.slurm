#!/bin/bash

#extractMaf.slurm
#Nate Edelman, 2/6/17
#Takes a hal file, a reference genome, and a scaffold name.
#returns a maf file of the alignment

#SBATCH -n 16
#SBATCH -t 2-00:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o extractMaf.out
#SBATCH -e extractMaf.err


halFile=$1
refGenome=$2
scafName=$3
mafFile=$(basename $halFile .hal)_$scafName.maf

hal2maf --noAncestors --refGenome $refGenome --refSequence $scafName $halFile $mafFile
