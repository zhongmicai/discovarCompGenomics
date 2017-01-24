#!/bin/bash

#halToMaf_orth.slurm

#SBATCH -n 16
#SBATCH -t 1-00:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o halToMaf_orth_trees.out
#SBATCH -e halToMaf_orth_trees.err

refGenome=$1
refTargets=$2
halFile=$3
mafFile=$4

hal2maf --noAncestors --printTree --onlyOrthologs --refTargets $refTargets --refGenome $refGenome $halFile $mafFile
