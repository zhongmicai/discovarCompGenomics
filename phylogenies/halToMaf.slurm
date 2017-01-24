#!/bin/bash

#halToMaf.slurm

#SBATCH -n 16
#SBATCH -t 1-00:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o halToMaf.out
#SBATCH -e halToMaf.err

refGenome=$1
refTargets=$2
halFile=$3
mafFile=$4

hal2maf --noAncestors --refTargets $refTargets --refGenome $refGenome $halFile $mafFile
