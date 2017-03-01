#!/bin/bash

#halToMaf.slurm

#SBATCH -n 16
#SBATCH -t 0-03:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o halToMaf.out
#SBATCH -e halToMaf.err

refGenome=$1
halFile=$2
mafFile=$3
scaffold=$4

hal2maf --noAncestors --refGenome $refGenome --refSequence $scaffold $halFile $mafFile
