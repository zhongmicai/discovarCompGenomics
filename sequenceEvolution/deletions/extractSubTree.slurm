#!/bin/bash

#extractSubTree.slurm

#SBATCH -n 4
#SBATCH -t 1-00:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o extractSubTree.out
#SBATCH -e extractSubTree.err

halFile=$1
root=$2
outputFile=$3

halExtract --root $root $halFile $outputFile
