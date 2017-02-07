#!/bin/bash

#phyloP_subTree.slurm


#SBATCH -n 16
#SBATCH -t 1-00:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o phyloP_subTree.out
#SBATCH -e phyloP_subTree.err

halPath=$1
refGenome=$2
mod=$3
outFile=$4
subTreeNode=$5

halPhyloP --subtree $subTreeNode $halPath $refGenome $mod $outFile
