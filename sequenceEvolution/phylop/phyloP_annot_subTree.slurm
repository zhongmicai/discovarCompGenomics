#!/bin/bash

#phyloP_annot_subTree.slurm


#SBATCH -n 16
#SBATCH -t 1-00:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o phyloP_annot_subTree.out
#SBATCH -e phyloP_annot_subTree.err

halPath=$1
refGenome=$2
mod=$3
outFile=$4
bedFile=$5
subTreeNode=$6

halPhyloP --refBed $bedFile --subtree $subTreeNode $halPath $refGenome $mod $outFile
