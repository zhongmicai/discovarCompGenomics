#!/bin/bash

#phyloP.slurm


#SBATCH -n 4
#SBATCH -t 0-10:00
#SBATCH --mem=5000
#SBATCH -p serial_requeue
#SBATCH -o phyloP_.out
#SBATCH -e phyloP.err

halPath=$1
refGenome=$2
mod=$3
seq=$4
overallOut=$5
outFile=phyloP_tmp_$seq.wig


halPhyloP --refSequence $seq $halPath $refGenome $mod $outFile

cat $outFile >> $overallOut
rm -f $outFile
