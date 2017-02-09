#!/bin/bash

#phyloFit.slurm
#Nate Edelman, 2/6/17
#Takes a maf file and gff file
#returns a neutral mod file made from the 4D sites found from the gff

#SBATCH -n 16
#SBATCH -t 2-00:00
#SBATCH --mem=10000
#SBATCH -p serial_requeue
#SBATCH -o phyloFit.out
#SBATCH -e phyloFit.err

mafFile=$1
gffFile=$2
outBase=$3
tree=$4

#msa_view $mafFile --in-format MAF --4d --features $gffFile > $outBase\_4d-codons.ss
msa_view $outBase\_4d-codons.ss --in-format SS --out-format SS --tuple-size 1 > $outBase\_4d-sites.ss
# phyloFit --tree $tree --msa-format SS \
#           --out-root $outBase\_nonconserved-4d $outBase\_4d-sites.ss
