#!/bin/bash

######################################################################################################################
## This is a wrapper script to run the data processing from the multiple alignment output of cactus to fasta files.
##
## Inputs: 1) hal file 2) gff file 3) reference genome name 4) min coverage 5) max gap length 6) min block length
## Outputs: 1) depth (number of genomes aligned) for each refernce base 2) bed file of blocks with good coverage
##      3) a sorted version of the good coverage blocks 4) the gff file converted to bed 5) the high coverage blocks within genes
##      6) the high coverage blocks outside of genes 7) maf alignments of (3) 8) maf alignments of (5) 9) maf alignments of (6)
##      10) fasta alignments of (3) 11) fasta alignments of (5) 12) fasta alignments of (6)
##
## Usage: coveragePipeline.sh <hal> <gff> <reference genome> <min coverage> <max gap length> <min block length> outputFolder
######################################################################################################################

########### load dependencies #############
module load halTools
module load bedTools
module load bedOps

########### read in Variables ############
halFile=$1
gffFile=$2
refGenome=$3
minCov=$4
maxGap=$5
minLength=$6
base=$(basename $halFile .hal)
gffBase=$(basename $gffFile .gff)

############ prepare data structure #########
mkdir -p outputFolder

############ outputs ###################
depthFile=outputFolder/$base\_alignmentDepth.wig
extractedBlocks=outputFolder/$base\_goodCoverage.bed
allBlocksSort=outputFolder/$base\_goodCoverage.sorted.bed
gffBed=outputFolder/$gffBase.bed
geneBlocks=outputFolder/genicBlocks.bed
interGeneBlocks=outputFolder/interGenicBlocks.bed
allMafs=outputFolder/$base\_allGoodCoverage.maf
geneMafs=outputFolder/$base\_genicBlocks.maf
interGeneMafs=outputFolder/$base\_interGenicBlocks.maf
allFasta=outputFolder/$base\_allGoodCoverage.fasta
geneFasta=outputFolder/$base\_genicBlocks.fasta
interGeneFasta=outputFolder/$base\_interGenicBlocks.fasta

############ script ################
#find depth at each reference base pair
halAlignmentDepth $halFile $refGenome > $depthFile

#find high-depth genomic segments
wigToBed.py $depthFile $extractedBlocks $minCov $maxGap $minLength

## extract interesting slices of the data in maf format
#extract all high-depth segments
sort -k1,1 -k2,2n $extractedBlocks > $allBlocksSort
hal2maf --refGenome=$refGenome --refTargets=$allBlocksSort $halFile > $allMafs
mafToFasta $allMafs $allFasta

#extract only protein coding genes with high depth
gff2bed < $gffFile > $gffBed #this converts the gff to a sorted maf
bedTools intersect -sorted -a $allBlocksSort -b $gffBed > $geneBlocks
hal2maf --refGenome=$refGenome --refTargets=$geneBlocks $halFile > $geneMafs
mafToFasta $geneMafs $geneFasta

#extract regions outside of known protein-coding genes
bedTools intersect -v -sorted -a $allBlocksSort -b $gffBed > $interGeneBlocks
hal2maf --refGenome=$refGenome --refTargets=$interGeneBlocks $halFile > $interGeneMafs
mafToFasta $interGeneMafs $interGeneFasta











