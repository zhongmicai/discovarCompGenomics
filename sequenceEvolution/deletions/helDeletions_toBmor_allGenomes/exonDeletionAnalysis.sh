#!/bin/bash

#take a bed file with deletion coordinates,
#and create protein alignments for all genes that have an exon > p percent deleted

#make sure that bedtools is available

deletionBed=$1
exonBed=$2
p=$3


exonsOfInterest=$(basename deletionBed .bed)\_deletedExons_$p\percent.bed
proteinNames=$(basename deletionBed .bed)\_deletedExons_$p\percent_proteinNames.txt

#intersectBed -f $p -a $exonBed -b $deletionBed > $exonsOfInterest

while read line;
do echo $line |awk -F "product=" '{print $2}' |cut -d ";" -f 1|cut -d "%" -f 1 >> $proteinNames
done < $exonsOfInterest
