#this is the master script for getting the alignment stats from a wig file.
#It includes basic stats like number of bases with each alignment depth
#and slightly more complex things like sizes of blocks of certain depth.
#I also want to distinguish between genic and non-genic regions.

export PATH=$PATH:~/Documents/Mallet_Lab/DiscovarCompGenomics/scripts

wigFile=$1
baseFile=$(basename $wigFile .wig)

#first, see what the highest depth is
highDepth=$(grep -v fixed $wigFile | awk '{for(i=1;i<=NF;i++) if($i>maxval) maxval=$i;}; END { print maxval;}' -)

#next, get the most basic statistics on number of bases with certain depths
wigAlignmentStats.py $wigFile

#Then, make bed files for a number of different alignment depth blocks
for i in $(seq 0 $highDepth);
do wigToBed.py $wigFile $baseFile\_atLeast$i.bed  $i 0 1;
done

#Make a bed file with all the gaps using bedtools complement
#requires a "genome file", which is just <scafname)\t<scafsize>.
#conveniently, that's just columns 1 and 3 of the atLeast0 file.
awk -v OFS="\t" '{print $1,$3}' $baseFile\_atLeast0.bed > $baseFile.genome

bedtools complement -i  $baseFile\_atLeast1.bed -g $baseFile.genome > $baseFile\_gaps.bed
