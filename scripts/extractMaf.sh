#extractMaf.sh
#Nate Edelman, 6/29/16
#Takes a hal file, a reference genome, and a scaffold name.
#returns a maf file of the alignment

module load hal

halFile=$1
refGenome=$2
scafName=$3
mafFile=$(basename $halFile .hal)_$scafName.maf

hal2maf --noAncestors --refGenome $refGenome --refSequence $scafName $halFile $mafFile
