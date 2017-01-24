#Alignment Stats  

Here, I just wand to get some basic statistics about what was aligned. I'll compare the alignment of melpomene discovar erato reference to melpomene reference eueides, and bombyx to mel ref. I'll also look at the regions where all 25 of the genomes align.

##Strategy:

for the pairwise comparisons, do halAlignmentDepth and wigToBed to get regions of Alignment

Then, extract exons, introns, and intergenic regions from Hmel2.gff.

To get the right bed files, do:

exons:
```shell
grep exon gffFile > exon.gff
```

introns:
```shell
awk '$3=="gene" {print}' gffFile.gff > gene.gff
subtractBed -a gene.gff -b exon.gff > intron.gff
```

intergenic:
```python
#make full scaffold bed file
from Bio import SeqIO
genome=SeqIO.parse(<genomeFile>,"fasta")
for record in genome:
   entry=[record.id,1,len(record)]
   intervals.append(entry)

out=open(<fullScaffolds.bed>,"w")
for e in intervals:
  out.write(e[0]+"\t"+str(e[1])+"\t"+str(e[2])+"\n")

```
```shell
subtractBed -a fullScaffolds.bed -b gene.gff > intergenic.bed
```

Then, for all types, get alignment with:

```shell
intersectBed -wo -a depthBedFile -b gffFile > outputFile
```
The wo option outputs the number of base pairs involved in the overlap. I'll then compare the number of bp overlapped to the total number of bp in the gff file.
