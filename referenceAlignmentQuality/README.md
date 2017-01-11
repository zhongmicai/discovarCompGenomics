#Assess alignment of reference genomes to conspecific DISCOVAR genomes

##Factors to consider: 
- Depth (Percent of genome aligned)
	- Break this down into genic and non-genic regions
- Compare scaffold/contig breaks in DISCOVAR vs Reference


##Scripts to run###

For coverage depth

Shell script:

```bash
#Find the alignment depth of the relevant comparisons
halAlignmentDepth --targetGenomes HmelDisco --outWiggle coverageDepth_melpomene_refToDisco.wig <halPath> HmelRef 
halAlignmentDepth --targetGenomes HmelRef --outWiggle coverageDepth_melpomene_discoToRef.wig <halPath> HmelDisco
halAlignmentDepth --targetGenomes HeraDisco --outWiggle coverageDepth_erato_refToDisco.wig <halPath> HeraRef 
halAlignmentDepth --targetGenomes HeraRef --outWiggle coverageDepth_erato_discoToRef.wig <halPath> HeraDisco

#take the coverage files and make bed files with all the covered regions. The numbers at the end are minCov, maxGap, minLength of blocks
wigToBed.py coverageDepth_melpomene_refToDisco.wig coverageDepth_melpomene_refToDisco.bed 1 0 1
wigToBed.py coverageDepth_melpomene_discoToRef.wig coverageDepth_melpomene_refToDisco.bed 1 0 1
wigToBed.py coverageDepth_erato_refToDisco.wig coverageDepth_melpomene_refToDisco.bed 1 0 1
wigToBed.py coverageDepth_melpomene_refToDisco.wig coverageDepth_melpomene_refToDisco.bed 1 0 1

#
```
R Script:

```R
#!/usr/bin/env Rscript
#install.packages("ggthemes")
library(ggplot2)
library(ggthemes)

=read.csv("results/expandedSubTree_Etal_depth.bed",sep="\t", header = F, col.names = c("scaffold","start","end"))
subTree$length=subTree$end-subTree$starts
subTreeNumBlocks=length(subTree$length)
subTreeTotal=sum(subTree$length)
'''


