#Assess alignment of reference genomes to conspecific DISCOVAR genomes

###This should give us some idea of how well the alignment worked, and possibly some idea of how good the DISCOVAR genomes are. In the case of melpomene, the DISCOVAR and reference individual are cousins. In the case of erato, I'm not sure how closely related the two individuals are, but we should find that out, and they are certainly more distantly related than the melpomene ones. 

##Factors to consider: 
- Depth (Percent of genome aligned)
	- Break this down into genic and non-genic regions
- Compare scaffold/contig breaks in DISCOVAR vs Reference


##Scripts to run###

For coverage depth

Shell script: (./getDepth.sh "hard-coded controller for 1:1 comparisons of reference genomes to conspecific DISCOVAR genomes")

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


