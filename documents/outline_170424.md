DISCOVAR paper plan outline

## Low-cost sequencing of high phyologenetic density reveals X

###Introduction

- Adaptive radiation, especially key innovations
- Heliconius, and why they are a good system
	- Lot of ecology pertaining to what makes Heliconius unique, but most molecular studies focused on color pattern differences within the clades
	- Using whole-genome data, we can try to find the genetic basis of the traits we already know about (pollen-feeding, intelligence, longevity, etc), and possibly find previously unknown pathways that are important
	- touch on hybridization - known to be occurring in modern times, but how prevalent has it been throughout the radiation?
- advantages of whole genome sequencing
	- identify large-scale chromosomal rearrangements
	- use as a resource for population-level work
- advantages of not using reference alignment
  - doesn't bias us towards things we already know about
- DISCOVAR, and why it is a good tool
- ProgressiveCactus, and why it is a good tool

###Materials and Methods
- DNA extraction
(supplementary – get extraction protocols from all labs involved)
- Sequencing and Assembly
- DISCOVAR
- TAGC scaffolding
- Multiple alignment (ProgressiveCactus)
- Data Analysis

###Results

**Table 1 - Figure 2 could all go in supplementary**

Table 1: Assembly statistics for all species
- from Davey: N50, length, largest scaffold, BUSCO

Figure 1: Example of a DISCOVAR assembly
- Use the figure design from Lepbase
- *use lepbase*

Figure S1: All of the DISCOVAR assemblies in Lepbase format
- *see above*

Figure 2: Alignment metrics
  - pairwise alignments
  - "fully aligned" sites
  - sensitivity to input order of tree


Figure 3: Phylogenetics - this could be the pretty butterfly figure
  - *CDS tree, 4D sites tree, non-coding tree, UCE tree, non-UCE tree, full coverage sites*
	- Follow up on why CDS tree is different from "fully aligned sites" tree
	- Evidence for hybridization towards the base of the phylogeny

Figure 4: Genome structure conservation
- We know there are no major rearrangements between melpomene and erato, but look across whole tree for major changes
	- Not sure if we have the power to really detect these within Heliconius, and John Davey already found "fusion points" where the 31
		Eueides chromosomes became 21 in Heliconius, but it could be nice to see if we can also confirm them with our data
	- maybe an argument could be made here that conservation of synteny allows hybridization? Or common hybridization would result in
	 	maintenance of synteny?
- If it works, this would be a good argument for creating independent genome assemblies as opposed to aligning resequence data


Figures 5-7: strategies for identifying Heliconius-specific genetic features

Novel genes that are highly conserved within Heliconius
- We're going to have a lot of regions that are heliconius-specific, just because they're more closely related, but if we intersect this list with a list of the most conserved genes maybe we'll find something interesting.

Sequence conservation changes
- Find regions that are highly conserved within Heliconius as opposed to the rest of the lepidoptera
- Find regions that were highly conserved among the lepidoptera but selection has relaxed in Heliconius
- both genic and intergenic

Deleted regions in Heliconius
- Find regions present in all lepidoptera, but absent in Heliconius

Figure 8: Make sense of the candidates
- Maybe for each of figures 5-7 we can have some top hits or top GO-categories. For Figure 8 we can see how they relate to one another: are they telling different parts of the same story or three different stories? Are there many more strong candidates among one type of identification method than the others?

###Conclusions
We used whole-genome assembly and alignment of several species within a single genus to identify genetic factors that led to their adaptive radiation. This includes traits like X, Y, and Z, as well as a general pattern of hybridization (maybe). We identified candidate genes for further study, and they look particularly promising for X. ...
