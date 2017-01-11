DISCOVAR paper plan outline

## Evolutionary dynamics across Heliconius

###Introduction 

- Adaptive radiation 
	- specifically, the idea of key innovations
- advantages of whole genome sequencing
	- identify large-scale chromosomal rearrangements
	- use as a resource for population-level work
- Heliconius, and why they are a good system
- DISCOVAR, and why it is a good tool

###Materials and Methods
- DNA extraction
(supplementary – get extraction protocols from all labs involved)
- Sequencing and Assembly
- DISCOVAR – ask Wiesenfeld or Jaffe? Alternatively site forum?
- Possibly TAGC scaffolding
- haplomerger
- Multiple alignment
- progressiveCactus
- use Kozak tree to prior topology
- Data Analysis

###Results
Table 1: Assembly statistics for all species
	- from Davey: N50, length, largest scaffold, BUSCO
  - *DATA AVAILABLE MAKE FIGURE*

Figure 1: Example of a DISCOVAR assembly
	- Use the figure design from Lepbase
  - *use lepbase*

Figure S1: All of the DISCOVAR assemblies in Lepbase format
  - *see above*

Figure S2: Alignment quality control for H. melpomene DISCO to Hmel2
  - *alignment depth* 
  - *break into genic vs intergenic*
  - *find syntenic breaks DISCOVAR vs reference* (haltools syntenyRates.py)

Figure 2: Whole genome tree
  - *Generate CDS tree, 4D sites tree, non-coding tree, UCE tree, non-UCE tree, full coverage sites*
    - *>= 10 starting trees *
    - *bootstraps until MRE returns converged result*

Figure 3: Genome structure conservation
- *Talk to Rebecca*
- show syntenic breaks
- This might be a messy figure with so many genomes/comparisons, so maybe we can just look for major differences and think of a way to display them well

Figure 4: Genome sequence conservation  
- *halphylop*
- *findRegionsExclusivelyInGroup*
- proportion conserved across group
- proportion conserved within heliconius
- proportion conserved within clades
- For statistical interest, compare random subsets *look into halRandGen*
- for all groups, look at coding vs non-coding

Figure 5: Evolutionary rates/Positive Selection
- *see above* 
- dN/dS along all branches	
- substitution rate
- comparison as above

Figure 6-X: Candidate loci exploration



###Next paper? - fine-scale phylogeny and introgression

- use SOM/HMM models or mosquito paper technique to get trees along genome 
-	How old are the trees?
-	How much of the genome does each tree represent?
-	How large are the genomic blocks that represent each tree type?


