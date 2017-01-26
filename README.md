# DiscovarCompGenomics

## Project to study Evolutionary dynamics across *Heliconius* using whole-genome DISCOVAR sequencing and progressiveCactus alignment


### Assembly quality
Goals:   
1) create Lepbase-esque circle plot of all genomes used in this study
  - [github repository](https://github.com/rjchallis/assembly-stats)  
  
2) evaluate completeness of DISCOVAR melpomene by comparing to Hmel2  

  - *alignment depth* 
  - *break into genic vs intergenic*
  - *find syntenic breaks DISCOVAR vs reference* (haltools syntenyRates.py)
  - *differences in alignment of erato ref and erato DISCO to Hmel2. Also melpomene ref and Disco to Erato ref*

### Alignment Quality



### Phylogeny

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

