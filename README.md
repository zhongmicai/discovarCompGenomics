# DiscovarCompGenomics

## Project to study Evolutionary dynamics across *Heliconius* using whole-genome DISCOVAR sequencing and progressiveCactus alignment


### [Assembly quality](AssemblyStats)
Goals:  

1) create Lepbase-esque circle plot of all genomes used in this study
  - [github repository](https://github.com/rjchallis/assembly-stats)  
  
2) [evaluate completeness of DISCOVAR melpomene by comparing to Hmel2](referenceAlignmentQuality)  

  - *alignment depth* 
  - *break into genic vs intergenic*
  - *find syntenic breaks DISCOVAR vs reference* (haltools syntenyRates.py)
  - *differences in alignment of erato ref and erato DISCO to Hmel2. Also melpomene ref and Disco to Erato ref*

### [Alignment Quality](AlignmentStats)
Ideas:    

  - pairwise alignment depth
  - all-Genomes alignment depth (which genome should it be based on?? melpomene?)
  - GO categories of genes that are conserved throughout (do we get all the housekeeping genes?)
  - exonic vs intronic vs intergenic


### [Phylogeny](phylogenies)

  - Paul
  - *Generate CDS tree, 4D sites tree, non-coding tree, UCE tree, non-UCE tree, full coverage sites* 
  
    - *>= 10 starting trees *
    - *bootstraps until MRE returns converged result*


### [Genome structure conservation](structuralEvolution)

  - Rebecca
  - show syntenic breaks
  - This might be a messy figure with so many genomes/comparisons, so maybe we can just look for major differences and think of a way to  display them well

### [Genome sequence conservation](sequenceEvolution)  

  - Nate
  - Deletions   
     - *findRegionsExclusivelyInGroup*
     - both deletions in heliconius and deletions in outgroups (duplications in heliconius)
     - Use cDNA-based annotations to analyze exons vs introns vs intergenic deletions
     - GO terms
  - sequence conservation  
    - *halphylop*
    - proportion conserved across group
    - proportion conserved within heliconius
    - proportion conserved within clades
    - Motifs???
    
  - For statistical interest, compare random subsets *look into halRandGen*
  - for all groups, look at coding vs non-coding
  
  - Evolutionary rates/Positive Selection
    - *see above* 
    - dN/dS along all branches	
    - substitution rate
    - comparison as above

### Candidate loci exploration


