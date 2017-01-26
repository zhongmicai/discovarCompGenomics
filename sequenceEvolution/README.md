###First aim: find heliconius-specific deletions

use this [script](deletions/findDeletions.slurm), which relies on the halTools script findRegionsExclusivelyInGroup.    
How I use it: 
```shell
findDeletions.slurm halFile Bombyx  all non-heliconius genomes
```
- I use Bombyx as reference because that genome is well-annotated, and if I'm really looking for all
regions that are in ALL non-heliconians it doesn't matter which one I use as reference.     

Problems so far:
- Don't know how to visualize! Annoying. It almost works with geneious, but transferring annotations isn't working.
- Need to group separate regions that are spaced by just a few nucleotides, but not entirely sure *how* close they need to be.

Interesting things so far:
- when looking at the intermediate hal file (most distant outgroup is monarch), there is one gene that is completely deleted in ALL heliconius! DPOGS210663, which is identified as similar to HSP70...should read up more about this.
- Looked at the full HAL file, and found only 16 regions >= 100bp that were present in ALL outgroups and absent in ALL heliconius. Went through and took IGV screenshots, as well as gathering info about what's known about the regions, though I still need to organize that in some reasonable way. Turns out all of them are in protein-coding regions, and for the most part they cover a full exon. In two cases, multiple exons were deleted, and in one case 5 (!) exons were deleted!! Definitely need to follow that one up. Should do some sort of more general analysis, like bedtools intersect of all the deletions to the bombyx gff file.'

scripts to try:


hal4dExtract    
halAlignedExtract   
halphylop   
halPhyloPTrain  
halRandGen    
halSnps   
halTreeMutations.py   
halTreeNIBackground.py    
halTreeNIConservation.py    
halTreeNITurnover.py    
halTreePhyloP.py    
constraintTurnover    
halContiguousRegions.py   
syntenyRates.py   
