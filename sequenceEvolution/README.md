###First aim: find heliconius-specific deletions

use this [script](findRegionsExclusivelyInGroup.slurm), which relies on the halTools script findRegionsExclusivelyInGroup.    
How I use it: findRegionsExclusivelyInGroup halFile Bombyx  all non-heliconius genomes
- I use Bombyx as reference because that genome is well-annotated, and if I'm really looking for all
regions that are in ALL non-heliconians it doesn't matter which one I use as reference.     

Problems so far:
- Don't know how to visualize! Annoying. It almost works with geneious, but transferring annotations isn't working.
- Need to group separate regions that are spaced by just a few nucleotides, but not entirely sure *how* close they need to be.

Interesting things so far:
- when looking at the intermediate hal file (most distant outgroup is monarch), there is one gene that is completely deleted in ALL heliconius! DPOGS210663, which is identified as similar to HSP70...should read up more about this.

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
