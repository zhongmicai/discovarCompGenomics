# DiscovarCompGenomics

Contents of folders:

Data
**Turns out the data is too big to easily put onto github...these files can be found on odyssey at /n/regal/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/data**

  - Hmel2.gff is the gff file with gene calls from the Hmel2 January 2016 release
  - subTree_18Genomes.hal is the output from progressive cactus of the tree with
      H. melpomene reference, H. melpomene discovar, H. cydno, H. timareta, and H. numata
  - subTree_18Genomes_Hmel201001.maf is the output of extractMaf.sh subTree_18Genomes.hal HmelRef Hmel201001
      i.e. it is a maf file with the alignment referenced on HmelRef (Hmel2) for the scaffold Hmel201001
  - progressiveCactusTree.txt is a newick format tree of the FULL DISCOVAR alignment which has not yet been run

Scripts
  - coveragePipeline.sh is a wrapper script to run the data processing from the hal output of cactus to fasta files.
      currently (June 30) it only gets as far as finding blocks of high alignment depth.
  - coveragePipeline.slurm is a command to run coveragePipeline.sh on the odyssey cluster
  - extractMaf.sh is a script that takes a hal alignment, a reference genome, and a sequence name and outputs a maf file.
      It uses the command hal2maf --noAncestors --refGenome $refGenome --refSequence $scafName $halFile $mafFile
  - mafToFasta.py takes a maf file and outputs the aligned sequences in fasta format
  - wigToBed.py is a script that takes a wig file, and extracts regions of interest based on depth,length, and gap length
      outputs a bed file indicating each high quality region

Documents
  - DISCOVAR paper plan outline.docx is the original plan I sent on 3/21/16
  - DISCOVAR_DC_outline.docx is the modified plan we came up with after meeting; this was uploaded on 4/17/16

Results
  *Nothing in here yet! a place holder for results when they come.
