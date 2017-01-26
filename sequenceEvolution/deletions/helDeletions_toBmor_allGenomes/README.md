##Regarding the gff gene file

I had some issues finding the right gene set, as well as figuring out why my genome file had scaffolds named differently than the ones on silkdb and ensembl. I finally found the naming conventions, and it's a simple transfer from one to the other. There are actually two transfers that must be made: one from scaffolds starting with NW (how our genome is named, and representing the refseq archive name) to ones starting with DF or BABH (the genome center name) these synonyms are found in this [file](bombyx_NWscaffold_names.txt) that I downloaded from [ncbi](ftp://ftp.ncbi.nih.gov/genomes/Bombyx_mori/scaffold_names). The second is from the DF or BABH names to Bm_scaf names. I got that translation from [kaikobase](http://sgp.dna.affrc.go.jp/KAIKObase/keyword_search.php?keyword=&field=all&chr_or_scaf=&chr_id=all&chr_start=&chr_end=&scaf_id=&scaf_start=&scaf_end=&graph_view=on&page=16). I used those tables to create a [python script](bombyxScaffoldTransfer.py) that changes the names of a gff file to the names associated with our genome. There were a few different options for gff files, but I think the best one is from [this paper](https://www.ncbi.nlm.nih.gov/pubmed/23821615), which uses cDNA and ESTs to find genes. There are two gene sets from that paper, one that relies on cDNA and the other on ESTs. They call them geneset A and B. There is also a geneset C, but those are not mapped to the genome so we won't have them. I combined genesets A and B into a single gff file, called [Bmor_cDNA.gff](Bmor_cDNA.gff), and then used my script to change the scaffold names to produce [Bmor_cDNA.gff_transfer.gff](Bmor_cDNA.gff_transfer.gff). I think this is the best data. I had already done some analysis with a gff file I got from lepbase, and I put all those files in the directory oldGFF.


###deletion analysis

Created the file heliconiusDeletions_bombyxRef_allGenomes.bed with the script findDeletions.slurm, with searching for regions that aligned in all lepidoptera except heliconius:

```shell
findDeletions.slurm ../../data/finalAssemblies_highQual_1kbFilter_161101.hal Bmor Pxyl,Bmor,Lacc,Ppol,Dple,Bany,Mcin,Avan,Etal > heliconiusDeletions_bombyxRef_allGenomes.bed
```

Then, divided deletions out among exons, introns, and intergenic regions. Interestingly, they are overwhelmingly in introns, at least by number of distinct deletions (1488 exonic, 6212 intergenic, 11038 intronic). Got these subsets with:

```shell
intersectBed -wo -a heliconiusDeletions_bombyxRef_allGenomes.bed -b <gff_of_interest> > output
```

To do:  
1) need to know how many intronic, exonic, intergenic bases were aligned among non-heliconius lepidopterans in the first place.   
2) get data on number of nucleotides deleted, not just number of deletions  
3) get info on genes hit by deletions