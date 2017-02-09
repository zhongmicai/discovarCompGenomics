##Conservation Analysis

In this section, I'l look for regions that have either accelerated evolution or become more conserved in *Heliconius*. The internal tool halPhloP should do this. halyPhyloP needs a neutral model of evolutions, and it used 4D sites to generate it. That script is halPhyloPTrain.py, and it takes a bed file to find the 4d sites. The bed file has to be of only exons and bed12, so I translated the Bmor_cDNA_exon.gff file with the line:
```shell
awk '{OFS="\t"} {print $1,$4,$5,$9,"0",$7,$4,$5,"0","1",$5-$4,"0"}' Bmor_cDNA_exon.gff > Bmor_cDNA_exon.bed
```
Then, I created phyloPTrain.slurm, which runs the command :
```
halPhyloPTrain.py --numProc 4 --noAncestors $halFile $refGenome $bedFile $outFile
```
and called it with
```
sbatch phyloPTrain.slurm /n/regal/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/data/finalAssemblies_highQual_1kbFilter_161101.hal Bmor Bmor_cDNA_exon.bed fullPhylogeny_neutralModel_BmorRef.mod
```
Did the same thing with monarch and bicyclus...not sure if it will make a difference but it might because we'll presumably get more regions aligned with more species? Need to check this with halAlignmentDepth.

Update 1/31
The halPhyloP script takes a VERY long time, so I'll try to parallelize it. The easiest way, I think, is to run halPhyloP for each chromosome (or at least small groups...there are 4782 bombyx chromosomes). Made a txt file with a list of all the bombyx scaffold names with
```
grep NW /n/holylfs/INTERNAL_REPOS/PCACTUS/edelman/genomes/1kbFilter/Bombyx_mori_ASM15162v1_-_scaffolds.fa_1kb.fa | awk '{print $1}'> bombyxGE1kbScafNames.txt
sed -i.bak 's/>/ /' bombyxGE1kbScafNames.txt
rm -f bombyxGE1kbScafNames.txt.bak
```
