##Conservation Analysis

In this section, I'l look for regions that have either accelerated evolution or become more conserved in *Heliconius*. The internal tool halPhloP should do this. halyPhyloP needs a neutral model of evolutions, and it used 4D sites to generate it. That script is halPhyloPTrain.py, and it takes a bed file to find the 4d sites. The bed file has to be of only exons and bed12, so I translated the Bmor_cDNA_exon.gff file with the line:
```shell
awk '{OFS="\t"} {print $1,$4,$5,"exon","0",$7,$4,$5,"0","1",$5-$4,"0"}' Bmor_cDNA_exon.gff > Bmor_cDNA_exon.bed
```
Then, I created phyloPTrain.slurm, which runs the command :
```
halPhyloPTrain.py --numProc 4 --noAncestors $halFile $refGenome $bedFile $outFile
```
and called it with
```
sbatch phyloPTrain.slurm /n/regal/mallet_lab/edelman/18Genomes/results/DiscovarCompGenomics/data/finalAssemblies_highQual_1kbFilter_161101.hal Bmor Bmor_cDNA_exon.bed fullPhylogeny_neutralModel.mod
```
Should think about the benefits of doing this with monarch as the reference...will it make a difference??