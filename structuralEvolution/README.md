update 2/10/17
extracted MAF alignments of each melpomene scaffold with the command:

```shell
scaffs=$(cat HmelGE1kbScafNames.txt);
for scaf in $scaffs;
do echo $scaf
sbatch halToMaf_scafs.slurm HmelRef <halFile> finalAssemblies_HmelRef_$scaf.maf $scaf
sleep 1
done
```

This referenced the script [halToMaf_scafs.slurm](melpomeneScaffoldMafs/halToMaf_scafs.slurm),
which calls the halTools tool hal2maf. 
