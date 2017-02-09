scafList=$(cat $1)
halFile=$2
refGenome=$3
mod=$4
overallOut=$5

for scaf in $scafList;
do echo $scaf;
sbatch phyloP_specSeqs.slurm $halFile $refGenome $mod $scaf $overallOut ;
sleep 1;
done
