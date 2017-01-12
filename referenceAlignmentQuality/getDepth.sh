#!/bin/bash

halPath=$1

#Find the alignment depth of the relevant comparisons
#This is hard-coded for the intermediate alignment downloaded 1/1//16, and to look only at 1:1 comparisons of the reference genomes to their conspecific DISCOVAR genomes.
sbatch getAlignmentDepth.slurm intermediateAlignment.hal HmelDisco HmelRef coverageDepth_melpomene_refToDisco.wig coverageDepth_melpomene_refToDisco.bed
sbatch getAlignmentDepth.slurm intermediateAlignment.hal HmelRef HmelDisco coverageDepth_melpomene_discoToRef.wig coverageDepth_melpomene_refToDisco.bed
sbatch getAlignmentDepth.slurm intermediateAlignment.hal HeraDisco HeraRef coverageDepth_erato_refToDisco.wig coverageDepth_melpomene_refToDisco.bed
sbatch getAlignmentDepth.slurm intermediateAlignment.hal HeraRef HeraDisco coverageDepth_erato_discoToRef.wig coverageDepth_melpomene_refToDisco.bed
