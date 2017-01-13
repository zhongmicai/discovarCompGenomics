#!/bin/bash

halPath=$1

#Find the alignment depth of the relevant comparisons
#This is hard-coded to look only at 1:1 comparisons of the reference genomes to their conspecific DISCOVAR genomes.
sbatch getAlignmentDepth.slurm $halPath HmelRef HmelDisco coverageDepth_melpomene_refToDisco.wig coverageDepth_melpomene_refToDisco.bed
sbatch getAlignmentDepth.slurm $halPath HmelDisco HmelRef coverageDepth_melpomene_discoToRef.wig coverageDepth_melpomene_discoToRef.bed
sbatch getAlignmentDepth.slurm $halPath HeraRef HeraDisco coverageDepth_erato_refToDisco.wig coverageDepth_melpomene_refToDisco.bed
sbatch getAlignmentDepth.slurm $halPath HeraDisco HeraRef coverageDepth_erato_discoToRef.wig coverageDepth_melpomene_discoToRef.bed
