#!/usr/bin/env Rscript
library(ggplot2)

bedFile=read.csv("results/finalAssemblies_expandedSubTree/finalAssemblies_expandedSubTree_alignmentDepth_gaps.bed",sep="\t", header = F, col.names = c("scaffold","start","end"))
bedFile$length=bedFile$end-bedFile$start

plot.new()
plot <- ggplot(bedFile,aes(x=length))
plot+geom_histogram(bins = 500)+scale_y_log10()+
  geom_vline(xintercept = mean(bedFile$length))+
  labs(title="gaps")

