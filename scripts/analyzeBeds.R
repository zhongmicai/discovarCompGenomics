#!/usr/bin/env Rscript
#install.packages("ggthemes")
library(ggplot2)
library(ggthemes)

bedFile=read.csv("data/finalData_subTree/finalAssemblies_expandedSubTree_alignmentDepth_gaps.bed",sep="\t", header = F, col.names = c("scaffold","start","end"))
bedFile$length=bedFile$end-bedFile$start


lengths <- ggplot(bedFile,aes(x=length))
lengths + geom_freqpoly(aes(y = ..count..),stat = "bin", bins=500) +scale_y_log10()+
  xlim(0,20000)+
  theme_stata()+
  labs(title="Gaps")


