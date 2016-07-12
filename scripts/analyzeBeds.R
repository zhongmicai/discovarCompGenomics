#!/usr/bin/env Rscript
library(ggplot2)

bedFile=read.csv("results/_greaterThan1.bed",sep="\t", header = F, col.names = c("scaffold","start","end"))
bedFile$length=bedFile$end-bedFile$start


plot <- ggplot(bedFile,aes(x=length))
plot+geom_histogram(bins = 500)+scale_y_log10()+
  geom_vline(xintercept = mean(bedFile$length))+
  labs(title="greaterThan1")

