#!/usr/bin/env Rscript
#install.packages("ggthemes")
library(ggplot2)
library(ggthemes)

subTree=read.csv("results/expandedSubTree_Etal_depth.bed",sep="\t", header = F, col.names = c("scaffold","start","end"))
subTree$length=subTree$end-subTree$start
subTreeNumBlocks=length(subTree$length)
subTreeTotal=sum(subTree$length)
subTree$alignment <- rep("subTree",nrow(subTree))

subTree_gap10=read.csv("results/expandedSubTree_Etal_gap10.bed",sep="\t", header = F, col.names = c("scaffold","start","end"))
subTree_gap10$length=subTree_gap10$end-subTree_gap10$start
subTree_gap10NumBlocks=length(subTree_gap10$length)
subTree_gap10Total=sum(subTree_gap10$length)
subTree_gap10$alignment <- rep("subTree_gap10",nrow(subTree_gap10))

smallAlign=read.csv("results/smallAligntrial_Hmel_Etal_alignmentDepth.bed",sep="\t", header = F, col.names = c("scaffold","start","end"))
smallAlign$length=smallAlign$end-smallAlign$start
smallAlignNumBlocks=length(smallAlign$length)
smallAlignTotal=sum(smallAlign$length)
smallAlign$alignment <- rep("smallAlign",nrow(smallAlign))

smallAlign_gap10=read.csv("results/smallAligntrial_Hmel_Etal_alignmentDepth_gap10.wig",sep="\t", header = F, col.names = c("scaffold","start","end"))
smallAlign_gap10$length=smallAlign_gap10$end-smallAlign_gap10$start
smallAlign_gap10NumBlocks=length(smallAlign_gap10$length)
smallAlign_gap10Total=sum(smallAlign_gap10$length)
smallAlign_gap10$alignment <- rep("smallAlign_gap10",nrow(smallAlign_gap10))

allAligns <- rbind(subTree, subTree_gap10, smallAlign, smallAlign_gap10)
notTinyAligns <- subset(allAligns, length>5000)


df <- data.frame(
  x = c(10000, 10000, 10000),
  y = c(100000,50000,30000),
  text=c("Five-way alignment blocks","Two=way alignment blocks", "Five-way alignment blocks with up to 10bp gap")
  #text = c(paste0("Five-way alignment: blocks= ", subTreeNumBlocks,", length= ",subTreeTotal), paste0("Two-way alignment: blocks= ", smallAlignNumBlocks,", length= ",smallAlignTotal))
  )

lengths <- ggplot(allAligns,aes(x=length))
lengths + geom_freqpoly(aes(y = ..count.., col=alignment),stat = "bin", bins=500) +
  #geom_text(data=df,aes(x=x,y=y,label = text, col=text))+
  scale_color_manual(values=c("blue","navyblue", "red","magenta"))+
  scale_y_log10()+
  xlim(0,25000)+
  theme_stata()+
  labs(title="blocks")

boxes <- ggplot(notTinyAligns,aes(factor(alignment),length))#, y=c(rep(1,length(subTree$length)))))
boxes+geom_boxplot(aes(fill=alignment))+
  scale_fill_manual(values=c("blue","navyblue", "red","magenta"))+
  scale_y_log10()



