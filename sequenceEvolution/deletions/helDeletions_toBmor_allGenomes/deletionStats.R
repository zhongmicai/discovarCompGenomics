#heliconiusDeletions.R

#!/usr/bin/env Rscript
library(ggplot2)
library(ggthemes)
library(rphast)


###Exons###

###Aligned regions###
alignedExons=read.csv("sequenceEvolution/deletions/helDeletions_toBmor_allGenomes/alignedRegionsOutgroups_exons.bed",sep="\t", header = F, col.names = c("alignmentScaf","alignmentStart","alignmentEnd","featureScaf","featureOrigin","type","featureStart","featureEnd","dot","strand","dot2","ID","overlapLength"))
geneName=c()
for (line in alignedExons$ID){
  geneName=c(geneName,strsplit(as.character(line),split="[.]|=|;")[[1]][2])
}
alignedExons$geneName=geneName
alignedExonsBases=sum(alignedExons$overlapLength)
numAlignedExonBases=length(alignedExons$overlapLength)
exonsAligned=length(unique(alignedExons$ID))
alignedExonGenesHit=length(unique(alignedExons$geneName))
alignedExonGenes=unique(alignedExons$geneName)


###Deletions###
exonDeletions=read.csv("sequenceEvolution/deletions/helDeletions_toBmor_allGenomes/heliconiusDeletions_bombyxRef_allGenomes_exon.bed",sep="\t", header = F, col.names = c("deletionScaf","deletionStart","deletionEnd","featureScaf","featureOrigin","type","featureStart","featureEnd","dot","strand","dot2","ID","overlapLength"))
geneName=c()
for (line in exonDeletions$ID){
  geneName=c(geneName,strsplit(as.character(line),split="[.]|=|;")[[1]][2])
}
exonDeletions$geneName=geneName
exonDeletedBases=sum(exonDeletions$overlapLength)
numExonDeletions=length(exonDeletions$overlapLength)
exonsHit=length(unique(exonDeletions$ID))
exonGenesHit=length(unique(exonDeletions$geneName))
exonGenes=unique(exonDeletions$geneName)
exonPerGeneNumHits=c()
exonPerGeneNumBases=c()
for (gene in exonGenes){
  subFrame=subset(exonDeletions, geneName==gene)
  exonPerGeneNumHits =c(exonPerGeneNumHits,nrow(subFrame))
  exonPerGeneNumBases =c(exonPerGeneNumBases,sum(subFrame$overlapLength))
}
exonDeletionStats=data.frame("gene"=exonGenes,"numDeletions"=exonPerGeneNumHits,"numBases"=exonPerGeneNumBases)
exonMultiHits=subset(exonDeletionStats,numDeletions>1)

df=cbind(c(exonGenes, rep("",length(alignedExonGenes)- length(exonGenes))),c(as.character(exonMultiHits$gene), rep("",length(alignedExonGenes)- length(exonMultiHits$gene))), alignedExonGenes)
write.csv(df,"exonDeletions.csv")

###Introns###

###Aligned regions###
alignedIntrons=read.csv("sequenceEvolution/deletions/helDeletions_toBmor_allGenomes/alignedRegionsOutgroups_introns.bed",sep="\t", header = F, col.names = c("alignmentScaf","alignmentStart","alignmentEnd","featureScaf","featureOrigin","type","featureStart","featureEnd","dot","strand","dot2","ID","overlapLength"))
geneName=c()
for (line in alignedIntrons$ID){
  geneName=c(geneName,strsplit(as.character(line),split="[.]|=|;")[[1]][2])
}
alignedIntrons$geneName=geneName
alignedIntronsBases=sum(alignedIntrons$overlapLength)
numAlignedIntronsBases=length(alignedIntrons$overlapLength)
intronsAligned=length(unique(alignedIntrons$ID))
alignedIntronGenesHit=length(unique(alignedIntrons$geneName))
alignedIntronGenes=unique(alignedIntrons$geneName)

###Deletions
intronDeletions=read.csv("sequenceEvolution/deletions/helDeletions_toBmor_allGenomes/heliconiusDeletions_bombyxRef_allGenomes_intron.bed",sep="\t", header = F, col.names = c("deletionScaf","deletionStart","deletionEnd","featureScaf","featureOrigin","type","featureStart","featureEnd","dot","strand","dot2","ID","overlapLength"))
geneName=c()
for (line in intronDeletions$ID){
  geneName=c(geneName,strsplit(as.character(line),split="[.]|=|;")[[1]][2])
}
intronDeletions$geneName=geneName
intronDeletedBases=sum(intronDeletions$overlapLength)
numIntronDeletions=length(intronDeletions$overlapLength)
intronsHit=length(unique(intronDeletions$ID))
intronGenesHit=length(unique(intronDeletions$geneName))
intronGenes=unique(intronDeletions$geneName)
intronPerGeneNumHits=c()
intronPerGeneNumBases=c()
for (gene in intronGenes){
  subFrame=subset(intronDeletions, geneName==gene)
  intronPerGeneNumHits =c(intronPerGeneNumHits,nrow(subFrame))
  intronPerGeneNumBases =c(intronPerGeneNumBases,sum(subFrame$overlapLength))
}
intronDeletionStats=data.frame("gene"=intronGenes,"numDeletions"=intronPerGeneNumHits,"numBases"=intronPerGeneNumBases)
intronMultiHits=subset(intronDeletionStats,numDeletions>1)

df=cbind(c(intronGenes, rep("",length(alignedIntronGenes)- length(intronGenes))),c(as.character(intronMultiHits$gene), rep("",length(alignedIntronGenes)- length(intronMultiHits$gene))), alignedIntronGenes)
write.csv(df,"intronDeletions.csv")

###Intergenic###
###Aligned Regions###
alignedIntergenic=read.csv("sequenceEvolution/deletions/helDeletions_toBmor_allGenomes/alignedRegionsOutgroups_intergenic.bed",sep="\t", header = F, col.names = c("deletionScaf","deletionStart","deletionEnd","featureScaf","featureStart","featureEnd","overlapLength"))
intergenicAlignedBases=sum(alignedIntergenic$overlapLength)

###Deletions###
intergenicDeletions=read.csv("sequenceEvolution/deletions/helDeletions_toBmor_allGenomes/heliconiusDeletions_bombyxRef_allGenomes_intergenic.bed",sep="\t", header = F, col.names = c("deletionScaf","deletionStart","deletionEnd","featureScaf","featureStart","featureEnd","overlapLength"))
intergenicDeletedBases=sum(intergenicDeletions$overlapLength)
numIntergenicDeletions=length(intergenicDeletions$overlapLength)
largeIGDels=subset(intergenicDeletions, overlapLength >22)


####Table###

deletionStats=data.frame("stats"=c("Deleted Bases","Pct of Aligned Bases", "Genes Hit","Pct of Aligned Genes","Genes Hit Multiple Times","Pct of Aligned Genes","Average Length of Deletion","Average Number Deleted Bases per Gene"), 
                         "exons"=c(exonDeletedBases,(exonDeletedBases/alignedExonsBases)*100,exonGenesHit,(exonGenesHit/alignedExonGenesHit)*100,length(which(exonDeletionStats$numDeletions >1)),(length(which(exonDeletionStats$numDeletions >1))/alignedExonGenesHit)*100,mean(exonDeletions$overlapLength),mean(exonPerGeneNumBases)),
                         "introns"=c(intronDeletedBases,(intronDeletedBases/alignedIntronsBases)*100,intronGenesHit,(intronGenesHit/alignedIntronGenesHit)*100,length(which(intronDeletionStats$numDeletions >1)),(length(which(intronDeletionStats$numDeletions >1))/alignedIntronGenesHit)*100,mean(intronDeletions$overlapLength),mean(intronPerGeneNumBases)),
                         "intergenic"=c(intergenicDeletedBases,(intergenicDeletedBases/intergenicAlignedBases)*100,NA,NA,NA,NA,NA,NA))









