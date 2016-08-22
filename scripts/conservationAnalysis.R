#ConservationAnalysis.R
install.packages("seqinr")
require(rphast)
library(ggplot2)
library(seqinr)


# coverage depth ----------------------------------------------------------
#make a pie chart of coverage depth
depths=read.delim("results/subTreeAlignmentDepth_coverageDistribution.tsv", header=F,col.names=c("coverage","values"))
depths$coverage=as.character(depths$coverage)
depths$coverage=factor(depths$coverage, levels=c("0","1","2","3","4"))


depth <- ggplot(data=depths, aes(x=coverage,y=values/sum(depths$values),fill=coverage))+
  geom_bar(width=1, stat="identity")+geom_label(label=round(depths$values/sum(depths$values),3)) +
  labs(title="Alignment Coverage",x="Coverage Depth",y="Fraction of Alignment")
depth
pie <- depth+coord_polar("y",start=0)
pie

#for exons only:
exonDepth=read.delim("results/exonDepth.bed", header=F,col.names=c("chrom","start","end","id","coverage"))
depthCounts=data.frame(cov=c("0","1","2","3","4"),depth=c(length(which(exonDepth$coverage==0)),length(which(exonDepth$coverage==1)),
                                                length(which(exonDepth$coverage==2)),length(which(exonDepth$coverage==3)),
                                                length(which(exonDepth$coverage==4))))
exonGraph<- ggplot(data=depthCounts, aes(x=cov,y=depth/sum(depthCounts$depth),fill=cov))+
  geom_bar(width=1, stat="identity")+geom_label(label=round(depthCounts$depth/sum(depthCounts$depth),3)) +
  labs(title="Exon Coverage",x="Coverage Depth",y="Fraction of Alignment")
exonGraph


# find and graph most conserved elements ----------------------------------

mafFile="data/subTree_18Genomes_Hmel201009.maf"
gffFile="data/Hmel2.gff"
scaffoldName="Hmel201009"
referenceName="HmelRef"
newickTree=butterflyTree
coverage=.55
leng=30

getConservedRegions(mafFile,gffFile,scaffoldName,referenceName,newickTree,30,.55)

getConservedRegions <- function(mafFile,gffFile,scaffoldName, referenceName, newickTree, leng, coverage){
  # read alignment
  align <- read.msa(mafFile)
  # read gene annotations from a gff file
  feats <- read.feat(gffFile)
  #The seqNames for the gff file are all the scaffold names (eg Hmel200001). They need to be HmelRef, but we don't want to mix up the scaffold positions. 
  #Therefore, must do one scaffold at a time.
  scaffoldFeats <- subset(feats, seqname==scaffoldName)
  scaffoldFeats$seqname <- referenceName
  
  #add introns
  fullScaffoldFeats <- add.introns.feat(scaffoldFeats)
  fullScaffoldFeats <- fullScaffoldFeats[fullScaffoldFeats$feature != "exon",]
  
  #make feature that represents whole chromosome
  wholeChrom <- feat(seq=referenceName, src=".", feature="all",
                     start=align$offset,
                     end=align$offset+ncol.msa(align, referenceName))
  
  # annotate intergenic regions
  intergenicFeats <- inverse.feat(fullScaffoldFeats, region.bounds=wholeChrom)
  intergenicFeats$feature <- "intergenic"
  fullScaffoldFeats <- rbind.feat(fullScaffoldFeats, intergenicFeats)
  
  #predict neutral model
  align4d <- get4d.msa(align, fullScaffoldFeats)
  neutralMod <- phyloFit(align4d, tree=newickTree, subst.mod="REV")
  
  #predict conserved model
  
  #predict conserved elements with phastCons
  pc <- phastCons(align, neutralMod, expected.length=leng,target.coverage=coverage,viterbi=TRUE)
  consElements <- pc$most.conserved
  output <- consElements
  output$seqname=scaffoldName
  write.feat(output,paste0(scaffoldName,"conservedElements.gff"))
  
  #number of conserved bases
  #coverage.feat(consElements)
  # this shows the fraction of bases covered by conserved elements
  #coverage.feat(consElements)/coverage.feat(wholeChrom)
  
  #For comparison, we will produce an alternative set of conservation scores using phyloP.
  pp <- phyloP(neutralMod, align, method="LRT", mode="CONACC")
  # the returned object is a data frame giving statistics for every base
  # in the alignment
  
  #Let us now plot the gene annotations, conserved elements, and conservation scores for a genomic segment of
  #interest. We will make use of functions in RPHAST that allow “tracks” to be defined and then plotted in a
  #browser-like display.
  par(mfrow=c(1,1))
  codingFeats <- fullScaffoldFeats[fullScaffoldFeats$feature=="CDS",]
  geneTrack <- as.track.feat(codingFeats, "genes", is.gene=TRUE)
  consElTrack <- as.track.feat(consElements, "phastCons most conserved", col="red")
  phastConsScoreTrack <- as.track.wig(wig=pc$post.prob.wig,
                                      name="phastCons post prob", col="red", ylim=c(0, 1))
  phyloPTrack <- as.track.wig(coord=pp$coord, score=pp$score, name="phyloP score",
                              col="blue", smooth=TRUE, horiz.line=0)
  plot.track(list(geneTrack, consElTrack, phastConsScoreTrack, phyloPTrack),
             xlim=c(0, 100000), cex.labels=1.25, cex.axis=1.25, cex.lab=1.5,main=paste(scaffoldName,"Conserved Elements"))
}
  
# define tree using Newick string
butterflyTree <- "(((HmelRef,HmelDisco),(Hcyd,Htim)),Hnum);"
lepTree <- "(Papilio_glaucus_v1x1 (Lerema_accius_v1x1 (Danaus_plexippus_v3 (Bicyclus_anynana_v1 (Melitaea_cinxia_v1 Heliconius_melpomene_v2)))));"
drosTree <- '((droGri2:0.183954,droVir3:0.093575):0.000000,(droMoj3:0.110563,((((droBip:0.034265,droAna3:0.042476):0.121927,(droKik:0.097564,((droFic:0.109823,(((dm3:0.023047,(droSim1:0.015485,droSec1:0.015184):0.013850):0.016088,(droYak2:0.026909,droEre2:0.029818):0.008929):0.047596,(droEug:0.102473,(droBia:0.069103,droTak:0.060723):0.015855):0.005098):0.010453):0.008044,(droEle:0.062413,droRho:0.051516):0.015405):0.046129):0.018695):0.078585,(droPer1:0.007065,dp4:0.005900):0.185269):0.068212,droWil1:0.259408):0.097093):0.035250);'

getConservedRegions("data/subTree_18Genomes_Hmel201009.maf","data/Hmel2.gff","Hmel201009","HmelRef",butterflyTree, 6,.50)




# extract genic regions, translate, and re-align --------------------------

extractFeatures <- function(mafFile,gffFile, scaffoldName, referenceName){
  align <- read.msa(mafFile)
  # read gene annotations from a gff file
  feats <- read.feat(gffFile)
  scaffoldFeats <- subset(feats, seqname==scaffoldName)
  scaffoldFeats$seqname <- referenceName
  scaffoldFeats <- scaffoldFeats[scaffoldFeats$feature == "exon",]
  return(extract.feature.msa(x = align,features = scaffoldFeats))
}

genes <- extractFeatures("data/subTree_18Genomes_Hmel201009.maf","data/Hmel2.gff","Hmel201009","HmelRef")

translated <- translate.msa(genes)
