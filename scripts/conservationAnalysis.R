#!/usr/bin/env Rscript

#ConservationAnalysis.R
require(rphast)


# find and graph most conserved elements ----------------------------------



getConservedRegions <- function(mafFile,gffFile,scaffoldName, referenceName, newickTree){
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
  pc <- phastCons(align, neutralMod,viterbi=TRUE, target.coverage=.25)# expected.length=12,target.coverage=0.525,viterbi=TRUE)
  consElements <- pc$most.conserved
  
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
  jpeg("trial.jpeg")
  plot.track(list(geneTrack, consElTrack, phastConsScoreTrack, phyloPTrack),
             xlim=c(0, 50000), cex.labels=1.25, cex.axis=1.25, cex.lab=1.5,main=paste(scaffoldName,"Conserved Elements"))
  dev.off()
}
  
# define tree using Newick string
butterflyTree <- "(((((HmelRef,HmelDisco),(Hcyd,Htim)),Hnum),Hera),Etal);"
subTree <- "(((HmelRef,HmelDisco),(Hcyd,Htim)),Hnum);"

getConservedRegions("results/finalAssemblies_expandedSubTree/finalAssemblies_expandedSubTree_Hmel201009.maf","data/Hmel2.gff","Hmel201009","HmelRef",butterflyTree)
mafFile="data/subTree_18Genomes_Hmel201001.maf"
gffFile="data/Hmel2.gff"
scaffoldName="Hmel201001"
referenceName="HmelRef"

# extract genic regions, translate, and re-align --------------------------

extractFeatures <- function(mafFile,gffFile, scaffoldName, referenceName){
  align <- read.msa(mafFile)
  # read gene annotations from a gff file
  feats <- read.feat(gffFile)
  scaffoldFeats <- subset(feats, seqname==scaffoldName)
  scaffoldFeats$seqname <- referenceName
  scaffoldFeats <- scaffoldFeats[scaffoldFeats$feature == "exon",]
  return(split.by.feature.msa(x=align,f=scaffoldFeats))
  #return(extract.feature.msa(x = align,features = scaffoldFeats))
}

genes <- extractFeatures("data/subTree_18Genomes_Hmel201001.maf","data/Hmel2.gff","Hmel201001","HmelRef")
mafFile <- "data/subTree_18Genomes_Hmel201001.maf"
gffFile <- "data/Hmel2.gff"
scaffoldName <- "Hmel201001"
referenceName <- "HmelRef"
