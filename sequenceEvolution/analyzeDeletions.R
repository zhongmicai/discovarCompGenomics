#!/usr/bin/env Rscript
#install.packages("ggthemes")
library(ggplot2)
library(ggthemes)
library(rphast)

helDeletions=read.csv("sequenceEvolution/helDeletions_toBmor_allGenomes/heliconiusDeletions_bombyxRef_allGenomes.bed",sep="\t", header = F, col.names = c("scaffold","start","end"))
helDeletions$length=helDeletions$end-helDeletions$start
helDeletionsNumBlocks=length(helDeletions$length)
helDeletionsTotal=sum(helDeletions$length)

notTinyDels <- subset(helDeletions, length>100)

ggplot(data=helDeletions, aes(x=length)) + 
  stat_density(aes(y=..count..)) + 
  scale_y_continuous(breaks=c(0,1,10,100,1000), trans="log1p") +
  ggtitle("Distribution of Deletion Size") +
  labs(y="count (log)")


DPSCF300401=read.msa("sequenceEvolution/DPOGS210663_orthologs.maf")
delRegion=sub.msa(DPSCF300401,seqs = c("Dple","Mcin","Avan","Etal" ,"Ldor","Hbes" ,"Hpar","Hnum" ,"Htim","Hcyd", "HmelRef","HmelDisco","Htel","Hhsa","HeraDisco","Hhim","HeraHhimHyb","Hsar" ,"Hdem" ,"HeraRef","Bany"),
        start.col = 186717,end.col=198609,refseq = "Dple")
delScaffold=delRegion=sub.msa(DPSCF300401,seqs = c("Dple","Mcin","Avan","Etal" ,"Ldor","Hbes" ,"Hpar","Hnum" ,"Htim","Hcyd", "HmelRef","HmelDisco","Htel","Hhsa","HeraDisco","Hhim","HeraHhimHyb","Hsar" ,"Hdem" ,"HeraRef","Bany"),
                              refseq = "Dple")
scaffoldFeats=read.feat("data/Dple_DPSCF300401_renamed.gff")

#add introns
fullScaffoldFeats <- add.introns.feat(scaffoldFeats)
fullScaffoldFeats <- fullScaffoldFeats[fullScaffoldFeats$feature != "exon",]

#make feature that represents whole chromosome
wholeChrom <- feat(seq="Dple", src=".", feature="all",
                   start=delScaffold$offset,
                   end=delScaffold$offset+ncol.msa(delScaffold, "Dple"))

# annotate intergenic regions
intergenicFeats <- inverse.feat(fullScaffoldFeats, region.bounds=wholeChrom)
intergenicFeats$feature <- "intergenic"
fullScaffoldFeats <- rbind.feat(fullScaffoldFeats, intergenicFeats)

#make neutral model

nymphalidTree <- "((((((((((Hpar,Hnum),Hbes),((Htim,Hcyd),(HmelRef,HmelDisco))),Ldor),(((((HeraRef,HeraDisco),(Hhim,HeraHhimHyb)),Hhsa),Htel),(Hsar,Hdem))),Etal),Avan),Mcin),Bany),Dple);"
align4d <- get4d.msa(delScaffold, fullScaffoldFeats)
neutralMod <- phyloFit(align4d, tree=smallTree, subst.mod="REV")

pc <- phastCons(delScaffold, neutralMod, viterbi=TRUE,expected.length=20,target.coverage=.01)
consElements <- pc$most.conserved
output <- consElements
output$seqname="DPSCF300401"
write.feat(output,paste0("DPSCF300401_conservedElements.gff"))

#For comparison, we will produce an alternative set of conservation scores using phyloP.
pp <- phyloP(neutralMod, delScaffold, method="LRT", mode="CON")

#Plot
par(mfrow=c(1,1))
codingFeats <- fullScaffoldFeats[fullScaffoldFeats$feature=="CDS",]
geneTrack <- as.track.feat(codingFeats, "genes", is.gene=TRUE)
consElTrack <- as.track.feat(consElements, "phastCons most conserved", col="red")
phastConsScoreTrack <- as.track.wig(wig=pc$post.prob.wig,
                                    name="phastCons post prob", col="red", ylim=c(0, 1))
phyloPTrack <- as.track.wig(coord=pp$coord, score=pp$score, name="phyloP score",
                            col="blue", smooth=TRUE, horiz.line=0)
#jpeg(paste0(scaffoldName,".jpeg"))
plot.track(list(geneTrack, consElTrack, phastConsScoreTrack, phyloPTrack),
           xlim=c(0, 100000), cex.labels=1.25, cex.axis=1.25, cex.lab=1.5,main=paste("DPSCF300401","Conserved Elements"))
#dev.off()


write.msa(delRegion,format="FASTA",file="deletedRegion_DPOGS210663.fa")
write.msa(delRegion,format="FASTA",file="deletionScaffold_DPSCF300401.fa")



