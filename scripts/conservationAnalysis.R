#ConservationAnalysis.R
require(rphast)
# extract alignment and annotation files from RPHAST package

# read alignment
align <- read.msa("data/subTree_18Genomes_Hmel201009.maf")
align
# read gene annotations from a UCSC "genepred" file
feats <- read.feat("data/Hmel2.gff")
# define tree using Newick string
butterflyTree <- "(((HmelRef,HmelDisco),(Hcyd,Htim)),Hnum);"

#The seqNames for the gff file are all the scaffold names (eg Hmel200001). They need to be HmelRef, but we don't want to mix up the scaffold positions. 
#Therefore, must do one scaffold at a time.

Hmel201009Feats <- subset(feats, seqname=="Hmel201009")
Hmel201009Feats$seqname <- "HmelRef"

#add introns
Hmel201009Feats <- add.introns.feat(Hmel201009Feats)
Hmel201009Feats <- Hmel201009Feats[Hmel201009Feats$feature != "exon",]

#make feature that represents whole chromosome
wholeChrom <- feat(seq="HmelRef", src=".", feature="all",
                   start=align$offset,
                   end=align$offset+ncol.msa(align, "HmelRef"))

# annotate intergenic regions
intergenicFeats <- inverse.feat(Hmel201009Feats, region.bounds=wholeChrom)
intergenicFeats$feature <- "intergenic"
Hmel201009Feats <- rbind.feat(Hmel201009Feats, intergenicFeats)

#predict neutral model
align4d <- get4d.msa(align, Hmel201009Feats)
neutralMod <- phyloFit(align4d, tree=butterflyTree, subst.mod="REV")
neutralMod

#predict conserved elements with phastCons
pc <- phastCons(align, neutralMod, expected.length=30,target.coverage=0.25,viterbi=TRUE, rho=.1)
consElements <- pc$most.conserved

#number of conserved bases
coverage.feat(consElements)

# this shows the fraction of bases covered by conserved elements
coverage.feat(consElements)/coverage.feat(wholeChrom)

#For comparison, we will produce an alternative set of conservation scores using phyloP.
pp <- phyloP(neutralMod, align, method="LRT", mode="CONACC")
# the returned object is a data frame giving statistics for every base
# in the alignment
names(pp)

#Let us now plot the gene annotations, conserved elements, and conservation scores for a genomic segment of
#interest. We will make use of functions in RPHAST that allow “tracks” to be defined and then plotted in a
#browser-like display.
par(mfrow=c(1,1))
codingFeats <- Hmel201009Feats[Hmel201009Feats$feature=="CDS",]
geneTrack <- as.track.feat(codingFeats, "genes", is.gene=TRUE)
consElTrack <- as.track.feat(consElements, "phastCons most conserved", col="red")
phastConsScoreTrack <- as.track.wig(wig=pc$post.prob.wig,
                                    name="phastCons post prob", col="red", ylim=c(0, 1))
phyloPTrack <- as.track.wig(coord=pp$coord, score=pp$score, name="phyloP score",
                            col="blue", smooth=TRUE, horiz.line=0)
plot.track(list(geneTrack, consElTrack, phastConsScoreTrack, phyloPTrack),
           xlim=c(0, 500000), cex.labels=1.25, cex.axis=1.25, cex.lab=1.5)

