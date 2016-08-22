#rPhast_tutorial.R
#from https://cran.r-project.org/web/packages/rphast/vignettes/vignette1.pdf
#additional info, especially for tuning expected.length and target.coverage: http://compgen.bscb.cornell.edu/phast/phastCons-HOWTO.html

install.packages("/Users/nbedelman/Downloads/rphast_1.6.4.tgz", repos=NULL)
library("rphast")

require(rphast)
# extract alignment and annotation files from RPHAST package
exampleArchive <- system.file("extdata","examples.zip",package="rphast")
unzip(exampleArchive, c("sol1.maf", "sol1.gp"))
# read alignment
align <- read.msa("sol1.maf")
# read gene annotations from a UCSC "genepred" file
feats <- read.feat("sol1.gp")
# define tree using Newick string
tomatoTree <- "((((tomato, potato), eggplant), pepper), petunia);"

#change alignment identifiers to common names
align$names <- c("tomato", "potato","pepper", "petunia", "eggplant")
#change id in features data to common name
feats$seqname <- "tomato"

#add introns
feats <- add.introns.feat(feats)
feats <- feats[feats$feature != "exon",]

# make a feature that represents the entire chromosome. We will ignore
# several thousand bases at the beginning of the reference genome for
# which no alignments are available by setting the "start" of the feature
# equal to the beginning of the aligned region
wholeChrom <- feat(seq="tomato", src=".", feature="all",
  start=align$offset,
  end=align$offset+ncol.msa(align, "tomato"))

# annotate intergenic regions
intergenicFeats <- inverse.feat(feats, region.bounds=wholeChrom)
intergenicFeats$feature <- "intergenic"
feats <- rbind.feat(feats, intergenicFeats)

#estimate neutral model from 4D sites using phyloFit
align4d <- get4d.msa(align, feats)
neutralMod <- phyloFit(align4d, tree=tomatoTree, subst.mod="REV")

#predict conserved elements with phastCons
pc <- phastCons(align, neutralMod, expected.length=6,target.coverage=0.125, viterbi=TRUE)
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

codingFeats <- feats[feats$feature=="CDS",]
geneTrack <- as.track.feat(codingFeats, "genes", is.gene=TRUE)
consElTrack <- as.track.feat(consElements, "phastCons most conserved", col="red")
phastConsScoreTrack <- as.track.wig(wig=pc$post.prob.wig,
                                    name="phastCons post prob", col="red", ylim=c(0, 1))
phyloPTrack <- as.track.wig(coord=pp$coord, score=pp$score, name="phyloP score",
                            col="blue", smooth=TRUE, horiz.line=0)
plot.track(list(geneTrack, consElTrack, phastConsScoreTrack, phyloPTrack),
           xlim=c(60000, 68000), cex.labels=1.25, cex.axis=1.25, cex.lab=1.5, main="Tomato tutorial")


# Now let us examine the predicted conserved elements in more detail. We will start by plotting their length
# distributions. We will plot distributions for all elements, and for the subsets of elements that primarily fall
# in coding or noncoding regions.
ce <- pc$most.conserved
plot(density.feat(ce), ylim=c(0, 0.018),
     main="Element Length by Type", xlab="Length",
     mgp=c(1.5,0.5,0),mar=c(2,2,2,2))
# obtain elements that overlap codingFeats by at least 50 percent
codingConsEle <- overlap.feat(ce, codingFeats, min.percent=0.5)
# obtain elements that overlap by less than 50 percent
noncodingConsEle <- overlap.feat(ce, codingFeats, min.percent=0.5,
                                 overlapping=FALSE)
lines(density.feat(codingConsEle), col="red")
lines(density.feat(noncodingConsEle), col="blue")
legend(c("All", "Coding", "Noncoding"), x="topright", inset=0.05,
       lty=1, col=c("black", "red", "blue"))

# Now let us examine the relationship between conserved elements and annotations of different types.
# First, we will plot the fold-enrichment of annotation types within conserved elements, as compared with the
# genomic segment as a whole. Second, we will plot the composition of conserved elements by annotation type,
# and compare it with the “background” composition for the entire region
par(mfrow=c(2, 2), cex.main=1.5, cex.lab=1.5, cex.axis=1.5, mar=c(5,5,4,2))
# look at fold-enrichment of each annotation type by conserved element
enrich <- enrichment.feat(ce, feats, wholeChrom)
col <- rainbow(nrow(enrich))
barplot(enrich$enrichment, col=col,
        main="Enrichment of\nConserved Elements",
        ylab="Fold Enrichment")
plot.new()
legend(x="center", legend=enrich$type, fill=col, cex=1.5)
# look at the composition of the conserved elements
comp <- composition.feat(ce, feats)
pie(comp$composition, col=rainbow(nrow(comp)), radius=1.0,
    main="Composition of\nConserved Elements", labels=NA)
# compare with background composition
comp <- composition.feat(wholeChrom, feats)
pie(comp$composition, col=rainbow(nrow(comp)), radius=1.0,
    main="Background\nComposition", labels=NA)


# Next, we will run phastCons again, but this time instead of fixing the the transition probabilities between
# the conserved and nonconserved states by specifying the “expected.length” and “target.coverage” parameters,
# we we will allow it to estimate these probabilities by maximum likelihood. It does this using an expectation
# maximization (EM) algorithm.
pcEM <- phastCons(align, neutralMod, viterbi=TRUE, estimate.transitions=TRUE)
