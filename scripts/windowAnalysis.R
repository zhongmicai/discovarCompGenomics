source("https://bioconductor.org/biocLite.R")
#install.packages("devtools")
library(devtools)
#install_github("bernatgel/karyoploteR")
#biocLite("karyoploteR")
library("karyoploteR")
library(ggplot2)

ABBABABA="introgression/melsil-erahim-sardem/ABBABABAsummary.csv"
TREES="introgression/melsil-erahim-sardem/filteredTrees.csv"

getChrom <- function (x) {
  chrom <- strsplit(x,"_")[[1]][1]
  return(chrom)
}

getEnd <- function (x) {
  end <- strsplit(x,"_")[[1]][2]
  return(end)
}

Hmel2pt5 <- read.delim("introgression/orderChroms/Hmel2.transitions.tsv") 
Hmel2pt5$Chromosome <- paste0("chr",Hmel2pt5$Chromosome)
Hmel2pt5Genome <- read.delim("introgression/orderChroms/Hmel2.transitions_fullChroms.bed", col.names=c("seqnames","start","end","chromosome"))
Hmel2pt5Genome$chr <-  unlist(strsplit(as.character(Hmel2pt5Genome$chromosome),"chr"))[c(FALSE,TRUE)]
Hmel2pt5Cyto <- data.frame(seqnames=Hmel2pt5$Chromosome, Hmel2pt5$ChromStart,Hmel2pt5$ChromEnd, strand=Hmel2pt5$Orientation, label=Hmel2pt5$Hmel2Scaffold)

Hmel2pt5Genome <- Hmel2pt5Genome[order(as.numeric(Hmel2pt5Genome[,5])),]

Hmel2pt5GenomeGrange <- makeGRangesFromDataFrame(Hmel2pt5Genome, seqnames.field = "seqnames",keep.extra.columns = TRUE)
Hmel2pt5CytoGrange <- makeGRangesFromDataFrame(Hmel2pt5Cyto)

trees <- read.csv(TREES, header=FALSE, col.names = c("region","tree","overallTree"))
trees$region <- as.character(trees$region)

Dstat <- read.csv(ABBABABA, header=FALSE,col.names = c("segment","alnLength","snps","biallelic","ABBA","BABA"))
Dstat$segment <- as.character(Dstat$segment)
Dstat$chrom <- as.character(unlist(lapply(Dstat$segment,getChrom)))
Dstat$end <- as.numeric(unlist(lapply(Dstat$segment,getEnd)))*25000+25000
Dstat$start <- Dstat$end-49999
Dstat$D <- (Dstat$ABBA-Dstat$BABA)/(Dstat$ABBA+Dstat$BABA)

allData <- cbind(Dstat,trees)

DstatLarge <- subset(allData,alnLength>20000 & alnLength<50000 & chrom!="chr0" & ABBA+BABA>20)

grangeInfo <- data.frame(seqNames=DstatLarge$chrom, start=DstatLarge$start, end=DstatLarge$end, value=DstatLarge$D, tree=DstatLarge$tree, ABBA=DstatLarge$ABBA, BABA=DstatLarge$BABA)
ABGranges <- makeGRangesFromDataFrame(df=grangeInfo,keep.extra.columns = TRUE)


alignPlot <- ggplot(data=DstatLarge)+
  geom_histogram(aes(x=alnLength),bins=100)
alignPlot

snpPlot <- ggplot(data=DstatLarge)+
  geom_histogram(aes(x=snps),bins=100)
snpPlot

biallelicSnpPlot <- ggplot(data=DstatLarge)+
  geom_histogram(aes(x=biallelic),bins=100)
biallelicSnpPlot

ABBAplusBABA <- ggplot(data=DstatLarge)+
  geom_histogram(aes(x=ABBA+BABA),bins=100)
ABBAplusBABA

ABBAandBABA <- ggplot(data=DstatLarge)+
  geom_histogram(aes(x=ABBA), binwidth=2, fill="blue")+
  geom_histogram(aes(x=BABA), binwidth=2, fill="orange", alpha=0.5)+
  labs(x="ABBA vs BABA")
ABBAandBABA

hist(Dstat$BABA, breaks=seq(min(Dstat$BABA), max(Dstat$BABA),2))

getMeanSE <- function(seqName){
  vec <- subset(grangeInfo, seqNames==seqName)$value
  return(mean_se(vec))
}

x <- as.character(unique(grangeInfo$seqNames))
summaryStats <- data.frame(seqNames=character(),start=numeric(),end=numeric(),Dmin=numeric(), Dmax=numeric(), ABBAmin=numeric(), 
ABBAmax=numeric(), BABAmin=numeric(), BABAmax=numeric())
for (chr in x){
  length=Hmel2pt5Genome[which(Hmel2pt5Genome$chromosome==chr),3]
  DmeanSE=mean_se(subset(grangeInfo, seqNames==chr)$value)
  AmeanSE=mean_se(subset(grangeInfo, seqNames==chr)$ABBA)
  BmeanSE=mean_se(subset(grangeInfo, seqNames==chr)$BABA)
  newVals <- data.frame(seqNames=chr, start=(length/2)-1000000, end=(length/2)+1000000, Dmin=DmeanSE$ymin, Dmax=DmeanSE$ymax,
                        ABBAmin=AmeanSE$ymin, ABBAmax=AmeanSE$ymax,BABAmin=BmeanSE$ymin, BABAmax=BmeanSE$ymax)
  summaryStats <- rbind(summaryStats,newVals)
}
summaryGranges <-  makeGRangesFromDataFrame(summaryStats,keep.extra.columns = TRUE)


totalMeanSE <- mean_se(DstatLarge$D)

####### Plot data ###########

for(chr in as.character(unique(grangeInfo$seqNames))){
  
  tiff(filename = paste0("melsil_erahim_sardem_","all",".tiff"), width=2000,height=1000)
  tiff(filename = paste0("mel_timcyd_parnum_","all",".tiff"), width=2000,height=1000)
  #Set up the plot
  kp <- plotKaryotype(genome=Hmel2pt5GenomeGrange,plot.type=4)#, chromosome = c(chr))

  #Plot all data as scatterplot
  dataMin=round(min(ABGranges$value)-0.1,1)
  dataMax=round(max(ABGranges$value)+0.1,1)
  kpDataBackground(kp, data.panel = 1, r0=0, r1=0.35)
  kpAxis(kp, ymin=dataMin, ymax=dataMax, r0=0, r1=0.35, col="gray50", cex=1, numticks = length(seq(dataMin,dataMax,0.1)))
  kpAbline(kp, h=seq(dataMin,dataMax,0.1), col="grey50",ymin=dataMin, ymax=dataMax, r0=0, r1=0.35)
  kpAbline(kp, h=0, col="red",ymin=dataMin, ymax=dataMax, r0=0, r1=0.35, lwd=2)
  kpPoints(kp, data=ABGranges,cex=.5,r0=0, r1=0.35,ymin=dataMin, ymax=dataMax)

  #Plot mean+/- SE for each chromosome
  sumMin=round(min(summaryGranges$Dmin)-0.1,1)
  sumMax=round(max(summaryGranges$Dmax)+0.1,1)
  kpDataBackground(kp, data.panel = 1, r0=0.45, r1=0.6)
  kpAxis(kp, ymin=sumMin, ymax=sumMax, r0=0.45, r1=0.6, col="gray50", cex=1, numticks = length(seq(sumMin,sumMax,0.1)))
  kpAbline(kp, h=seq(sumMin,sumMax,0.1), col="gray50",ymin=sumMin, ymax=sumMax, r0=0.45, r1=0.6)
  kpAbline(kp, h=0, col="red",ymin=sumMin, ymax=sumMax, r0=0.45, r1=0.6,lwd=2)
  kpRect(kp, data=summaryGranges,y0=(summaryGranges$Dmin+summaryGranges$Dmax)/2,y1=summaryGranges$Dmax,r0=0.45, r1=0.6,ymin=sumMin, ymax=sumMax)
  kpRect(kp, data=summaryGranges,y1=(summaryGranges$Dmin+summaryGranges$Dmax)/2,y0=summaryGranges$Dmin,r0=0.45, r1=0.6,ymin=sumMin, ymax=sumMax)  
  kpAbline(kp, h=totalMeanSE$y,r0=0.45, r1=0.6,ymin=sumMin, ymax=sumMax, col="blue", lty="dashed")
  
  
  #Plot mean +/- SE for numABBA and numBABA for each chromosome
  ABmin=round(min(summaryGranges$ABBAmin, summaryGranges$BABAmin)-5,-1)
  ABmax=round(max(summaryGranges$ABBAmax, summaryGranges$BABAax)+5,-1)
  kpDataBackground(kp, data.panel = 1, r0=0.65, r1=0.95)
  kpAxis(kp, ymin=ABmin, ymax=ABmax, r0=0.65, r1=0.95, col="gray50",cex=1, numticks = length(seq(ABmin,ABmax,10)))
  kpAbline(kp, h=seq(ABmin,ABmax,10), col="gray50",ymin=ABmin, ymax=ABmax, r0=0.65, r1=0.95)
  kpRect(kp, data=summaryGranges,y0=(summaryGranges$ABBAmin+summaryGranges$ABBAmax)/2,y1=summaryGranges$ABBAmax,ymin=ABmin, ymax=ABmax, r0=0.65, r1=0.95, col="orange")
  kpRect(kp, data=summaryGranges,y1=(summaryGranges$ABBAmin+summaryGranges$ABBAmax)/2,y0=summaryGranges$ABBAmin,ymin=ABmin, ymax=ABmax, r0=0.65, r1=0.95, col="orange") 
  kpRect(kp, data=summaryGranges,y0=(summaryGranges$BABAmin+summaryGranges$BABAmax)/2,y1=summaryGranges$BABAmax,ymin=ABmin, ymax=ABmax, r0=0.65, r1=0.95, col="purple")
  kpRect(kp, data=summaryGranges,y1=(summaryGranges$BABAmin+summaryGranges$BABAmax)/2,y0=summaryGranges$BABAmin,ymin=ABmin, ymax=ABmax, r0=0.65, r1=0.95, col="purple") 
  
  
  #Plot trees
  kpDataBackground(kp, data.panel = 1, r0=.36, r1=.40, ymin=0, ymax=1, col="black")
  kpRect(kp, subset(ABGranges,tree=="Expected"),y0=0,y1=1,r0=.36, r1=.40, col="royalblue", border=NA)
  kpRect(kp, subset(ABGranges,tree=="Pardalinus/Numata+Melpomene"),y0=0,y1=1,r0=.36, r1=.40, col="yellow", border=NA)
  
  
  
  dev.off()

}
