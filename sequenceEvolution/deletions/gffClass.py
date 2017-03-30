#!/usr/bin/env python

class gffEntry(object):
    def __init__(self,gffLine):
        self.entry=gffLine
        atts=gffLine.split("\t")
        self.chrom=atts[0]
        self.origin=atts[1]
        self.type=atts[2]
        self.start=atts[3]
        self.end=atts[4]
        self.strand=atts[6]
        self.fullDescription=atts[8]
    def getEntry(self):
        return self.entry
    def getChrom(self):
        return self.chrom
    def getOrigin(self):
        return self.origin
    def getType(self):
        return self.type
    def getStart(self):
        return self.start
    def getEnd(self):
        return self.end
    def getStrand(self):
        return self.strand
    def getFullDescription(self):
        return self.fullDescription
    # def getGenBank(self):
    #     desc=self.getFullDescription()
    #     desc=desc.replace(",",";")
    #     genBank=desc.split("Genbank:")[1].split(";")[0]
    #     return genBank
    def getGeneId(self):
        desc=self.getFullDescription()
        desc=desc.replace(",",";")
        desc=desc.replace(":",";")
        try:
            geneID=desc.split("commonID=")[1].split(";")[0].strip()
            return geneID
        except IndexError:
            print(desc)

    def writeBed(self):
        out=open(self.getGeneId()+".bed","w")
        out.write('''%s\t%s\t%s\t%s\t%s\t%s\n''' % (self.getChrom(),self.getStart(),self.getEnd(),self.getGeneId(),"0",self.getStrand()))
        out.close()

class fullGene(object):
    def __init__(self,gffEntries):
        self.allEntries=gffEntries
        self.allExons=[]
        self.allCDS=[]
        self.allmRNA=[]
        self.gene=[]
        for entry in gffEntries:
            if entry.getType()=="exon":
                self.allExons.append(entry)
            if entry.getType()=="CDS":
                self.allCDS.append(entry)
            if entry.getType()==("mRNA"):
                self.allmRNA.append(entry)
            if entry.getType()=="gene":
                self.gene.append(entry)
    def getAllEntries(self):
        return self.allEntries
    def getAllExons(self):
        return self.allExons
    def getAllCDS(self):
        return self.allCDS
    def getAllmRNA(self):
        return self.allmRNA
    def getGene(self):
        return self.gene
    def getGeneId(self):
        return self.getGene()[0].getGeneId()
    def getProteinIDs(self):
        ids=[]
        prots=self.getAllmRNA()
        for prot in prots:
            desc=prot.getFullDescription()
            desc=desc.replace(",",";")
            i=desc.split("translation_stable_id=")[1].split(";")[0].strip()
            ids.append(i)
        return ids


def makeGffEntries(gffFile):
    gff=open(gffFile,"r")
    out=[]
    for line in gff:
        out.append(gffEntry(line))
    return out

def makeFullGeneEntries(gffEntryList):
    geneDict={}
    for entry in gffEntryList:
        try:
            geneDict[entry.getGeneId()].append(entry)
        except KeyError:
            geneDict[entry.getGeneId()]=[entry,]
    fullGenes=[]
    for gene in geneDict.keys():
        fullGenes.append(fullGene(geneDict[gene]))
    return fullGenes
