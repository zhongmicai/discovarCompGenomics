#!/usr/bin/env python

##program/scripts to go from hal file to protein alignements of genes with deleted exons

import argparse
import os
import copy
from Bio import SeqIO
from gffClass import *
from scaffoldLengths import *

parser = argparse.ArgumentParser(description='Find and Validate Genes of Interest')

subparsers=parser.add_subparsers(help='sub-command help')

#find Deletions parser
parser_findDeletions=subparsers.add_parser('findDeletions', help='findDeletions help')
parser_findDeletions.add_argument('alignment', type=str,
                                    help='path to alignment hal file')
parser_findDeletions.add_argument('refGenome', type=str,
                                    help='name of reference genome')
parser_findDeletions.add_argument('ingroupList',type=str,
                                    help='comma-delimeted list of ingroup genome names')
parser_findDeletions.add_argument('outputFile', type=str,
                                    help='output file name')
parser_findDeletions.add_argument('-i','--intersectFile', nargs='+', type=str,
                                    help='bed or gff file(s) to intersect with can list \
                                    multiple files separated by spaces')
parser_findDeletions.add_argument('-wo', "--writeOverlap", action='store_true',
                                    help='If using -i. Option from bedtools; output entries from both \
                                    bed files (intersect and deletions) as well as overlap length. \
                                    Default is to only output entry from intersectFile. This is a flag.')
parser_findDeletions.add_argument('-f', '--percentOverlap', type=float,
                                    help='If using -i. Option from bedtools; require -f fraction overlap of the entry \
                                    in the intersectFile (input as eg .1 for 10 pct of the feature covered). \
                                    Default is to output all intersects.')
parser_findDeletions.set_defaults(which='findDeletions')


### Extract protein alignment from bed ###

parser_extractProteinAlignments=subparsers.add_parser('protAlign', help='protAlign help')
parser_extractProteinAlignments.add_argument('regionsOfInterest', type=str,
                                                help='gff file with coordinates and names of *genic* regions of interest')
parser_extractProteinAlignments.add_argument('refGFF', type=str,
                                                   help='full GFF file from reference species')
parser_extractProteinAlignments.add_argument('targetGFF', type=str,
                                                   help='comma separated list of GFF file(s) of genes from target species. ')
parser_extractProteinAlignments.add_argument('refDB', type=str,
                                                   help='BLAST database of proteins from reference species')
parser_extractProteinAlignments.add_argument('targetDB', type=str,
                                                   help='comma separated list of BLAST database(s) of proteins from target species. ')
parser_extractProteinAlignments.add_argument('speciesOfInterest', type=str,
                                                   help='comma separated list of species in format that they appear in progressive cactus alignment. Reference genome first; \
                           others in same order as listed in targetDB')
parser_extractProteinAlignments.add_argument('halFile', type=str,
                                                help='hal file with pregressive cactus alignment')
parser_extractProteinAlignments.set_defaults(which='protAlign')

### parse arguments ###

arguments=parser.parse_args()

##### The protAlign Command #####
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

def ucscToBed(scaf,start,end):
    scafLength=getScaffoldSize(scaf,scaffoldSizeDict)
    newStart=scafLength-int(end)+1
    newEnd=scafLength-int(start)+1
    return (newStart,newEnd)

def getHitCoordinates(mafFile):
    m=open(mafFile,"r")
    allCoords={}
    for line in m:
      if line.startswith("s"):
        atts=line.split()
        chrom=atts[1]
        start=atts[2]
        length=atts[3]
        strand=atts[4]
        try:
          allCoords[chrom][1].append(int(start))
          allCoords[chrom][1].append(int(start)+int(length))
        except KeyError:
          allCoords[chrom]=[strand,[int(start),int(start)+int(length)]]
    out={}
    for key in allCoords.keys():
        endName=key.index(".")
        scaf=key[endName+1:]
        start=min(allCoords[key][1])
        end=max(allCoords[key][1])
        strand=allCoords[key][0]
        if strand =='-':
            newCoords=ucscToBed(scaf,start,end)
            start=newCoords[0]
            end=newCoords[1]
        out[key]=[scaf,start,end,strand]
    print(mafFile,out)
    return out

def getNucAlignment(bedCoordinates, refGenome, halFile, speciesList):
    '''take bed coordinates and a species list; return a maf format alignment of those species.'''
    mafFile=bedCoordinates[:-3]+"maf"
    cmd='''hal2maf --noAncestors --maxRefGap 1000 --refGenome %s --targetGenomes %s --refTargets %s %s %s''' % (refGenome,speciesList,bedCoordinates,halFile,mafFile)
    print(cmd)
    os.system(cmd)

def makePositiveGff(gffFile):
    o=open(gffFile,"r")
    p=open(gffFile[:-5]+"_positive.gff","w")
    for line in o:
        new=line.replace("-","+")
        p.write(new)
    o.close()
    p.close()
    return gffFile[:-5]+"_positive.gff"

def writeBed(name,attList):
    bed=open(name+".bed","a+")
    scaf=attList[0]
    start=attList[1]
    end=attList[2]
    strand=attList[3]
    bed.write('''%s\t%s\t%s\t%s\t%s\t%s\n''' % (scaf,start,end,name,"0",strand))
    bed.close()

def getHitIntersects(hitDict,data,refGenome,gene):
    specs=[]
    for hit in hitDict.keys():
        species=hit.split(".")[0]
        if (species not in specs) and (species != refGenome):
            specs.append(species)
        if species != refGenome:
            #print(hitDict[hit])
            name=gene+"_"+species
            writeBed(name,hitDict[hit])
            intersectCMD='''intersectBed -wa -f 2E-9 -a %s -b %s >> %s ''' % (data[species][2],name+".bed",name+"_intersect.bed")
            print(intersectCMD)
            os.system(intersectCMD)
    return specs

def getProtIDs(gffEntry,gffGenes):
    protIds=[]
    geneID=gffEntry.getGeneId()
    for gene in gffGenes:
        if gene.getGeneId() == geneID:
            protIds+=gene.getProteinIDs()
    #print(protIds)
    return protIds

def getGeneEntry(gffEntry,gffGenes):
    for gene in gffGenes:
        try:
            if gene.getGeneId()==gffEntry.getGeneId():
                return gene.getGene()[0]
        except IndexError:
            #print(gene.getAllEntries()[0].getGeneId() + "has no gene entry?")
            pass

def getProteinSequences(protIDList,protFasta,outputHandle,species):
    prots=SeqIO.parse(open(protFasta, "r"),"fasta")
    for prot in prots:
        if prot.name in protIDList:
            hit=prot
            hit.name=prot.name+"|"+species
            hit.id=prot.id+"|"+species
            SeqIO.write([hit,],outputHandle,"fasta")


def runProtAlign(queryGFF,refGFF,refDB,targetDB, targetGFF,speciesOfInterest,halFile):
    '''take a gff file, locations to reference and target protein databases, and a list of species of interest
    return protein alignments of each protein in the gff between all species in list'''
    #First, set up the environment
    #os.system('source /n/sw/progressiveCactus-latest/progressiveCactus/environment')
    #os.system('module load bedtools')
    #os.system('module load mafft')
    #os.system('export PYTHONPATH=/n/home09/nedelman/natesbiopython/lib/python3.4/site-packages/:$PYTHONPATH')

    #Read in data and organize it into a dictionary
    species=speciesOfInterest.split(",")
    refGenome=species[0]
    dbs=[refDB,]+(targetDB.split(","))
    gffs=[refGFF,]+(targetGFF.split(","))
    refGffEntries=makeGffEntries(refGFF)
    refGffGenes=makeFullGeneEntries(refGffEntries)
    bedtoolsGffs=[]
    #fullGeneData=[]
    for gff in gffs:
        bedtoolsGffs.append(makePositiveGff(gff))
        # gffEntries=makeGffEntries(gff)
        # fullGeneData.append(makeFullGeneEntries(gffEntries))
    data={}
    for s in range(len(species)):
        data[species[s]]=[dbs[s],gffs[s],bedtoolsGffs[s]]

    #Go through the unique genes of interest, get their alignments, and find the corresponding protein(s) from each species of interest
    regionsOfInterest=makeGffEntries(queryGFF)
    geneNames=[]
    genesOfInterest=[]
    for region in regionsOfInterest:
        if region.getGeneId() not in geneNames:
            geneNames.append(region.getGeneId())
            genesOfInterest.append(getGeneEntry(region,refGffGenes))
    for gene in genesOfInterest:
        geneFasta=open(gene.getGeneId()+"_proteins.fa","a+")
        gene.writeBed()
        getNucAlignment(gene.getGeneId()+".bed", refGenome, halFile, speciesOfInterest)
        hitDict=getHitCoordinates(gene.getGeneId()+".maf")
        hitSpecies=getHitIntersects(hitDict,data,refGenome,gene.getGeneId())
        os.system('''msa_view %s.maf --in-format MAF --missing-as-indels --seqs %s > %s''' % (gene.getGeneId(),speciesOfInterest,gene.getGeneId()+"_nucAlignment.fa"))
        for species in hitSpecies:
            hitFile=gene.getGeneId()+"_"+species+"_intersect.bed"
            hitEntries=makeGffEntries(hitFile)
            hitGenes=makeFullGeneEntries(hitEntries)
            specProtIDs=[]
            for hit in hitGenes:
                specProtIDs+=(hit.getProteinIDs())
                print hit.getGeneId(), hit.getProteinIDs()
            print(species, specProtIDs)
            getProteinSequences(specProtIDs,data[species][0],geneFasta,species)

        refProtIDs=getProtIDs(gene,refGffGenes)
        print("ref",refProtIDs)
        refProtSequence=getProteinSequences(refProtIDs,refDB,geneFasta,refGenome+"_ref")
        geneFasta.close()
        #os.system('''clustalw2 -infile=%s -type=protein -outfile=%s -PIM -outorder=alignment''' % (geneFasta, gene.getGeneId()+"_proteinAlign.out"))
        #os.system()#MAFFT alignment




if arguments.which=='protAlign':
    runProtAlign(arguments.regionsOfInterest,arguments.refGFF,arguments.refDB,arguments.targetDB,arguments.targetGFF,arguments.speciesOfInterest,arguments.halFile)

##### The findDeletions Command ######
def getDeletions(halFile,refGenome,ingroupGenomes,outputFile):
    print ("finding regions exclusively in:", ingroupGenomes)
    print ("using ", refGenome, "as the reference...")
    cmd='''sbatch findDeletions.slurm %s %s %s %s ''' % (halFile,refGenome,ingroupGenomes, outputFile)
    print (cmd)
    os.system(cmd)

def getIntersections(deletionFile,intersectFile, options):
    print ('''finding intersections betweeen %s and %s ...''' % (deletionFile,intersectFile))
    baseDeletions=os.path.basename(deletionFile)[:-3]
    baseIntersects=os.path.basename(intersectFile)[:-3]
    output='''%s_%s_intersects.bed''' % (deletionFile,intersectFile)
    cmd='''intersectBed %s -a %s -b %s > %s''' % (options, intersectFile,deletionFile,output)
    print (cmd)
    os.system(cmd)

if arguments.which=='findDeletions':
    getDeletions(arguments.alignment,arguments.refGenome, arguments.ingroupList,arguments.outputFile)
    if arguments.i:
        optionList=["-wa",""]
        if arguments.wo:
            optionList[0]="-wo"
        if arguments.f:
            optionList[1]='''-f %s''' % (arguments.f)
        options='''%s %s''' % (optionList[0],optionList[1])
        for intersect in arguments.i:
            getIntersections(arguments.outputFile, i,options)
