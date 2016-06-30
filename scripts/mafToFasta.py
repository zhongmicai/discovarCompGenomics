#!/bin/env python

'''Takes a maf file and returns a fasta
Usage: mafToFasta.py <inputMaf> outputFasta '''

import sys

class mafAttributes(object):
    def __init__(self,mafLine):
        line=mafLine.split()
        try:
            self.type=line[0]
            try:
                self.label=line[1]
            except ValueError:
                self.label="NA"
            try:
                self.start=int(line[2])
            except ValueError:
                self.start=0
            try:
                self.end=self.start+int(line[3])
            except ValueError:
                self.end=0
            try:
                self.strand=line[4]
            except IndexError:
                self.strand="NA"
            try:
                self.seq=line[6].strip("\n")
            except IndexError:
                self.seq="NA"
        except IndexError:
            self.type="NA"
            self.label="NA"
            self.start="NA"
            self.end="NA"
            self.length="NA"
            self.strand="NA"
            self.seq="NA"

    def getType(self):
        return self.type
    def getLabel(self):
        return self.label
    def getStart(self):
        return self.start
    def getEnd(self):
        return self.end
    def getLength(self):
        return (self.getEnd()-self.getStart())
    def getStrand(self):
        return self.strand
    def getSeq(self):
        return self.seq

def mafToFasta(mafFile, outFile):
    '''  Takes a maf file
        outputs a fasta alignment file.'''
    output=open(outFile, "w")
    with open(mafFile, "r") as file:
        for line in file:
            if "#" in line:
                output.write(line)
            else:
                maf=mafAttributes(line)
                if maf.getType() =='s':
                    chrom=maf.getLabel()[:-3]
                    scaf=maf.getLabel()
                    chromStart=maf.getStart()
                    chromEnd=maf.getEnd()
                    length=maf.getLength()
                    seq=maf.getSeq()
                    output.write('''> %s_%s_%s_%s\n''' % (scaf, chromStart, chromEnd, length))
                    output.write(seq)
                    output.write("\n")
                elif line != "\n":
                    output.write("\n")
                else:
                    pass


mafInput=sys.argv[1]
fastaOutput=sys.argv[2]

mafToFasta(mafInput, fastaOutput)
