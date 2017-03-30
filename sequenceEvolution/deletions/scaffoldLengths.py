#!/usr/bin/env python
from Bio import SeqIO


genomeFiles=["/n/holylfs/INTERNAL_REPOS/PCACTUS/edelman/genomes/1kbFilter/Hmel2.fa_1kb.fa","/n/holylfs/INTERNAL_REPOS/PCACTUS/edelman/genomes/1kbFilter/Bombyx_mori_ASM15162v1_-_scaffolds.fa_1kb.fa",\
"/n/holylfs/INTERNAL_REPOS/PCACTUS/edelman/genomes/1kbFilter/Bicyclus_anynana_v1.2_-_scaffolds.fa_1kb.fa"]

def scafSizeHash(scaf,numSlots):
    ords=[ord(i) for i in scaf]
    hashKey=sum(ords)%numSlots
    return hashKey

def makeDict(files,hash,numSlots):
    newDict={}
    for i in range(numSlots):
        newDict[i]={}
    for file in files:
        f=SeqIO.parse(open(file,"r"),"fasta")
        for record in f:
            scafName=record.name.split()[0]
            hashKey=scafSizeHash(scafName,numSlots)
            newDict[hashKey][scafName]=len(record)
    return newDict

scaffoldSizeDict=makeDict(genomeFiles,scafSizeHash,75)

def getScaffoldSize(scaf,scafDict):
    hashKey=scafSizeHash(scaf,75)
    size=scafDict[hashKey][scaf]
    return size
