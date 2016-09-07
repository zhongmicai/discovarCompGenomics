#!/usr/bin/env python

import sys
sys.path.insert(1,"/Users/edelman/Documents/programming/alignio-maf")
sys.path.insert(1,"/Users/edelman/Documents/programming/bcbb/gff")
from Bio import AlignIO
from Bio.AlignIO import MafIO
from BCBio import GFF

#mafFile=sys.argv[1]
#gffFile=sys.argv[2]
#sequence=sys.argv[3]

gff = "data/Hmel2.gff"
maf="data/subTree_18Genomes_Hmel201001.maf"
gff_handle = open(gff)
sequence="HmelRef.Hmel201001"

mafIndex=MafIO.MafIndex("trial.mafIndex",maf,sequence)

limit_info = dict(
        gff_id = [sequence.split(".")[1]],
        gff_type = ["exon"])


mafs=[]
for rec in GFF.parse(gff_handle, limit_info=limit_info):
    for i in range(len(rec.features)):
        search=mafIndex.search([rec.features[i].location.start],[rec.features[i].location.end])
        for align in sea
gff_handle.close()

m=mafs[0]

for i,m in enumerate(mafs):
    theseMafs=[align for align in m]
    AlignIO.write(theseMafs,open("Hmel201001_exonMafs"+str(i)+".maf","w"),"maf")

def getGenes(maf,gff):
    
    feats=GFF.parse()
