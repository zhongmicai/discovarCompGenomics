#!/bin/env python

import sys
from Bio import SeqIO

genome=sys.argv[1]
outFile=sys.argv[2]

g=SeqIo.parse(open(genome,"r"),"fasta")
out=open(outFile,"w")
startPos=0
for record in g:
    out.write('''%s\t%i\n''' % (record.id, len(record+startPos))
    startPos+=len(record)
out.close()
