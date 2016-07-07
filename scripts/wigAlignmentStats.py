#!/usr/bin/env python

#extract statistics on the alignment depth in wig format.

import sys
import numpy

wigFile=sys.argv[1]
baseFile=wigFile[:-4]
coverageOut=baseFile+".tsv"

def countWig(inputWig):
    f=open(inputWig, "r")
    coverageDict={}
    for line in f:
        if "fixed" not in line:
            try:
                coverageDict[line]+=1
            except KeyError:
                coverageDict[line]=1
    out=open(coverageOut, "w")
    for key in coverageDict.keys():
        out.write('''%s\t%s\n''' % (key.strip("\n"), coverageDict[key]))
    out.close()
    #return coverageDict


countWig(wigFile)
