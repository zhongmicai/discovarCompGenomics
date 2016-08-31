#!/usr/bin/env python

import argparse
import os

parser = argparse.ArgumentParser(description='Find intervals present in outgroups but absent \
in ingroups.')
parser.add_argument('--ingroups', type=str, nargs='+',
                    help='bed files of ingroups')
parser.add_argument('--outgroups', type=str, nargs='+',
                    help='bed files of outgroups')
parser.add_argument('--feat', type=str,
                    help='feature (gff or bed) file to compare deletions to')
parser.add_argument('--overlap', default='1.0',type=str,
                    help='amount of overlap necessary to report deleted feature')
parser.add_argument('-o', type=str,
                    help='prefix of output files')
args = parser.parse_args()

allSpecs = args.outgroups+args.ingroups

#Run bedtools multiintersct on all the input files
minterCommand="bedtools multiinter -i"
for spec in allSpecs:
    minterCommand += ''' %s ''' % (spec)
minterCommand+="> intersect_tmp.out"

print minterCommand
os.system(minterCommand)

#Parse the output to get only those intervals that are present in the outgroups
#but not in the ingroups
numOuts=len(args.outgroups)
correctList=''
for num,s in enumerate(args.outgroups):
    correctList+= str(num+1)+","
correctList=correctList[:-1]

parsingCommand='''awk '$4==%s' intersect_tmp.out | awk '$5=="%s"' > %s''' % \
(str(numOuts),correctList,args.o+"_allDeletions.bed")

print parsingCommand
os.system(parsingCommand)

#If provided, compare the deleted regions with a set of features (genes, exons, etc)
if args.feat:
    interCommand='''bedtools intersect -f %s -a %s -b %s > %s ''' % \
    (args.overlap, args.feat, args.o+"_allDeletions.bed", args.o+"_deletedFeatures.bed" )
print interCommand
os.system(interCommand)

#clean up
os.system("rm -f intersect_tmp.out")
