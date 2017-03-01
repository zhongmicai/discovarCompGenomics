#!/usr/bin/env python

##program/scripts to go from hal file to protein alignements of genes with deleted exons

import argparse
import os

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
parser_findDeletions.add_argument('-h', '--help', type=str,
                                    help='If using -i. Option from bedtools; require -f fraction overlap of the entry \
                                    in the intersectFile (input as eg .1 for 10 pct of the feature covered). \
                                    Default is to output all intersects.')
parser_findDeletions.set_defaults(which='findDeletions')


### Extract protein alignment from bed ###
parser_extractProteinAlignments=subparsers.add_parser('protAlign', help='protAlign help')
parser_extractProteinAlignments.add_argument('bedFile', type=str,
                                                help='bed file with ')


### parse arguments ###

arguments=parser.parse_args()

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
