#!/usr/bin/env python
# coding: utf-8

"""
- NanoSum ver0.1.1
- Written by Zilin Ren
- Function: make a summary for fastq file and bam file
"""

import argparse

from myLib.myFq import myFqVis
from myLib.myBam import myBamVis
from myLib.myBed import myBedVis



def run_myFqVis(args):
    sampleID, fqPath = args.s, args.fq
    myFqVis(sampleID, fqPath)
    return 

def run_myBamVis(args):
    sampleID, bamPath, samtools = args.s, args.b, args.st
    myBamVis(sampleID, bamPath, samtools)
    return 

def run_myBedVis(args):
    sampleID, bedPath, bamPath, samtools =  args.s, args.bd, args.b, args.st
    # myBamVis(sampleID, bamPath, samtools)
    myBedVis(sampleID, bedPath, bamPath, samtools)
    return 


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Visualization on FASTQ/BAM/BAM+BED', 
                                     usage='''python visNano.py <options> [<args>]

Available options are:    
    fq    Visualization on FASTQ
    bam   Visualization on BAM
    bed   Visualization on BAM+BED
    ''')
    subparsers = parser.add_subparsers(help='sub-command help')
    subparsers.required = True


    fq_parser = subparsers.add_parser('fq', help='Visualization on FASTQ')
    fq_parser.add_argument('--fq', help='path of fastq file (not *.gz file)', default = None)
    fq_parser.add_argument('--s', help='sample ID', default = None)
    fq_parser.set_defaults(func=run_myFqVis)


    bam_parser = subparsers.add_parser('bam', help='Visualization on BAM')
    bam_parser.add_argument('--s', help='sample ID', default = None)
    bam_parser.add_argument('--b', help='path of bam/sam file', default = None)
    bam_parser.add_argument('--st', help='path of SAMTOOLS', default = None)
    bam_parser.set_defaults(func=run_myBamVis)


    bed_parser = subparsers.add_parser('bed', help='Visualization on BAM+BED')
    bed_parser.add_argument('--s', help='sample ID', default = None)
    bed_parser.add_argument('--b', help='path of bam/sam file', default = None)
    bed_parser.add_argument('--bd', help='path of bed file', default = None)
    bed_parser.add_argument('--st', help='path of SAMTOOLS', default = None)
    bed_parser.set_defaults(func=run_myBedVis)

    ## run
    args = parser.parse_args()
    args.func(args)