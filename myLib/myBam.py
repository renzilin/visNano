#!/usr/bin/env python
# coding: utf-8

import os
import gzip
import pathlib

import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm


def myBamVis(sampleID, bamPath, samtools):
    def myFetchBam(bamfile, samtools):
        # print("## USE samtools to extract the reads covering the REGION;\n"+
        # "## the command: %s view -F 2304 %s %s" % (samtools, bamfile, coordinate))

        file = os.popen('%s view %s' % (samtools, bamfile)).readlines()

        FLAG = []
        RNAME = []
        POS = []
        # SEQ = []
        SEQLen = []
        for line in tqdm(file, desc="Reading SAM/BAM: "):
            lnLst  = line.strip().split()
            FLAG.append( lnLst[1] )
            RNAME.append( lnLst[2] )
            POS.append( lnLst[3] )
            # SEQ.append( lnLst[9] )
            SEQLen.append( len(lnLst[9]) )
            
        return FLAG, RNAME, POS, SEQLen

    def myBamPlot(FLAG, RNAME, POS, SEQLen, sampleID):
        fig = plt.figure(figsize=(12, 10))
        plt.figure(1)

        ## Aligned Reads
        ax1 = plt.subplot(121)

        validRNAME = ['chr%s' % i for i in range(1,23)] + ['chrX', 'chrY', '*']
        RNAMETab          = pd.Series(RNAME).value_counts().rename_axis("Chrom").reset_index(name='Reads').sort_values(by=['Chrom'])
        RNAMETab['Other'] = RNAMETab['Chrom'].apply(lambda x: 1 if x in validRNAME else 2)
        otherSum          = RNAMETab.loc[RNAMETab['Other'] == 2, 'Reads'].sum()

        count = []
        for i in validRNAME:
            count += RNAMETab['Reads'][RNAMETab['Chrom'] == i].tolist()

        validRNAME+= ['others']
        count     += [otherSum]

        rects = plt.barh(validRNAME[::-1], count[::-1], height=0.5)
        plt.xlim(0, max(count)*1.15)

        for rect in rects:
            wd = rect.get_width()
            plt.text(wd, rect.get_y() + 0.25, '{:,}'.format(wd), va='center')

        plt.xlabel('Reads No.')
        plt.ylabel('Chromosome')
        plt.title("Aligned Reads (except *, others): {:,}".format(sum(count[:24])))

        ax2 = plt.subplot(122)
        alignedLen = [SEQLen[i] for i in range(len(SEQLen)) if RNAME[i] in validRNAME[:-2]]
        plt.hist(alignedLen, bins=1000, orientation='horizontal')
        plt.ylim(0, 2000)
        plt.ylabel('Length (bp)')
        plt.xlabel('Reads No.')
        plt.title("Total reads: {:,}".format(len(alignedLen)))
        plt.suptitle(sampleID)
        outdir = pathlib.Path('LOG-visnano').mkdir(parents=True, exist_ok=True)
        plt.savefig("LOG-visnano/summary-on-sam-%s.jpg" % (sampleID), dpi=300)
        # plt.show()
        return 

    FLAG, RNAME, POS, SEQLen = myFetchBam(bamPath, samtools)
    myBamPlot(FLAG, RNAME, POS, SEQLen, sampleID)

    return 