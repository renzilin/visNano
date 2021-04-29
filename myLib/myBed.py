# coding: utf-8

import os
import gzip
import pathlib

import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

def myBedVis(sampleID, bedPath, bamPath, samtools):
    def myFetchBam(bamfile, samtools, coordinate):
        file = os.popen('%s view %s %s' % (samtools, bamfile, coordinate)).readlines()
        # FLAG = []
        # RNAME = []
        # POS = []
        # SEQ = []
        SEQLen = []
        for line in file:
            lnLst  = line.strip().split()
            # FLAG.append( lnLst[1] )
            # RNAME.append( lnLst[2] )
            # POS.append( lnLst[3] )
            # SEQ.append( lnLst[9] )
            SEQLen.append( len(lnLst[9]) )

        return SEQLen

    def myFetchBed(bedPath):
        regionNameLst = []
        coordinateLst = []
        with open(bedPath, 'r') as file:
            for line in file:
                lnLst = line.strip().split(',')
                regionName = lnLst[0]
                chrom      = lnLst[1]
                start      = lnLst[2]
                end        = lnLst[3]

                regionNameLst += [regionName]
                coordinateLst += ['%s:%s-%s' % (chrom, start, end)]
        return regionNameLst, coordinateLst

    def myLociPlot(regionNameLst, lenLst, cntLst, sampleID):
        fig = plt.figure(figsize=(12, 10))
        plt.figure(1)

        ## Aligned Reads
        ax1 = plt.subplot(121)
        plt.boxplot(lenLst, vert=False, sym='+')
        plt.yticks([i for i in range(1, len(lenLst)+1)], regionNameLst)
        plt.xlabel('Regions')
        plt.ylabel('Read Length (bp)')


        ax2   = plt.subplot(122)
        rects = plt.barh(regionNameLst, cntLst, height=0.5)
        plt.ylim(-0.5, len(cntLst)-0.3)
        plt.xlim(0, max(cntLst)*1.15)
        plt.xlabel("Reads No.")
        ply.ylabel("Regions")

        for rect in rects:
            wd = rect.get_width()
            plt.text(wd, rect.get_y() + 0.25, '{:,}'.format(wd), va='center')

        
        outdir = pathlib.Path('LOG-visnano').mkdir(parents=True, exist_ok=True)
        plt.savefig("LOG-visnano/summary-on-bed-%s.jpg" % (sampleID), dpi=300)
        
        return



    regionNameLst, coordinateLst = myFetchBed(bedPath)

    cntLst  = []
    lenLst  = []

    for ind, regionName in tqdm(enumerate(regionNameLst), desc="Prcessing region: "):
        coordinate = coordinateLst[ind]
        SEQLen = myFetchBam(bamPath, samtools, coordinate)

        cntLst.append(len(SEQLen))
        lenLst.append(SEQLen)

    myLociPlot(regionNameLst, lenLst, cntLst, sampleID)

    return
    
    