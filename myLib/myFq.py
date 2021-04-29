#!/usr/bin/env python
# coding: utf-8

import os
import gzip
import pathlib

import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

def myFqVis(sampleID, fastqPath):
    def myFqInfo(fastqPath):
        def QClass(Qlst, L):
            Larr = np.array(L)
            for n, i in enumerate(range(5, 81, 5)):
                Qlst[n] += len(Larr[(Larr > i - 5)&(Larr <= i)])
            return Qlst
        
        readLength = []
        readScore  = [] 
        baseScore  = [0 for i in range(16)]
        if '*.gz' not in fastqPath:
            with open(fastqPath, 'r') as f:
                for lnCnt, ln in tqdm(enumerate(f), desc="Processing FASTQ: "):
                    if lnCnt % 4 == 1:
                        readLength.append( len(ln.strip()) )
                    if lnCnt % 4 == 3:
                        L = [ord(i) for i in ln.strip()]
                        readScore.append(sum(L)/len(L))
                        baseScore = QClass(baseScore, L)
        
        if '*.gz' in fastqPath:
            with gzip.open(fastqPath, 'r') as f:
                for lnCnt, ln in tqdm(enumerate(f), desc="Processing FASTQ: "):
                    if lnCnt % 4 == 1:
                        readLength.append( len(ln.strip()) )
                    if lnCnt % 4 == 3:
                        L = [ord(i) for i in ln.strip()]
                        readScore.append(sum(L)/len(L))
                        baseScore = QClass(baseScore, L)
        
        return baseScore, readLength, readScore

    def myFqPlot(baseScore, readLen, readScore, sampleID):
        fig = plt.figure(figsize=(12, 10))
        plt.figure(1)

        ## base score
        ax1 = plt.subplot(221)
        plt.bar([i for i in range(16)], baseScore)
        plt.xticks([i + 0.5 for i in range(1, 15, 2)], ['Q%s' % i for i in range(10, 80, 10)])
        plt.xlabel('Base Quality')
        plt.ylabel('Bases No.')
        plt.title("Total bases: {:,}".format(sum(baseScore)))

        ax2 = plt.subplot(222)
        plt.hist(readLen, bins=1000)
        plt.xlim(0, 2000)
        plt.xlabel('Length (bp)')
        plt.ylabel('Reads No.')
        plt.title("Total reads: {:,}".format(len(readLen)))

        ax3 = plt.subplot(223)
        readLenArr   = np.array(readLen)
        readScoreArr = np.array(readScore)

        readQ20 = readLenArr[readScoreArr < 20]
        plt.hist( readQ20, bins=1000, label='<Q20 ({:,} reads)'.format(len(readQ20)) )
        readQ30 = readLenArr[(readScoreArr < 30) & (readScoreArr >= 20)]
        plt.hist( readQ30, bins=1000, label='>=Q20, <Q30 ({:,} reads)'.format(len(readQ30)) )
        readQ50 = readLenArr[(readScoreArr < 50) & (readScoreArr >= 30)]
        plt.hist( readQ50, bins=1000, label='>=Q30, <Q50 ({:,} reads)'.format(len(readQ50)) )
        readQ50_ = readLenArr[(readScoreArr >= 50)]
        plt.hist( readQ50_, bins=1000, label='>=Q50 ({:,} reads)'.format(len(readQ50_)) )
        plt.xlim(0, 2000)
        plt.xlabel('Length (bp)')
        plt.ylabel('Reads No.')
        plt.legend()


        ax4 = plt.subplot(224)
        plt.boxplot(readLen, sym = '+')
        plt.xlabel('Reads')
        plt.ylabel('Length (bp)')
        plt.title("Max Length: {:,}".format(max(readLen)))

        plt.suptitle(sampleID)
        outdir = pathlib.Path('LOG-visnano').mkdir(parents=True, exist_ok=True)
        plt.savefig("LOG-visnano/summary-on-fq-%s.jpg" % (sampleID), dpi=300)
        # plt.show()
        return 

    baseScore, readLen, readScore = myFqInfo(fastqPath)
    print("Plotting...")
    myFqPlot(baseScore, readLen, readScore, sampleID)
    return 