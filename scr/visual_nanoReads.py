#!/usr/bin/env python
# coding: utf-8

# ## 1. visualization on nanopore fastq data
# 
# ### 1.1 work on each line
# - read info
# - each read length -> 1.2 total # of bases, total # of reads, length distribution, max length, min length
# - each read nucletide percentage -> 1.2 CG content (all output is space-consuming)
# - head/tail quality -> 1.2 all reads head/tail score distribution
# - mean read quality -> 1.2 distribution of mean score, # of reads with poor score
# 
# ### 1.2 work on total
# 1. total # of bases
# 2. total # of reads
# 3. length distribution
# 4. max/min length
# 5. all reads head/tail score distribution
# 6. length distribution with mean score
# 7. distribution of mean score
# 8. CG content
# 
# 
# ## 2. visualization on nanopore bam data

# In[1]:


import os
import re

import argparse

import numpy as np
import pandas as pd

import collections as cl

from tqdm import tqdm

import seaborn as sns
import matplotlib as mpl

mpl.rcParams.update({'font.size': 11,
                     'axes.titlesize': 15,     # title font size
                     'axes.labelsize': 13,     # xlabel/ylabel font size
                     'axes.linewidth': 1.5,    # axis width
                     'xtick.labelsize': 11,    # xticks font size
                     'ytick.labelsize': 11,    # yticks font size
                     'xtick.major.size': 6,    # xticks length
                     'xtick.major.width': 1.5, # xticks width
                     'ytick.major.size': 6,    # yticks length
                     'ytick.major.width': 1.5, # yticks width
                     })

import matplotlib.pyplot as plt


# In[ ]:


parser = argparse.ArgumentParser(description='Make a summary for long reads. Still working on')
parser.add_argument('--fq', help='summary about fastq file (not *.gz file)')
args = parser.parse_args()

fastq_path = args.fq


# In[2]:


## work on each line
def get_info(line):
    read_info = re.search('@([^\s]*).*', line)
    read_name = read_info.group(1)
    return read_name


def get_length(line):
    return len(line)


def score_dist(line, flank_len):
    ## base score
    read_score = np.array([ord(i) for i in line])
    
    ## read mean score
    read_mean_score = np.mean(read_score)
    
    ## head & tail
    head_score = read_score[0:flank_len]
    tail_score = read_score[-flank_len:]
    
    return read_mean_score, head_score, tail_score

def get_CGcont(line):
    count_dict = cl.Counter(line)
    return count_dict


# In[3]:


## work on total
def print_info(read_lengths, read_mean_scores):
    return 1


# In[4]:


def length_plot(read_lengths, read_mean_scores, file_tag, OUTDIR):
    def length_dist_plot(length_list, tag, OUTDIR):

        fig = plt.figure(figsize=(6.2, 4.2))
        sns.distplot(length_list, kde=False, bins=int(len(length_list)/200))
        plt.xlabel('The length of reads')
        plt.ylabel('The number of reads')
        plt.title('The distribution of %s read length (%s reads)' % (tag, "{:,}".format(len(read_lengths))))
        plt.annotate('Max length: %s\nMin length: %s' % ("{:,}".format(max(length_list)), "{:,}".format(min(length_list))), 
                     (.7,.7), xycoords = 'axes fraction')

        fig.savefig('%s/%s_read_length_dist.pdf' % (OUTDIR, tag), dpi=300)
        plt.close(fig)

        return

    def score_dist_plot(read_mean_scores, file_tag, OUTDIR):

        fig = plt.figure(figsize=(6.2, 4.2))
        sns.boxplot(y = read_mean_scores, width=.3)

        plt.xlabel(file_tag)
        plt.ylabel('Mean score of reads')
        plt.title('The distribution of mean read score (%s reads)' % "{:,}".format(len(read_mean_scores)))

        fig.savefig('%s/read_score_dist.pdf' % OUTDIR, dpi=300)
        plt.close(fig)

        return

    def length_score_plot(plotDF, OUTDIR):
        jp = sns.jointplot("read length", "read score", data=plotDF, kind='kde', space=0)
        jp.savefig('%s/read_length_score.pdf' % OUTDIR, dpi=300)
        return
    
    read_lengths     = np.array(read_lengths)
    read_mean_scores = np.array(read_mean_scores)

    total_bases = sum(read_lengths)
    total_reads = len(read_lengths)
    print('the total bases is %s' % "{:,}".format(total_bases))
    print('the total reads is %s' % "{:,}".format(total_reads))
    
    ## all length dist
    length_dist_plot(read_lengths, 'all', OUTDIR)
    

    ## N20 length dist
    total_N20_reads   = read_lengths[read_mean_scores > 20]
    length_dist_plot(total_N20_reads, 'N20', OUTDIR)

    ## N30 length dist
    total_N30_reads   = read_lengths[read_mean_scores > 30]
    length_dist_plot(total_N30_reads, 'N30', OUTDIR)

    ## N50 length dist
    total_N50_reads   = read_lengths[read_mean_scores > 50]
    length_dist_plot(total_N50_reads, 'N50', OUTDIR)
    
    ## length and score
    plotDF = pd.concat([pd.DataFrame(read_lengths, columns=['read length']), pd.DataFrame(read_mean_scores, columns=['read score'])], axis=1)
    length_score_plot(plotDF, OUTDIR)
    
    ## read score dist
    score_dist_plot(read_mean_scores, file_tag, OUTDIR)
    
    return


# In[5]:


def plot_head_tail(read_head_scores, read_tail_scores, OUTDIR):
    
    flierprops = dict(marker='.', markerfacecolor='k', markersize=5,
                  linestyle='none', markeredgecolor='none')

    ## head
    fig = plt.figure(figsize=(6.2, 4.2))

    pd.DataFrame(read_head_scores, 
                 columns=[str(i+1) for i in range(len(read_head_scores[0]))]).boxplot(grid=False, flierprops=flierprops)
    plt.ylim(10,90)
    plt.title('The distribution of reads score (head)')
    plt.xlabel('base coordinate')
    plt.ylabel('phred score')
    
    fig.savefig('%s/head_score_dist.pdf' % OUTDIR, dpi=300)
    plt.close(fig)
    
    ## tail
    fig = plt.figure(figsize=(6.2, 4.2))

    pd.DataFrame(read_tail_scores, 
                 columns=[str(i+1) for i in range(len(read_head_scores[0]))]).boxplot(grid=False, flierprops=flierprops)
    plt.ylim(10,90)
    plt.title('The distribution of reads score (tail)')
    plt.xlabel('base coordinate')
    plt.ylabel('phred score')

    fig.savefig('%s/tail_score_dist.pdf' % OUTDIR, dpi=300)
    plt.close(fig)
    
    return


# In[6]:


def CGcontent_plot(read_CGconts, OUTDIR):
    CGtable = read_CGconts[0]
    for i in read_CGconts[1:]:
        CGtable += i
    
    
    CG_print = 'The CG content is %.3f %%' % (100*(CGtable['C'] + CGtable['G'])/(CGtable['C'] + CGtable['G'] + CGtable['A'] + CGtable['T']))

    print(CG_print)
    
    CGtable_ = []
    for i in CGtable:
        CGtable_.append([i, CGtable[i]])
    
    CGtable_ = []
    for i in CGtable:
        CGtable_.append([i, CGtable[i]])

    CGtable_df = pd.DataFrame(CGtable_, columns=['base', 'count'])
    CGtable_df['count_log'] = np.log10(CGtable_df['count'])

    fig = plt.figure(figsize=(6.2, 4.2))
    sns.barplot(x="base", y="count_log", data=CGtable_df)
    plt.xlabel('base')
    plt.ylabel('log10(Counts)')
    plt.title('The counts of bases (%s base)' % "{:,}".format(sum(CGtable_df['count'])))
    plt.ylim(0, max(CGtable_df['count_log'])*1.25)
    plt.annotate(CG_print, (.5, .9), xycoords = 'axes fraction')
    # 
    fig.savefig('%s/bases_dist.pdf' % (OUTDIR), dpi=300)
    plt.close(fig)
    
    return


# In[7]:



def main_func(fastq_path, flank_len = 20):

    read_names   = []
    read_lengths = []
    read_CGconts = []
    read_mean_scores = []
    read_head_scores = []
    read_tail_scores = []

    with open(fastq_path, 'r') as file:
        for line_num, line in tqdm(enumerate(file), ascii = True, 
                                   desc = 'handling each line', ncols=100):
            clean_line = line.strip()

            if line_num % 4 == 0:
                next
    #             read_names.append(get_info(clean_line))
            elif line_num % 4 == 1:
                read_lengths.append(get_length(clean_line))
                read_CGconts.append(get_CGcont(clean_line))
            elif line_num % 4 == 3:
                score_tuples = score_dist(clean_line, flank_len)
                read_mean_scores.append(score_tuples[0])
                read_head_scores.append(score_tuples[1])
                read_tail_scores.append(score_tuples[2])

                
    file_tag = re.search('.*/([^/.]*)', fastq_path).group(1)

    if not os.path.exists('nanopore_vis_logs'):
        os.mkdir('nanopore_vis_logs')

    if not os.path.exists('nanopore_vis_logs/%s' % file_tag):
        os.mkdir('nanopore_vis_logs/%s' % file_tag)

    OUTDIR = 'nanopore_vis_logs/%s' % file_tag

    length_plot(read_lengths, read_mean_scores, file_tag, OUTDIR)
    plot_head_tail(read_head_scores, read_tail_scores, OUTDIR)
    CGcontent_plot(read_CGconts, OUTDIR)
    
    return


# In[9]:


## input1: combined fastq
# fastq_path = '../data/fastq_pass/barcode01.fastq'
if fastq_path is not None:
    main_func(fastq_path)




## TODO: input2: barcode dir
# fastq_dir = '../data/fastq_pass/barcode01'


# In[ ]:





# In[ ]:




