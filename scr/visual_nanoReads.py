#!/usr/bin/env python
# coding: utf-8

"""

- Visual_nanoReads ver1.0
- Written by Ren Zilin
- Function: make a summary for fastq file and bam file

"""


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
#
# 1. alignment distribution
# 2. alignment read length
# 3. alignment distribution in region
#

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
import matplotlib.pyplot as plt
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


import pysam


# In[15]:


parser = argparse.ArgumentParser(description='Make a summary for long reads. Still working on')
parser.add_argument('--fq', help='summary about fastq file (not *.gz file)', default = None)
parser.add_argument('--bam', help='summary about bam file', default = None)
parser.add_argument('--bed', help='summary about bed file', default = None)
args = parser.parse_args()

fastq_path = args.fq
bam_path = args.bam
bed_path = args.bed


# In[27]:


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

def cal_CGcont(count_dict):
    return (count_dict['C'] + count_dict['G'])/(count_dict['A'] + count_dict['T']
                                                + count_dict['C'] + count_dict['G'])
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

        fig.savefig('%s/%s_read_length_dist.png' % (OUTDIR, tag), dpi=300)
        plt.close(fig)

        return

    def score_dist_plot(read_mean_scores, file_tag, OUTDIR):

        fig = plt.figure(figsize=(6.2, 4.2))
        sns.boxplot(y = read_mean_scores, width=.3)

        plt.xlabel(file_tag)
        plt.ylabel('Mean score of reads')
        plt.title('The distribution of mean read score (%s reads)' % "{:,}".format(len(read_mean_scores)))

        fig.savefig('%s/read_score_dist.png' % OUTDIR, dpi=300)
        plt.close(fig)

        return

    def length_score_plot(plotDF, OUTDIR):
        jp = sns.jointplot("read length", "read score", data=plotDF, kind='kde', space=0)
        jp.savefig('%s/read_length_score.png' % OUTDIR, dpi=300)
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

    fig.savefig('%s/head_score_dist.png' % OUTDIR, dpi=300)
    plt.close(fig)

    ## tail
    fig = plt.figure(figsize=(6.2, 4.2))

    pd.DataFrame(read_tail_scores,
                 columns=[str(i+1) for i in range(len(read_head_scores[0]))]).boxplot(grid=False, flierprops=flierprops)
    plt.ylim(10,90)
    plt.title('The distribution of reads score (tail)')
    plt.xlabel('base coordinate')
    plt.ylabel('phred score')

    fig.savefig('%s/tail_score_dist.png' % OUTDIR, dpi=300)
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
    fig.savefig('%s/bases_dist.png' % (OUTDIR), dpi=300)
    plt.close(fig)

    return


# In[8]:
def extract_bam(bam_path):
    bam_list = []
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    for line in bamfile.fetch():
        if not line.query:  ## second alignment is omitted
            continue
        bam_list.append([line.query_name, line.query_length, line.reference_name,
                         line.reference_start, line.reference_end, line.query_alignment_length,
                         line.flag, line.get_tag('AS'), line.mapping_quality,
                         cal_CGcont(get_CGcont(line.query))])
    bamfile.close()

    ## data frame
    bam_dataframe = pd.DataFrame(bam_list)
    bam_dataframe.columns = ['read_name', 'read_length', 'reference_name', 'ref_start',
                             'ref_end', 'query_alignment_length', 'FALG', 'DP_score', 'MAPQ',
                             'CGcont']

    return bam_dataframe

# In[9]:
def align_dist_plot(bam_info, OUTDIR):

    def create_chr_index(reference_name):
        search_dict = {}
        for i in range(22):
            search_dict['chr%s' % (i+1)] = i
        search_dict['chrX'] = 22
        search_dict['chrY'] = 23
        search_dict['chrM'] = 24
        record_notin= 25

        sort_index = []
        sort_map_ind = []
        sort_map_chr = []

        for i in reference_name:
            if i not in search_dict:
                search_dict[i] = record_notin
                record_notin += 1
            sort_index.append(search_dict[i])

        sort_map_ind = list(range(record_notin))
        sort_map_chr = list(search_dict.keys())

        return sort_index, sort_map_ind, sort_map_chr


    sort_index, sort_map_ind, sort_map_chr = create_chr_index(bam_info.loc[:,'reference_name'])

    bam_info['sort_index'] = sort_index


    ## alignment length dist
    fig = plt.figure(figsize=(6, 6))
    sns.boxplot(x="sort_index", y="read_length", data=bam_info)
    plt.xlabel('chromosome')
    plt.ylabel('length of reads')
    plt.title('The distribution of alignment length')
    plt.xticks(sort_map_ind, sort_map_chr, rotation='vertical')
    fig.tight_layout()

    fig.savefig('%s/align__length_dist.png' % (OUTDIR), dpi=300)
    plt.close(fig)


    ## alignemnt count dist
    align_info = pd.DataFrame(bam_info['sort_index'].value_counts())
    align_info.reset_index(inplace=True)
    align_info.columns = ['sort_index', 'counts']

    ##
    fig = plt.figure(figsize=(6, 6))
    sns.barplot(x="sort_index", y="counts", data=align_info)
    plt.xlabel('chromosome')
    plt.ylabel('num of reads')
    plt.title('The distribution of alignment (%s reads)' % "{:,}".format(sum(align_info['counts'])))
    plt.xticks(sort_map_ind, sort_map_chr, rotation='vertical')
    fig.tight_layout()

    fig.savefig('%s/align_count_dist.png' % (OUTDIR), dpi=300)
    plt.close(fig)

    ## alignemnt log_count dist
    align_info = pd.DataFrame(bam_info['sort_index'].value_counts())
    align_info.reset_index(inplace=True)
    align_info.columns = ['sort_index', 'counts']

    ##
    fig = plt.figure(figsize=(6, 6))
    sns.barplot(x="sort_index", y="counts", data=align_info)
    plt.xlabel('chromosome')
    plt.ylabel('num of reads (log)')
    plt.title('The distribution of alignment (%s reads)' % "{:,}".format(sum(align_info['counts'])))
    plt.xticks(sort_map_ind, sort_map_chr, rotation='vertical')
    plt.yscale('log')
    fig.tight_layout()

    fig.savefig('%s/align_log_count_dist.png' % (OUTDIR), dpi=300)
    plt.close(fig)

    return


# In[24]:
def region_func(region_name, region_interest, bam_info, OUTDIR):

    ##

    def region_plots(bam_info, info_dfs, OUTDIR):
        ## region dist
        region_dist = pd.DataFrame(info_dfs['region_index'].value_counts())
        region_dist.reset_index(inplace=True)
        region_dist.columns = ['regions', 'counts']
        region_dist['sort_ind'] = region_dist['regions'].str.extract('(\d+)', expand=False).astype(int)
        region_dist.sort_values(by=['sort_ind'], inplace=True, ascending=True)
        region_dist.drop('sort_ind', axis=1, inplace=True)

        ##
        fig = plt.figure(figsize=(6, 6))
        xbar = np.arange(len(region_dist))
        plt.bar(xbar, region_dist['counts'])
        plt.xticks(xbar, region_dist['regions'].values, rotation=90)

        # for i, v in enumerate(region_dist['counts'].values):
        #     plt.annotate('%s / %s' % (v, len(bam_info)), (i, v), ha='center')

        plt.yscale('log')
        plt.xlabel('region interested')
        plt.ylabel('num of reads (log scale)')
        plt.title('The number of reads in regions')
        fig.tight_layout()

        fig.savefig('%s/region_align_dist_logscale.png' % (OUTDIR), dpi=300)
        plt.close(fig)


        ##
        fig = plt.figure(figsize=(6, 6))
        xbar = np.arange(len(region_dist))
        plt.bar(xbar, region_dist['counts'])
        plt.xticks(xbar, region_dist['regions'].values, rotation=90)

        plt.xlabel('region interested')
        plt.ylabel('num of reads')
        plt.title('The number of reads in regions')
        fig.tight_layout()

        fig.savefig('%s/region_align_dist.png' % (OUTDIR), dpi=300)
        plt.close(fig)

        ## region length boxplot
        fig = plt.figure(figsize=(6, 6))
        sns.boxplot(data=info_dfs, x='region_index', y='read_length', fliersize=.5)
        plt.xlabel('region interested')
        plt.ylabel('length of reads')
        plt.title('The length of reads in regions')
        plt.xticks(rotation=90)
        fig.tight_layout()

        fig.savefig('%s/region_length_dist.png' % (OUTDIR), dpi=300)
        plt.close(fig)

        ## region MAPQ boxplot
        fig = plt.figure(figsize=(6, 6))
        sns.boxplot(data=info_dfs, x='region_index', y='DP_score', fliersize=.5)
        plt.xlabel('region interested')
        plt.ylabel('DP alignment score')
        plt.title('The Alignment Score in regions')
        plt.xticks(rotation=90)
        fig.tight_layout()

        fig.savefig('%s/region_AS_dist.png' % (OUTDIR), dpi=300)
        plt.close(fig)

        ## region DP alignment score boxplot
        fig = plt.figure(figsize=(6, 6))
        sns.boxplot(data=info_dfs, x='region_index', y='MAPQ', fliersize=.5)
        plt.xlabel('region interested')
        plt.ylabel('Mapping Quality')
        plt.title('The MAPQ in regions')
        plt.xticks(rotation=90)
        fig.tight_layout()

        fig.savefig('%s/region_MAPQ_dist.png' % (OUTDIR), dpi=300)
        plt.close(fig)

                ## region CG content boxplot
        fig = plt.figure(figsize=(6, 6))
        sns.boxplot(data=info_dfs, x='region_index', y='CGcont', fliersize=.5)
        plt.xlabel('region interested')
        plt.ylabel('CG concent')
        plt.title('The CG content in regions')
        plt.xticks(rotation=90)
        fig.tight_layout()

        fig.savefig('%s/region_CGcontent_dist.png' % (OUTDIR), dpi=300)
        plt.close(fig)

        return

    ##

    def each_region_check(region_name, region_interest, bam_info):
        region_search  = re.search('([^:]*):([0-9]*)-([0-9]*)', region_interest)
        chrom_interest = region_search.group(1)
        start_interest = int(region_search.group(2))
        end_interest = int(region_search.group(3))

        region_info = bam_info.loc[(bam_info['reference_name'] == chrom_interest)
                                   & (bam_info['ref_start'] <= start_interest)
                                   & (bam_info['ref_end'] >= end_interest)].copy(deep=True)

        region_info['region_index'] = region_name

        return region_info

    ##

    info_dfs = []
    for i, each_region in enumerate(region_interest):
        info_dfs.append(each_region_check(region_name[i], each_region, bam_info))

    info_dfs = pd.concat(info_dfs)

    region_plots(bam_info, info_dfs, OUTDIR)

    return


# In[21]:
def bed_file_reader(bed_path):
    region_name = []
    region_interest = []
    with open(bed_path, 'r') as file:
        for i, line in enumerate(file):
            line_list = line.strip().split()
            if len(line_list) < 4:
                region_name.append('loc_%s' % i)
            else:
                region_name.append('loc_%s/%s' % (i, line_list[3]))
            region_interest.append('%s:%s-%s' % (line_list[0],line_list[1],line_list[2]))
    return region_name, region_interest




# In[7]:



def main_func_fq(fastq_path, flank_len = 20):

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

    if not os.path.exists('Visual_nanoReads_logfqs'):
        os.mkdir('Visual_nanoReads_logfqs')

    if not os.path.exists('Visual_nanoReads_logfqsVisual_nanoReads_logfqs/%s' % file_tag):
        os.mkdir('Visual_nanoReads_logfqs/%s' % file_tag)

    OUTDIR = 'Visual_nanoReads_logfqs/%s' % file_tag

    length_plot(read_lengths, read_mean_scores, file_tag, OUTDIR)
    plot_head_tail(read_head_scores, read_tail_scores, OUTDIR)
    CGcontent_plot(read_CGconts, OUTDIR)

    return

def main_func_bam(bam_path, bed_path):

    file_tag = re.search('.*/([^/.]*)', bam_path).group(1)

    if not os.path.exists('Visual_nanoReads_logbams'):
        os.mkdir('Visual_nanoReads_logbams')

    if not os.path.exists('Visual_nanoReads_logbams/%s' % file_tag):
        os.mkdir('Visual_nanoReads_logbams/%s' % file_tag)

    OUTDIR = 'Visual_nanoReads_logbams/%s' % file_tag


    bam_info = extract_bam(bam_path)

    ## alignment distribution
    align_dist_plot(bam_info, OUTDIR)

    if bed_path != None:
        ## alignment in concerned region
        region_name, region_interest = bed_file_reader(bed_path)
        region_func(region_name, region_interest, bam_info, OUTDIR)

    return






# In[22]:


## input1: combined fastq
# fastq_path = './data/fastq_pass/barcode01.fastq'
if fastq_path is not None:
    main_func_fq(fastq_path)

## TODO: input2: barcode dir
# fastq_dir = '../data/fastq_pass/barcode01'


# In[22]:
if bam_path is not None:
    print(bam_path)
    print(bed_path)
    main_func_bam(bam_path, bed_path)

"""
bam_path = './alignment/bams/barcode01.bam'
bed_path = './data/M20022-20200318/M20022-20200318/loci-gene.bed'



file_tag = re.search('.*/([^/.]*)', bam_path).group(1)

if not os.path.exists('Visual_nanoReads_logbams'):
    os.mkdir('Visual_nanoReads_logbams')

if not os.path.exists('Visual_nanoReads_logbams/%s' % file_tag):
    os.mkdir('Visual_nanoReads_logbams/%s' % file_tag)

OUTDIR = 'Visual_nanoReads_logbams/%s' % file_tag


# In[13]:


bam_info = extract_bam(bam_path)

##
align_dist_plot(bam_info, OUTDIR)

##
region_name, region_interest = bed_file_reader(bed_path)
region_func(region_name, region_interest, bam_info, OUTDIR)


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:
"""
