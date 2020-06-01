# Visual_nanoReads
A tool for making summary of long reads. This tools can handle FASTQ/BAM files.

## Functions:
- FASTQ:
    - total # of bases
    - total # of reads
    - length distribution in FASTQ
    - all reads head/tail score distribution (default: 20 bp)
    - length distribution with mean phred score
    - distribution of read mean score
    - CG content
- BAM:
    - distribution of alignment length (read length)
    - distribution of alignment counts (read count)
    - distribution of alignment in specific region (REQUIRE a bed file with 4 cols: CHROM START END NAME)


## Usage:
```bash
python XXX/scr/visual_nanoReads.py --fq $FASTQ

python XXX/scr/visual_nanoReads.py --bam $BAM

python XXX/scr/visual_nanoReads.py --bam $BAM --bed $BED
```

## Install:
```
python packages:
    - numpy
    - pandas
    - tqdm
    - seaborn
    - matplotlib
    - pysam
```

## NOTE:
This tool is still under development.
If you want to add some functions, please send a email to me. zilin.bj@gmail.com
