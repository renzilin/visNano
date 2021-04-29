# visNano 0.1.1
A tool for making summary of long reads. This tools can handle FASTQ/BAM files.

## Functions:
1. Work on FASTQ file

![image1](https://github.com/renzilin/visNano/blob/master/figs/summary-on-fq-F36.jpg)

2. Work on BAM file

![image2](https://github.com/renzilin/visNano/blob/master/figs/summary-on-sam-F36.jpg)

3. Work on BAM+BED file

![image3](https://github.com/renzilin/visNano/blob/master/figs/summary-on-bed-F36.jpg)



## Usage:
```bash
python visNano.py -h
```

Or,

Please see the examply.ipynb



## Install:

```bash
conda create -n visnano
conda activate visnano
conda install numpy pandas tqdm matplotlib samtools
```

```
python packages:
    - numpy
    - pandas
    - tqdm
    - matplotlib
    - samtools
```

## NOTE:
This tool is still under development.
If you want to add some functions, please send a [email](zilin.ren@outlook.com) to me
