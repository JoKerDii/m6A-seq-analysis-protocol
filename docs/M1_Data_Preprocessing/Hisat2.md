# HISAT2

[HISAT2](http://daehwankimlab.github.io/hisat2/manual/) is a fast and sensitive alignment program for mapping next-generation sequencing reads to reference genome(s) [1,2].



## Install HISAT2

```shell
# Download and extract the latest version
$ wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.0.4-Linux_x86_64.zip
$ unzip hisat2-2.0.4-Linux_x86_64.zip

# Append to PATH environment variable
$ export PATH=$PATH:/path/to/hisat2-2.0.4

# Verify installation
$ hisat2 --help
$ hisat2 --version
```



## Sequence Alignment

### 1. Build indexes

You may need reference sequence, and gene annotation to build indexes. 

```shell
# You can either download HISAT2 indexes from its website: http://daehwankimlab.github.io/hisat2/download/
$ wget https://genome-idx.s3.amazonaws.com/hisat/hg19_genome.tar.gz
$ tar -zxvf hg19_genome.tar.gz

# or download reference sequence and gene annotation from Illumina iGenome before building index by hisat2-build
$ wget https://s3.amazonaws.com/local-run-manager/genomes/Homo_sapiens.zip
$ unzip Homo_sapiens.zip
# genome in: /path/to/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
$ mkdir /path/to/index
$ cd /path/to/index
$ hisat2-build -p 20 /path/to/genome.fa genome
```

 `hisat2-build` generates eight `.ht2` files, from `genome.1.ht2` to `genome.8.ht2`.

### 2. Run HISAT2 to get SAM output

```shell
$ hisat2 -p 10 -x genome -U /path/to/SRR5179431_trimmed.fq -S SRR5179431_trimmed.sam --summary-file SRR5179431_trimmed_summary.txt
```

**Note 1**:

| Common Arguments and Options | Description |
| ------ | ----------- |
| -x <hisat2-idx> | The basename of the index for the reference genome. |
| -U <r> | Comma-separated list of files containing unpaired reads to be aligned. |
| -S <hit> | File to write SAM alignments to. |
| --summary-file    | Print alignment summary to this file. |
| -p/--threads NTHREADS | Launch NTHREADS parallel search threads (default: 1). |

**Note 2**: HISAT2 does not require a GFF or GTF file but these annotation files will be needed at a later stage of analysis.

### Outputs

Two files should be generated: one SAM file ("SRR5179431_trimmed.sam") and a text file ("SRR5179431_trimmed_summary.txt").

```shell
# To see the first few lines of the output SAM file
$ head SRR5179431_trimmed.sam
```

```markdown
@HD     VN:1.0  SO:unsorted
@SQ     SN:chrM LN:16571
@SQ     SN:chr1 LN:249250621
@SQ     SN:chr2 LN:243199373
@SQ     SN:chr3 LN:198022430
@SQ     SN:chr4 LN:191154276
@SQ     SN:chr5 LN:180915260
@SQ     SN:chr6 LN:171115067
@SQ     SN:chr7 LN:159138663
@SQ     SN:chr8 LN:146364022
```

```shell
# To see the first few lines of the output text file
$ head SRR5179431_trimmed_summary.txt
```

```markdown
9758377 reads; of these:
  9758377 (100.00%) were unpaired; of these:
    1191861 (12.21%) aligned 0 times
    7985060 (81.83%) aligned exactly 1 time
    581456 (5.96%) aligned >1 times
87.79% overall alignment rate
```



# Reference

[1] D. Kim, J. M. Paggi, C. Park, C. Bennett, and S. L. Salzberg, "Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype," Nature Biotechnology, vol. 37, no. 8, pp. 907-915, 2019/08/01 2019, doi: 10.1038/s41587-019-0201-4. [[paper](https://pubmed.ncbi.nlm.nih.gov/31375807/)]

[2] HISAT2 Manual on Github: http://daehwankimlab.github.io/hisat2/manual/