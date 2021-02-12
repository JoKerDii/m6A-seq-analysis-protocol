# Trim Galore

Trim Galore is a Perl wrapper around Cutadapt and FastQC to consistently apply adapter and quality trimming to FastQ files.



## Install Trim Galore

Before installation, ensure [Cutadapt](https://github.com/marcelm/cutadapt) and [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) are already installed.

```shell
# Check the version of cutadapt
$ cutadapt --version

# Check the version of FastQC
$ fastqc -v
```

Install the latest version of Trim Galore from [Github](https://github.com/FelixKrueger/TrimGalore) or [project website](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/):

```shell
# Install Trim Galore
$ curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz
$ tar xvzf trim_galore.tar.gz

# Verify installation
$ trim_galore -v
```

```markdown
# Usage
trim_galore [options] <filename(s)>
```



## Adaptive Quality and Adapter Trimming

First, low-quality base calls are trimmed off from the 3' end of the reads before adapter removal. Next, adapter sequences from the 3’ end of reads are detected and removed by cutadapt. Lastly, trimmed short sequences (default: < 20bp) are filtered.

```shell
$ trim_galore -o /path/to/TrimGlore_result/ SRR5179431.fastq
```

**Note 1**: `-o` (or `--outdir`) will create all output files in the specified output directory. 

**Note 2**: Please refer to the [User Guide](https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md) about other options for specific purposes and conditions. 



## Outputs

Trim Galore produced two output files for each FastQ file: one text file ("SRR5179431.fastq_trimming_report.txt") and a trimmed FastQ file ("SRR5179431_trimmed.fq"). 

### 1. The text file

The text file provides a summary of running parameters.

```shell
# To see the first few lines of the text file
$ head SRR5179431.fastq_trimming_report.txt
```

```markdown
SUMMARISING RUN PARAMETERS
==========================
Input filename: /data/zhendi/protocol/SRR5179431.fastq
Trimming mode: single-end
Trim Galore version: 0.6.4_dev
Cutadapt version: 3.2
Number of cores used for trimming: 1
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33
(base) zhendi@amax:/data/zhendi/protocol/trim_galore_result$ head SRR5179431.fastq_trimming_report.txt

SUMMARISING RUN PARAMETERS
==========================
Input filename: /data/zhendi/protocol/SRR5179431.fastq
Trimming mode: single-end
Trim Galore version: 0.6.4_dev
Cutadapt version: 3.2
Number of cores used for trimming: 1
Quality Phred score cutoff: 20
Quality encoding type selected: ASCII+33

```

### 2. The trimmed FastQ

The trimmed FastQ file can be used for further analysis.

# Reference

Trim Galore User Guide on Github: https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md

