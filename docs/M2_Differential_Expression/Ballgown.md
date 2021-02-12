# Ballgown

[Ballgown](http://www.bioconductor.org/packages/release/bioc/html/ballgown.html) is an R/Bioconductor package for flexible, isoform-level differential expression analysis of RNA-seq data [1]. Ballgown's data structures make it easy to use table-based packages like limma ([Smyth (2005)](https://github.com/alyssafrazee/ballgown/blob/master)), limma Voom ([Law et al. (2014)](http://dx.doi.org/10.1186/gb-2014-15-2-r29)), DESeq ([Anders & Huber (2010)](http://dx.doi.org/10.1186/gb-2010-11-10-r106)), DEXSeq ([Anders et al. (2012)](http://dx.doi.org/10.1101/gr.133744.111)), or edgeR ([Robinson et al. (2010)](http://dx.doi.org/10.1093/bioinformatics/btp616)) for differential expression analysis.

Ballgown requires three pre-processing steps:

1. RNA-Seq reads should be aligned to a reference genome. (**HISAT2**)
2. A transcriptome should be assembled, or a reference transcriptome should be downloaded. (**StringTie**)
3. Expression for the features (transcript, exon, and intron junctions) in the transcriptome should be estimated in a Ballgown readable format. (**StringTie**)



## 1. Convert SAM to BAM Using Samtools

Samtools is a suite of programs for interacting with high-throughput sequencing data [2]. SAM files produced by HISAT2 can be converted to BAM by samtools.

```shell
# Download and extract samtools
$ wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 -O samtools-1.11.tar.bz2
$ tar -xjvf samtools-1.11.tar.bz2

# Append to PATH environment variable
$ export PATH=$PATH:/path/to/samtools-1.11

# Verify installation
$ samtools --version
```

```shell
# Convert SAM to BAM
$ samtools view -S -b SRR5179431_trimmed.sam > SRR5179431_trimmed.bam

# Check the BAM
$ samtools view SRR5179431_trimmed.bam | head

# Sort the alignment
$ samtools sort SRR5179431_trimmed.bam -o SRR5179431_trimmed_s.bam

# Check the order
$ samtools view SRR5179431_trimmed_s.bam | head

# Index the sorted BAM
$ samtools index SRR5179431_trimmed_s.bam
```



## 2. Assemble Transcriptome and Estimate Feature Expression with StringTie

[StringTie](https://ccb.jhu.edu/software/stringtie/) is a highly efficient assembler for RNA-Seq alignments using a novel network flow algorithm [3]. It can simultaneously assemble and quantify expression levels for the features of the transcriptome in a Ballgown readable format. StringTie's output can be processed by specialized software like Ballgown( [Alyssa et al. (2014)](https://www.biorxiv.org/content/10.1101/003665v1.abstract)), Cuffdiff ([Cole et al. (2010)](https://www.nature.com/articles/nbt.1621)) or other programs (DESeq2 ([Anders & Huber (2010)](http://dx.doi.org/10.1186/gb-2010-11-10-r106)), edgeR ([Robinson et al. (2010)](http://dx.doi.org/10.1093/bioinformatics/btp616)), etc).

### Install StringTie

```shell
# Download and extract StringTie
$ wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.1.4.Linux_x86_64.tar.gz
$ tar xvfz stringtie-2.1.4.Linux_x86_64.tar.gz

# Append to PATH environment variable
$ 
# Verify installation
$ 
```



## 3. Ballgown Differential Expression Analysis

https://rnabio.org/module-03-expression/0003/03/01/Differential_Expression/







# Reference

[1] A. C. Frazee, G. Pertea, A. E. Jaffe, B. Langmead, S. L. Salzberg, and J. T. Leek, "Flexible isoform-level differential expression analysis with Ballgown," bioRxiv, p. 003665, 2014, doi: 10.1101/003665. [[paper](https://www.biorxiv.org/content/10.1101/003665v1.abstract)]

[2] H. Li, B. Handsaker, A. Wysoker, T. Fennell, J. Ruan, N. Homer et al., "The Sequence Alignment/Map format and SAMtools," (in eng), Bioinformatics, vol. 25, no. 16, pp. 2078-9, Aug 15 2009, doi: 10.1093/bioinformatics/btp352.[[paper](https://pubmed.ncbi.nlm.nih.gov/19505943/)]

[3] M. Pertea, G. M. Pertea, C. M. Antonescu, T.-C. Chang, J. T. Mendell, and S. L. Salzberg, "StringTie enables improved reconstruction of a transcriptome from RNA-seq reads," Nature Biotechnology, vol. 33, no. 3, pp. 290-295, 2015/03/01 2015, doi: 10.1038/nbt.3122. [[paper](https://www.nature.com/articles/nbt.3122)]