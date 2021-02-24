# Ballgown

[Ballgown](http://www.bioconductor.org/packages/release/bioc/html/ballgown.html) is an R/Bioconductor package for flexible, isoform-level differential expression analysis of RNA-seq data [1]. Ballgown's data structures make it easy to use table-based packages like limma ([Smyth (2005)](https://github.com/alyssafrazee/ballgown/blob/master)), limma Voom ([Law et al. (2014)](http://dx.doi.org/10.1186/gb-2014-15-2-r29)), DESeq ([Anders & Huber (2010)](http://dx.doi.org/10.1186/gb-2010-11-10-r106)), DEXSeq ([Anders et al. (2012)](http://dx.doi.org/10.1101/gr.133744.111)), or edgeR ([Robinson et al. (2010)](http://dx.doi.org/10.1093/bioinformatics/btp616)) for differential expression analysis.

Ballgown requires three pre-processing steps:

1. RNA-Seq reads should be aligned to a reference genome. (**HISAT2**)
2. A transcriptome should be assembled, or a reference transcriptome should be downloaded. (**StringTie**)
3. Expression for the features (transcript, exon, and intron junctions) in the transcriptome should be estimated in a Ballgown readable format. (**StringTie**)



## Differential Expression Analysis

An example of the working directory:

```markdown
homo_tables/
    SRR5978827/
        e2t.ctab
        e_data.ctab
        i2t.ctab
        i_data.ctab
        t_data.ctab
    SRR5978828/
        e2t.ctab
        e_data.ctab
        i2t.ctab
        i_data.ctab
        t_data.ctab
    ...
```

Create an R script `load_bg.R` for loading ballgown object:

```R
library(ballgown)
bg <- ballgown(dataDir="/path/to/homo_tables", samplePattern='SRR', meas='all')
save(bg, file='bg.rda')
```

Then run this script in terminal

```shell
$ R CMD BATCH load_bg.R
```

Differential expression analysis with Ballgown:

```R
# Set directory
setwd("/path/to/ballgown/")

# Load R packages
library(ballgown)
library(genefilter)

# Read in the data from StringTie
load("bg.rda")

# Load all attributes and gene names
bg_table = texpr(bg, 'all')
bg_gene_names = unique(bg_table[, 9:10])

# Add pData
# pData should hold a data frame of phenotype information for the samples in the experiment, 
# and be added during ballgown object construction. It can also be added later
group <- c(rep("iSLK-KSHV_BAC16-48hr-input", 3), rep("iSLK-uninf-input", 3))
pData(bg) = data.frame(id=sampleNames(bg), group=group)

# Perform differential expression (DE) analysis with no filtering
results_transcripts = stattest(bg, feature="transcript", covariate="group", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg, feature="gene", covariate="group", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_gene_names, by.x=c("id"), by.y=c("gene_id"))

# Filter low-abundance genes. Here we remove all transcripts with a variance across the samples of less than one
bg_filt = subset(bg,"rowVars(texpr(bg)) > 1", genomesubset=TRUE)

# Load all attributes including gene name
bg_filt_table = texpr(bg_filt , 'all')
bg_filt_gene_names = unique(bg_filt_table[, 9:10])

# Perform DE analysis now using the filtered data
results_transcripts = stattest(bg_filt, feature="transcript", covariate="group", getFC=TRUE, meas="FPKM")
results_genes = stattest(bg_filt, feature="gene", covariate="group", getFC=TRUE, meas="FPKM")
results_genes = merge(results_genes, bg_filt_gene_names, by.x=c("id"), by.y=c("gene_id"))

# Identify the significant genes with p-value < 0.05
sig_transcripts = subset(results_transcripts, results_transcripts$pval<0.05)
sig_genes = subset(results_genes, results_genes$pval<0.05)
```

## Results

DE results for dataset GSE93676:

```R
# To see the first few lines of significant genes and transcripts:
head(sig_transcripts)
#       feature  id           fc        pval      qval
# 4  transcript  40 1.237362e-02 0.047948681 0.7342077
# 13 transcript  68 3.643648e-13 0.041800485 0.7342077
# 35 transcript 133 2.283918e+17 0.009168513 0.7342077
# 51 transcript 167 6.634653e-10 0.034473953 0.7342077
# 70 transcript 244 2.241706e-16 0.048610761 0.7342077
# 76 transcript 276 2.637682e-18 0.002536443 0.7342077

head(sig_genes)
#              id feature           fc       pval      qval gene_name
# 1     100128071    gene 1.445238e-59 0.04327419 0.7012609   FAM229A
# 16       119710    gene 1.378979e-51 0.01122029 0.6171722  C11orf74
# 57         6289    gene 1.718444e-71 0.01229523 0.6171722      SAA2
# 74  MSTRG.10000    gene 1.942633e+58 0.01965612 0.6250552     NGLY1
# 97  MSTRG.10028    gene 5.764237e+05 0.04138674 0.6962953   PDCD6IP
# 141 MSTRG.10096    gene 1.998718e-66 0.02531592 0.6422143     PTH1R
```

DE results for dataset GSE47217:

```R
head(sig_transcripts)
#        feature   id        fc        pval      qval
# 57  transcript  479 0.9339480 0.036290263 0.7220780
# 66  transcript  522 0.4784699 0.002714966 0.5205885
# 86  transcript  685 0.8718279 0.008432731 0.5794108
# 107 transcript  847 0.8360476 0.031839607 0.7176861
# 117 transcript  908 1.0784123 0.004445226 0.5428939
# 137 transcript 1051 0.4439950 0.036100181 0.7220780

head(sig_genes)
#             id feature        fc       pval      qval gene_name
# 1   MSTRG.1000    gene 1.0657124 0.04960470 0.5469039  Ppp1r15b
# 36 MSTRG.10096    gene 0.9003980 0.03641091 0.5417405     Lpin2
# 55 MSTRG.10153    gene 0.9036117 0.02278910 0.5227679     Birc6
# 56 MSTRG.10153    gene 0.9036117 0.02278910 0.5227679         .
# 85  MSTRG.1024    gene 1.2618979 0.03277910 0.5417405  Tmem183a
# 94  MSTRG.1027    gene 0.9601584 0.04558797 0.5458907   Adipor1
```



# Reference

[1] A. C. Frazee, G. Pertea, A. E. Jaffe, B. Langmead, S. L. Salzberg, and J. T. Leek, "Flexible isoform-level differential expression analysis with Ballgown," bioRxiv, p. 003665, 2014, doi: 10.1101/003665. [[paper](https://www.biorxiv.org/content/10.1101/003665v1.abstract)]
