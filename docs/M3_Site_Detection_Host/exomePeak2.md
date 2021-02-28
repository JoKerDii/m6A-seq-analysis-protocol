# exomePeak2

exomePeak2 is an R/Bioconductor package which provides bias-aware quantification and peak detection for Methylated RNA immunoprecipitation sequencing data (MeRIP-Seq) [1]. We are going to use this package for peak calling, finding enriched m6A sites.



## Installation

```R
if(!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("exomePeak2")
```



## Peak Calling

```R
library(exomePeak2)
set.seed(1)
root = "/path/to/homo_result"
setwd(root)

f1 = file.path(root, "SRR5978834_sorted.bam")
f2 = file.path(root, "SRR5978835_sorted.bam")
f3 = file.path(root, "SRR5978836_sorted.bam")
IP_BAM = c(f1,f2,f3) 

f1 = file.path(root, "SRR5978827_sorted.bam")
f2 = file.path(root, "SRR5978828_sorted.bam")
f3 = file.path(root, "SRR5978829_sorted.bam")
INPUT_BAM = c(f1,f2,f3)

exomePeak2(bam_ip = IP_BAM,
           bam_input = INPUT_BAM,
           genome = "hg19",
           library_type = "1st_strand",
           paired_end = FALSE)
```

```R
library(exomePeak2)
set.seed(1)
root = "/path/to/hhv8_result"
setwd(root)

f1 = file.path(root, "SRR5978834_sorted.bam")
f2 = file.path(root, "SRR5978835_sorted.bam")
f3 = file.path(root, "SRR5978836_sorted.bam")
IP_BAM = c(f1,f2,f3) 

f1 = file.path(root, "SRR5978827_sorted.bam")
f2 = file.path(root, "SRR5978828_sorted.bam")
f3 = file.path(root, "SRR5978829_sorted.bam")
INPUT_BAM = c(f1,f2,f3)

GENE_ANNO_GTF = file.path("/path/to/sequence.gff3")

exomePeak2(bam_ip = IP_BAM,
           bam_input = INPUT_BAM,
           gff_dir = GENE_ANNO_GTF,
           library_type = "1st_strand",
           paired_end = FALSE)
```

```R
library(exomePeak2)
set.seed(1)
root = "/path/to/mm10_result"
setwd(root)

f1 = file.path(root, "SRR866997_sorted.bam")
f2 = file.path(root, "SRR866999_sorted.bam")
f3 = file.path(root, "SRR867001_sorted.bam")
IP_BAM = c(f1,f2,f3) 

f1 = file.path(root, "SRR866998_sorted.bam")
f2 = file.path(root, "SRR867000_sorted.bam")
f3 = file.path(root, "SRR867002_sorted.bam")
INPUT_BAM = c(f1,f2,f3)

exomePeak2(bam_ip = IP_BAM,
           bam_input = INPUT_BAM,
           genome = "mm10",
           paired_end = FALSE)
```



An output folder named `exomePeak2_output` will be created in the working directory containing. The most important two files "Mod.bed" and "Mod.csv" will be used in further analysis.

```markdown
- exomePeak2_output
	- LfcGC.pdf
	- RunInfo.txt
	- Mod.bed
	- Mod.csv
	- Mod.rds
	- ADDInfo
		- ADDInfo_SizeFactors.csv
		- ADDInfo_GLM_allDesigns.csv
		- ADDInfo_ReadsCount.csv
		- ADDInfo_RPKM.csv
```



# Reference

[1] Zhen Wei (2020). exomePeak2: Bias Awared Peak Calling and Quantification for MeRIP-Seq. R package version 1.0.0. [[bioc](http://www.bioconductor.org/packages/release/bioc/html/exomePeak2.html)]

