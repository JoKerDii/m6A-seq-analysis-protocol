# exomePeak2

exomePeak2 is an R/Bioconductor package which provides bias-aware quantification and peak detection for Methylated RNA immunoprecipitation sequencing data (MeRIP-Seq) [1].



## Peak Calling and Differential Modification Analysis

```R
library(exomePeak2)
set.seed(1)
root = "/path/to/data/homo"
setwd(root)

f1 = file.path(root, "SRR5978834_trimmed_s.bam")
f2 = file.path(root, "SRR5978835_trimmed_s.bam")
f3 = file.path(root, "SRR5978836_trimmed_s.bam")
IP_BAM = c(f1,f2,f3) 

f1 = file.path(root, "SRR5978827_trimmed_s.bam")
f2 = file.path(root, "SRR5978828_trimmed_s.bam")
f3 = file.path(root, "SRR5978829_trimmed_s.bam")
INPUT_BAM = c(f1,f2,f3)

f1 = file.path(root, "SRR5179446_trimmed_s.bam")
f2 = file.path(root, "SRR5179447_trimmed_s.bam")
f3 = file.path(root, "SRR5179448_trimmed_s.bam")
TREATED_INPUT_BAM = c(f1,f2,f3)

f1 = file.path(root, "SRR5978869_trimmed_s.bam")
f2 = file.path(root, "SRR5978870_trimmed_s.bam")
f3 = file.path(root, "SRR5978871_trimmed_s.bam")
TREATED_IP_BAM = c(f1,f2,f3)

exomePeak2(bam_ip = IP_BAM,
           bam_input = INPUT_BAM,
           bam_treated_input = TREATED_INPUT_BAM,
           bam_treated_ip = TREATED_IP_BAM,
           genome = "hg19",
           paired_end = FALSE)
```

An output folder named `exomePeak2_output` will be created in the working directory containing:

```markdown
- exomePeak2_output
	- LfcGC.pdf
	- RunInfo.txt
	- DiffMod.bed
	- DiffMod.csv
	- DiffMod.rds
	- ADDInfo
		- ADDInfo_SizeFactors.csv
		- ADDInfo_GLM_allDesigns.csv
		- ADDInfo_ReadsCount.csv
		- ADDInfo_RPKM.csv
```



# Reference

[1] Zhen Wei (2020). exomePeak2: Bias Awared Peak Calling and Quantification for MeRIP-Seq. R package version 1.0.0. [[bioc](http://www.bioconductor.org/packages/release/bioc/html/exomePeak2.html)]

