## m6A-seq Data Analysis (exomePeak2)

Here we are going to use exomePeak2 for peak calling and differential methylation analysis. The BED file as output will be used for further analysis.

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

f1 = file.path(root, "SRR5179446_sorted.bam")
f2 = file.path(root, "SRR5179447_sorted.bam")
f3 = file.path(root, "SRR5179448_sorted.bam")
TREATED_INPUT_BAM = c(f1,f2,f3)

f1 = file.path(root, "SRR5978869_sorted.bam")
f2 = file.path(root, "SRR5978870_sorted.bam")
f3 = file.path(root, "SRR5978871_sorted.bam")
TREATED_IP_BAM = c(f1,f2,f3)

exomePeak2(bam_ip = IP_BAM,
           bam_input = INPUT_BAM,
           bam_treated_input = TREATED_INPUT_BAM,
           bam_treated_ip = TREATED_IP_BAM,
           genome = "hg19",
           library_type = "1st_strand",
           paired_end = FALSE)
```

```bash
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

f1 = file.path(root, "SRR866992_sorted.bam")
f2 = file.path(root, "SRR866994_sorted.bam")
f3 = file.path(root, "SRR866996_sorted.bam")
TREATED_INPUT_BAM = c(f1,f2,f3)

f1 = file.path(root, "SRR866991_sorted.bam")
f2 = file.path(root, "SRR866993_sorted.bam")
f3 = file.path(root, "SRR866995_sorted.bam")
TREATED_IP_BAM = c(f1,f2,f3)

exomePeak2(bam_ip = IP_BAM,
           bam_input = INPUT_BAM,
           bam_treated_input = TREATED_INPUT_BAM,
           bam_treated_ip = TREATED_IP_BAM,
           genome = "mm10",
           library_type = "1st_strand",
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

