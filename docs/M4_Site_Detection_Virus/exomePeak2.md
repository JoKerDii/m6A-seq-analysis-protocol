## Peak Calling (exomePeak2)

We are going to use exomePeak2 for peak calling to find enriched m6A sites on HHV8 transcripts. The BED file as output could be used for further analysis.

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

An output folder named `exomePeak2_output` will be created in the working directory containing:

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



