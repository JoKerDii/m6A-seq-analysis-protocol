## exomePeak2

Here we use exomePeak2 to conduct reference based quantification and differential analysis.

### Download and Convert Basic Site Information

Download single base RNA modification annotation of m6A on human genome from [m6A-Atlas database](http://10.7.6.58/m6A-Atlas/download.html), which should be tabular data in a txt file. Convert tabular data to genomic ranges by:

```R
library(readr)
library(GenomicRanges)
# Downloaded from m6A-Atlas database
my.file="/path/to/home/tangyujiao/big/download/m6A_H.sapiens_basical_information.txt"

# Load txt file
m6A_basic_info <- read_table2(my.file, col_names = FALSE)
m6A_basic_info <- m6A_basic_info[, c(1:11)]
colNames<- c("ID","chr","start","end", "num1","strand", "LOC", "ENS","RNA", "Gene", "seq")
colnames(m6A_basic_info) <- colNames

# Concert to Grange object
mod_annot <- makeGRangesFromDataFrame(m6A_basic_info,
                         keep.extra.columns=FALSE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field="chr",
                         start.field="start",
                         end.field="end",
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)

# Save Grange to rds
saveRDS(mod_annot, "/path/to/mod_annot.rds")
```



### Reference-based Analysis

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

f2 = "/path/to/mod_annot.rds"
MOD_ANNO_GRANGE <- readRDS(f2)
exomePeak2(bam_ip = IP_BAM,
           bam_input = INPUT_BAM,
           genome = "hg19",
           library_type = "1st_strand",
           paired_end = FALSE,
           mod_annot = MOD_ANNO_GRANGE)
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



