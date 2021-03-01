# m6AFunViewer

Funm6AViewer is a novel platform to identify, prioritize, and visualize the functional gene interaction networks driven by dynamic m6A RNA methylation unveiled from a case control study [1].

Install m6AFunViewer following the instruction on its [Github](https://github.com/NWPU-903PR/Funm6AViewer/).

```R
BiocManager::install(c("GenomicFeatures", "GenomicAlignments", "Rsamtools", "Guitar", "trackViewer", "DESeq2", "apeglm", "STRINGdb",    
                       "TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"), version = "3.10")
devtools::install_github("NWPU-903PR/Funm6AViewer")
```

```R
library(Funm6AViewer)

setwd("/path/to/homo_result/")
ip_bams <- c("SRR5978869_sorted.bam","SRR5978870_sorted.bam", "SRR5978871_sorted.bam", "SRR5978834_sorted.bam", "SRR5978835_sorted.bam", "SRR5978836_sorted.bam")
input_bams <- c("SRR5179446_sorted.bam","SRR5179447_sorted.bam", "SRR5179448_sorted.bam", "SRR5978827_sorted.bam", "SRR5978828_sorted.bam", "SRR5978829_sorted.bam")
sample_condition <- c("treated", "treated", "treated", "untreated", "untreated", "untreated")

grlist <- makegrreadsfrombam(IP_bams = ip_bams,    
                             Input_bams = input_bams,    
                             condition = sample_condition,    
                             minimal_alignment_MAPQ = 30,    
                             txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,    
                             savepath = "/path/to/Funm6A/")    
deinfo <- getdeinfo(grlist = grlist,    
                    txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,    
                    savepath = savepath)

datapath <- "/path/to/Funm6AViewer_data" # downloaded required data
enrich_input_directory <- "/path/to/Funm6AViewer_data"

setwd("/path/to/Funm6A/")
# Saved from outputs of exomePeak2
dminfo <- "dminfo.xlsx"
dminfo <- read.table(dminfo, header = TRUE, stringsAsFactors = FALSE)
deinfo <- "deinfo.xlsx"
deinfo <- read.table(deinfo, header = TRUE, stringsAsFactors = FALSE)

savepath <- "/path/to/Funm6A/"
siggene <- c("CCNT1", "MYC", "BCL2")
permutime <- 1000
re <- funm6aviewer(dminfo, deinfo, grlist, intrested_gene =  siggene, permutime = permutime, version = "10", datapath = datapath, enrich_input_directory = enrich_input_directory, savepath = savepath)
```



# Reference

[1] S.-Y. Zhang, S.-W. Zhang, X.-N. Fan, T. Zhang, J. Meng, and Y. Huang, "FunDMDeep-m6A: identification and prioritization of functional differential m6A methylation genes," Bioinformatics, vol. 35, no. 14, pp. i90-i98, 2019, doi: 10.1093/bioinformatics/btz316. [[paper](https://academic.oup.com/bioinformatics/article/35/14/i90/5529234?login=true)] 