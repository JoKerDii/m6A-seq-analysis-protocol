# m6AFunViewer

Funm6AViewer is a novel platform to identify, prioritize, and visualize the functional gene interaction networks driven by dynamic m6A RNA methylation unveiled from a case control study [1].

Install m6AFunViewer following the instruction on its [Github](https://github.com/NWPU-903PR/Funm6AViewer/).

```R
BiocManager::install(c("GenomicFeatures", "GenomicAlignments", "Rsamtools", "Guitar", "trackViewer", "DESeq2", "apeglm", "STRINGdb",    
                       "TxDb.Hsapiens.UCSC.hg19.knownGene", "org.Hs.eg.db"), version = "3.10")
devtools::install_github("NWPU-903PR/Funm6AViewer")
```





# Reference

[1] S.-Y. Zhang, S.-W. Zhang, X.-N. Fan, T. Zhang, J. Meng, and Y. Huang, "FunDMDeep-m6A: identification and prioritization of functional differential m6A methylation genes," Bioinformatics, vol. 35, no. 14, pp. i90-i98, 2019, doi: 10.1093/bioinformatics/btz316. [[paper](https://academic.oup.com/bioinformatics/article/35/14/i90/5529234?login=true)]