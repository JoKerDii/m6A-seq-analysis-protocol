# DAVID

The **D**atabase for **A**nnotation, **V**isualization and **I**ntegrated **D**iscovery (**DAVID** ) is a website providing a comprehensive set of functional annotation tools for investigators to understand biological meaning behind large list of genes [1]. We are going to use this tool to discover enriched gene functions.



## Prepare and Upload Gene Lists 

From the output file ("Mod.csv") generated from exomePeak2, we should first remove duplicated values in "geneID" column and copy all the unique IDs to a txt file ("geneID.txt"). Then upload txt files to DAVID website with Identifier as "ENTREZ_GENE_ID", species as "Homo sapiens", and "Gene List" selected. Submit all lists and wait for results.



## Analyze Results

Open "Functional Annotation Chart" and click on "Download File" to download the txt file containing results. 

![streme_motif_finding](../assets/images/M3/david_chart.png)

Import txt file into R, analyze results and plot Top 30 terms.

```R
library(readr)
library(dplyr)
chart = read_tsv("chart.txt")

generateFigure = function(chart, num, term = "GOTERM_BP_DIRECT"){

  p = selectPvalue(chart)
  frame = as.data.frame(chart %>%
    filter(Category == term) %>%
    select(c("Term", "%", "PValue")) %>%
    rename("Ratio"=`%`) %>%
    mutate(Term = as.factor(gsub("^.*?~", "",Term)), 
           Ratio = Ratio / 100))[1:30,] # make sure no less than 30 terms in total
  
  fig2 = frame %>%
    ggplot(aes(x=reorder(Term, -PValue),y=-log(PValue),fill = Ratio)) +
    geom_bar(stat="identity") +
    coord_flip() + 
    xlab("Gene Ontology: Biological Process") 
    
  return(fig2)
}

generateFigure(chart, 30, term)
```



![GO_bar_plot](../assets/images/M3/GO_barplot_mod.png)

# Reference

[1] X. Jiao, B. T. Sherman, D. W. Huang, R. Stephens, M. W. Baseler, H. C. Lane et al., "DAVID-WS: a stateful web service to facilitate gene/protein list analysis," (in eng), Bioinformatics (Oxford, England), vol. 28, no. 13, pp. 1805-1806, 2012, doi: 10.1093/bioinformatics/bts251.[[paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3381967/?report=abstract)]