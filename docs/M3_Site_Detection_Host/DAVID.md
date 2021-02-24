# DAVID



## Prepare and Upload Gene Lists 

From the output file ("DiffMod.csv") generated from exomePeak2, we remove duplicated values in "geneID" column and copy all the unique IDs to a txt file ("geneID.txt"). We select all methylation sites as background data, and genes with Log2(FC) > 0 and adjusted p-value < 0.05 as foreground data. 

Upload txt files to DAVID website with Identifier as "ENTREZ_GENE_ID", species as "Homo sapiens", and "Gene List" selected. Submit all lists and wait for results.



## Analyze Results

Open "Functional Annotation Chart" and click on "Download File" to download the txt file containing results. 

![streme_motif_finding](../assets/images/M3/david_chart.png)