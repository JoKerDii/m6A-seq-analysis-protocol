# Obtaining GEO data Using SRA Toolkit

The Sequence Read Archive (SRA) is a publicly accessible archive for high throughput sequencing data. The [SRA Toolkit](https://ncbi.github.io/sra-tools/) from NCBI is a collection of tools for using data in the INSDC SRA. It takes the following steps to download data from SRA:



## Install and Config SRA Toolkit

```shell
# Download and extract the latest version
$ wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.9/sratoolkit.2.10.9-ubuntu64.tar.gz
$ tar -vxzf sratoolkit.2.10.9-ubuntu64.tar.gz

# Append the path to your PATH environment variable:
$ export PATH=$PATH:/path/to/sratoolkit.2.10.9-ubuntu64/bin

# Verify the installation
$ which fastq-dump
```



## Download Data from SRA

1. Access the GEO summary page by searching "GSE93676" on [GEO website](https://www.ncbi.nlm.nih.gov/geo/).

![GEOwebsite](../assets/images/M1/GEOwebsite.png)

2. Find a link for "SRA" under the heading "Relations".

![GEORelations](../assets/images/M1/GEORelations.png)

3. Click on the link (SRP096845) which sends you to a page of all the biological samples with specific runs and files in this study.

![AllRuns](../assets/images/M1/AllRuns.png)

4. To find files of interest in one comprehensive list, navigate to the bottom of the page then click: "send to" > "Run Selector" > "go". Use "Filter List" to narrow down the choices. 

![RunSelector](../assets/images/M1/RunSelector.png)

5. Extract FastQ files from SRA-accession using SRA-Toolkit

```bash
#!/bin/bash
cd /path/to/raw_data/homo/ 
for s in SRR5978827 SRR5978828 SRR5978829 SRR5978834 SRR5978835 SRR5978836 SRR5978869 SRR5978870 SRR5978871 SRR5179446 SRR5179447 SRR5179448
do 
prefetch $s
fastq-dump $s
wait
done
```

```bash
#!/bin/bash
cd /path/to/raw_data/mm10/ 
for s in SRR866997 SRR866998 SRR866999 SRR867000 SRR867001 SRR867002 SRR866991 SRR866992 SRR866993 SRR866994 SRR866995 SRR866996
do 
prefetch $s
fastq-dump $s
wait
done
```



# Reference

SRA-tools wiki on Github: https://github.com/ncbi/sra-tools/wiki

