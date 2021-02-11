# Obtaining GEO data Using SRA Toolkit

The Sequence Read Archive (SRA) is a publically accessible archive for high throughput sequencing data. The [SRA Toolkit](https://ncbi.github.io/sra-tools/) from NCBI is a collection of tools for using data in the INSDC SRA. It takes the following steps to download data from SRA:



## Install and Config SRA Toolkit

For Ubuntu

1.  Fetch the tar file from the canonical location at NCBI (December 16, 2020, version 2.10.9)

```shell
wget --output-document sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.9/sratoolkit.2.10.9-ubuntu64.tar.gz
```

2. Extract the contents of the tar file:

```shell
tar -vxzf sratoolkit.2.10.9-ubuntu64.tar.gz
```

3. Append the path to the binaries to your PATH environment variable:

```shell
export PATH=$PATH:$PWD/sratoolkit.2.10.9-ubuntu64/bin
```

4. Verify the configuration

```shell
which fastq-dump
```

This should produce the similar output as follows:

```shell
/home/User/sratoolkit.2.10.9-ubuntu64/bin/fastq-dump
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

5. Extract FASTQ file from SRA-accession using SRA-Toolkit

```shell
prefetch SRR5179431
fastq-dump SRR5179431
```

**Note 1**: you can download the text file (.txt) of the total or selected accession list by clicking on "Accession List".

**Note 2**: SRA toolkit has an option `--split-files` for downloading paired end data to two FastQ files.



# Reference

SRA-tools wiki on Github: https://github.com/ncbi/sra-tools/wiki

