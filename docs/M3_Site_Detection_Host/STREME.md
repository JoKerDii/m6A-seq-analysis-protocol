# STREME

[STREME](https://meme-suite.org/meme/doc/streme.html) can find motifs in large sequence datasets and report accurate significance estimates for each motif that it discovers [1]. Specifically, STREME works by discovering ungapped motifs that are enriched in the sequences or relatively enriched in them compared to the control sequences.



## Convert BED to FASTA

### 1. Install bedtools

```shell
# Download and config
$ wget https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools-2.30.0.tar.gz
$ tar -zxvf bedtools-2.30.0.tar.gz
$ cd bedtools2
$ make
$ export PATH=$PATH:/path/to/bedtools2/bin

# Verify installation
$ bedtools --version
```

### 2. Convert BED to FASTA

```shell
$ bedtools getfasta -fi /path/to/genome.fa -bed /path/to/Mod.bed -fo mod.fa
```



## Install STREME

```shell
# Download from: https://meme-suite.org/meme/doc/download.html
$ wget https://meme-suite.org/meme/meme-software/5.3.2/meme-5.3.2.tar.gz
$ tar zxf meme-5.3.2.tar.gz
$ cd meme-5.3.2
$ ./configure --prefix=/path/to/meme-5.3.2 --enable-build-libxml2 --enable-build-libxslt
$ make
$ make test
$ make install
$ export PATH=/path/to/meme-5.3.2/bin:/path/to/meme-5.3.2/libexec/meme-5.3.2:$PATH

# Verify installation
$ meme -version
$ streme -version
```



## Motif Discovery

```shell
# Usage of command-line version: https://meme-suite.org/meme/doc/streme.html?man_type=web
$ streme --dna --objfun de --o /path/to/streme_result --p Mod.fa
```



## Outputs

```markdown
 - streme_result
 	- streme.html
 	- streme.txt
 	- streme.xml
```



# Reference

[1] T. L. Bailey, "STREME: Accurate and versatile sequence motif discovery," bioRxiv, p. 2020.11.23.394619, 2020, doi: 10.1101/2020.11.23.394619. [[paper](https://www.biorxiv.org/content/10.1101/2020.11.23.394619v1.full)]