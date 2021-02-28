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

export PATH=$PATH:/home/zhendi/bedtools2/bin

# Verify installation
$ bedtools --version
```

### 2. Convert BED to FASTA

```shell
$ bedtools getfasta -fi /path/to/genome.fa -bed /path/to/Mod.bed -fo mod.fa

bedtools getfasta -s -fi /data/zhendi/wei/star-genome-ref/g2/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -bed /data/zhendi/protocol/exomePeak2/data/homo/exomePeak2_output_peakcalling_1strand/Mod.bed -split -fo mod.fa &
bedtools getfasta -s -fi /data/zhendi/wei/star-genome-ref/g2/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -bed /data/zhendi/protocol/exomePeak2/data/homo/exomePeak2_output_diff_1strand/DiffMod.bed -split -fo Diffmod1.fa &


# /data/zhendi/wei/star-genome-ref/g2/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
# /data/zhendi/protocol/exomePeak2/data/homo/exomePeak_output_peakcalling_2strand*/Mod.bed
bedtools getfasta -fi /data/zhendi/wei/star-genome-ref/g2/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -bed /data/zhendi/protocol/exomePeak2/data/homo/exomePeak2_output_peakcalling_2strand*/Mod.bed -fo mod.fa

bedtools getfasta -fi /data/zhendi/wei/star-genome-ref/g2/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa -bed /data/zhendi/protocol/exomePeak2/data/homo/exomePeak2_output_diff_2strand*/DiffMod.bed -fo Diffmod.fa
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

export PATH=/home/zhendi/meme-5.3.2/bin:/home/zhendi/meme-5.3.2/libexec/meme-5.3.2:$PATH

# Verify installation
$ meme -version
$ streme -version
```



## Motif Discovery

```shell
# Usage of command-line version: https://meme-suite.org/meme/doc/streme.html?man_type=web
$ streme --dna --objfun cd --minw 5 --maxw 10 --o /path/to/streme_result --p Mod.fa
# not sure whether --objfun de or --objfun cd, the former requires equal length of peaks
streme --dna --objfun cd --minw 5 --maxw 10 --o /data/zhendi/protocol/streme_result/mod_result4 --p mod.fa

streme --dna --objfun cd --minw 5 --maxw 10 --o /data/zhendi/protocol/streme_result/diffmod_result4 --p Diffmod.fa
```



## Outputs

```markdown
 - streme_result
 	- streme.html
 	- streme.txt
 	- streme.xml
```

Transfer the HTML file to local place by *FileZilla* (mac) or *WinSCP* (win), and open the file in browser. A screenshot of part of the HTML file is shown below.

![streme_motif_finding](../assets/images/M3/motif.png)

# Reference

[1] T. L. Bailey, "STREME: Accurate and versatile sequence motif discovery," bioRxiv, p. 2020.11.23.394619, 2020, doi: 10.1101/2020.11.23.394619. [[paper](https://www.biorxiv.org/content/10.1101/2020.11.23.394619v1.full)]