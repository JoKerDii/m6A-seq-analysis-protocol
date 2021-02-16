# STREME

[STREME](https://meme-suite.org/meme/doc/streme.html) can find motifs in large sequence datasets and report accurate significance estimates for each motif that it discovers [1].

## Installation

```shell
# Download software from: https://meme-suite.org/meme/doc/download.html
# Instruction: https://meme-suite.org/meme/doc/install.html?man_type=web
$ wget https://meme-suite.org/meme/meme-software/5.3.2/meme-5.3.2.tar.gz
$ tar zxf meme-5.3.2.tar.gz
$ cd meme-5.3.2
$ ./configure --prefix=$HOME/meme --enable-build-libxml2 --enable-build-libxslt
$ make
$ make test
$ make install
$ export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.3.2:$PATH

```







# Reference

[1] T. L. Bailey, "STREME: Accurate and versatile sequence motif discovery," bioRxiv, p. 2020.11.23.394619, 2020, doi: 10.1101/2020.11.23.394619. [[paper](https://www.biorxiv.org/content/10.1101/2020.11.23.394619v1.full)]