# Homer

Homer is a software for motif discovery and next-gen sequencing analysis. We are going to use this tool to find m6A motifs within peaks that has been found.



## Install Homer

```shell
# Download and install following the instruction: http://homer.ucsd.edu/homer/introduction/install.html
$ mkdir homer
$ cd homer
$ wget http://homer.ucsd.edu/homer/configureHomer.pl
$ perl configureHomer.pl -install

# Append to PATH environment variable
$ export PATH=$PATH:/path/to/homer/.//bin/

# Verify installation
$ findMotifs.pl
```



## Find Motifs

```shell
$ findMotifsGenome.pl Mod.bed /path/to/hg19_genome.fa /path/to/MotifOutput -rna -p 10 -len 5,6
```

Note that Homer can analyze strand-specific genomic regions for motifs by running `findMotifsGenome.pl`  with an `rna` option. Also note that filtering out the peaks longer than 1000bp can improve the reliability of motifs found by Homer.

The figure below shows the enriched motifs in peaks on hg19 transcripts.

![homer_homo_motifs](../assets/images/M3/hg19_motif_homer.png)

The figure below shows the enriched motifs in peaks on mm10 transcripts.

![homer_mm10_motifs](../assets/images/M3/mm10_motif_homer.png)



# Reference

[1] S. Heinz, C. Benner, N. Spann, E. Bertolino, Y. C. Lin, P. Laslo et al., "Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities," Molecular cell, vol. 38, no. 4, pp. 576-589, 2010. [[paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4763482/)]