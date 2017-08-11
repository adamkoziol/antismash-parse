## antismash-parse

### Introduction

Parses results from antiSMASH. Hits are often annotated vaguely. This script performs remote blastp
analyses, and associates ORFs with the original strains, and descriptions extracted from the BLAST results

### Installation

` git clone https://github.com/adamkoziol/antismash-parse.git`

` cd antismash-parse`

`python3 setup.py install`

### Requirements

- Unix-like environment (tested in Ubuntu)
- biosynML.xml file created by antiSMASH
- list of clusters/ORFs of interest. This is created by running the following command in '
your antiSMASH_run folder: `find -name "*_gene_info.txt" -exec grep "ctg*" {} \; | cut -f1 | uniq > contigs.txt`
- A folder containing all the GenBank files (.gbk) except *...final.gbk from the antiSMASH_run folder

### Arguments

- path: contains file of contigs
- -f: name of file of contigs
- -c: name and path of folder containing .gbk files

### Running

#### Example command to run the script

`parseresults.py <PATH> -f contigs.txt -c <PATH/TO/CLUSTERS`

#### Usage

```
usage: parseresults.py [-h] [-f FILE] -c CLUSTERPATH path

Parses anti-SMASH output file biosyn.xml to extract protein sequencesof genes
of interest

positional arguments:
  path                  Specify input directory

optional arguments:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  Name of file of contig numbers. This is created by
                        running the following command in your antiSMASH_run
                        folder: find -name "*_gene_info.txt" -exec grep "ctg*"
                        {} \; | cut -f1 | uniq > contigs.txt
  -c CLUSTERPATH, --clusterpath CLUSTERPATH
                        Path to .gbk cluster files created by antiSMASH.
                        Required to associate renamed contigs with original
                        files
```
