![SKA](https://github.com/simonrharris/SKA/blob/master/images/ska.png)

## Contents
* [Introduction](#introduction)
* [Installation](#installation)
* [Temporary files](#temporary-files)
* [Usage](#usage)
* [License](#license)

## Introduction
SKA (Split Kmer Analysis) is a toolkit for prokaryotic (and small, haploid eukaryotic) DNA sequence analysis using split kmers. A split kmer is a pair of kmers in a DNA sequence that are separated by a single base. Split kmers allow rapid comparison and alignment of small genomes.

## Installation
SKA can be installed by cloning this repository and running make
```
git clone https://github.com/simonrharris/SKA
cd SKA
make
```
The executable will be compiled into a directory named bin. You can either add this bin directory to your path or move the executable into a path directory.

## Usage
```
    ska <subcommand>

    Subcommands:
    align       Reference-free alignment from split kmer files
    compare     Compare two split kmer files
    fasta       Create split kmer file from fasta file(s)
    fastq       Create split kmer file from fastq file(s)
    map         Align split kmer file(s) against a reference fasta file
    summary     Print split kmer file summary statistics
    version     Print the version and citation for ska
    weed        Weed kmers from a split kmer file
```
Please read the [SKA wiki page](https://github.com/simonrharris/SKA/wiki) for full usage instructions.
## License
SKA is free software, licensed under [GPLv3](https://github.com/simonrharris/SKA/blob/master/LICENSE).
