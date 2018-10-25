![SKA](https://github.com/simonrharris/SKA/blob/master/images/ska.png)

## Contents
* [Introduction](#introduction)
* [Installation](#installation)
* [Requirements](#requirements)
* [Usage](#usage)
* [License](#license)

## Introduction
SKA (Split Kmer Analysis) is a toolkit for prokaryotic (and any other small, haploid) DNA sequence analysis using split kmers. A split kmer is a pair of kmers in a DNA sequence that are separated by a single base. Split kmers allow rapid comparison and alignment of small genomes, and is particulalry suited for surveillance or outbreak investigation. SKA can produce split kmer files from fasta format assemblies or directly from fastq format read sequences, cluster them, align them with or without a reference sequence and provide various comparison and summary statistics. Currently all testing has been carried out on high-quality Illumina read data, so results for other platforms may vary.

## Installation
SKA can be installed by cloning this repository and running make
```
git clone https://github.com/simonrharris/SKA
```
Or by Downloading and unpacking the latest [release](https://github.com/simonrharris/SKA/releases).

Then simply navigate into the SKA directory and run make
```
cd SKA
make
```
The executable will be compiled into a directory named bin. You can either add this bin directory to your path or move the executable into a path directory.
```
sudo make install
```
will move the executable to /usr/local/bin.

## Requirements
SKA simply requires GNU make and a version of g++ which supports C++11.

## Usage
```
ska <subcommand>

Subcommands:
align		Reference-free alignment from a set of split kmer files
alleles		Create a merged split kmer file for all sequenes in one or
		more multifasta files
annotate	Locate/annotate split kmers in a reference fasta/gff file
compare		Print comparison statistics for a query split kmer file
		against a set of subject split kmer files
distance	Pairwise distance calculation and clustering from split kmer
		files
fasta		Create a split kmer file from fasta file(s)
fastq		Create a split kmer file from fastq file(s)
humanise	Print kmers from a split kmer file in human readable format
info		Print some information about one or more kmer files
map		Align split kmer file(s) against a reference fasta file
merge		Merge split kmer file(s) into one file
summary		Print split kmer file summary statistics
type		Type split kmer files using a set of allele files
unique		Output kmers unique to a set of split kmer files
version		Print the version and citation for ska
weed		Weed kmers from a split kmer file
```
Please read the [SKA wiki page](https://github.com/simonrharris/SKA/wiki) for full usage instructions.
## License
SKA is free software, licensed under [GPLv3](https://github.com/simonrharris/SKA/blob/master/LICENSE).
