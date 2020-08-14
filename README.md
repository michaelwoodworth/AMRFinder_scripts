# AMRFinder scripts
AMRFinder result analytic workflow

## Motivation
[AMRFinder plus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) is a tool that takes a set of gene sequence data and returns an output table of high confidence matches to reference sequences or hidden markov models (HMM).  This workflow documents input parameters for upstream steps to produce AMRFinder results for reproducibility as well as custom python scripts written to process output tables for downstream analysis and visualization.

## Overall workflow:
1. Prepare reads for assembly
2. Assembly
3. Predict genes
4. Annotate genes with AMRFinder
5. Post-processing
	1. filter for complete genes ([00_amrfinder_filter.py](https://github.com/michaelwoodworth/AMRFinder_scripts/blob/master/01_amrfinder_binary_matrix.py))
	2. create binary presence/absence matrix ([01_amrfinder_binary_matrix.py](https://github.com/michaelwoodworth/AMRFinder_scripts/blob/master/01_amrfinder_binary_matrix.py))
	3. estimate relative abundance of identified AMR genes
6. Analysis in R
7. Produce heatmaps in R

***

## Workflow detail

### 1. Prepare reads for assembly
Since I work with human-associated metagenomes, in addition to trimming with trimmomatic, I usually attempt to remove reads aligning to reference human genomes.  This can be efficiently done in one step with the Huttenhower lab tool [kneaddata](https://huttenhower.sph.harvard.edu/kneaddata/).

### 2. Assembly
I primarily assemble genomes and metagenomes with [SPAdes](https://cab.spbu.ru/software/spades/).

### 3. Predict genes
Gene prediction at the metagenome, metagenome-assembled genome (MAG), and genome level are performed using [Prodigal](https://github.com/hyattpd/Prodigal)

### 4. Annotate genes with AMRFinder
Genes predicted in step 3 are then analyzed with [AMRFinder plus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/).

### 5. Post-processing

1. 00_amrfinder_filter.py - this python script was written because AMRFinder produces some hits that are incomplete genes that may reduce confidence of your results.  This may be fine in an exploratory analysis, but for my purposes, I prefer to filter only hits that AMRFinder classifies as "ALLELE", "EXACT", "BLASTX", "HMM", which is the default usage.

This script allows users to also include AMRFinder hits that were partial but located at the end of a contig sequence, which could be consistent with a sequencing/assembly issue of a gene that may be complete in host.  This option is flagged with the -m add_partial_end option.

```console
usage: 00_amrfinder_filter.py [-h] -i  -o  [-m]

Filter AMRFinder Plus results for high confidence matches.

This script filters AMRFinder output tables for matches, with
default criteria focused on high quality & complete matches.
e.g. >90% identity, >90% match length.

Parameter options allow user to vary match length, identity,
and inclusion of partial matches at the end of a contig.

optional arguments:
  -h, --help        show this help message and exit
  -i , --in_file    Please specify AMRFinder output tsv file name & path.
  -o , --out_file   Please specify AMRFinder filtered filename & path.
  -m , --method     Please specify filtered AMRFinder output tsv file name &
                    path. Select from: complete -or- add_partial_end
```

2. 01_amrfinder_binary_matrix.py - this python script searches for .tsv files in an input directory and produces a binary presence/absence matrix for all genes across all samples coded as 0 for absent and 1 as present.

```console
usage: 01_amrfinder_binary_matrix.py [-h] -i INPUT -o OUTPUT [-v]

Create summary matrix of AMRFinder Plus results for plots & analysis.

This script sumamrizes filtered tables from 00_amrfinder_filter.

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Please specify input directory path.
  -o OUTPUT, --output OUTPUT
                        Please specify output filename & path.
  -v, --verbose         Increase output messaging detail, print results.
```

3. 02_amrfinder_estimate_gene_RPKM.py - this python script normalizes AMRFinder gene hit relative abundance across sampes using reads per kilobase per million mapped reads [RPKM](https://sites.google.com/site/wiki4metagenomics/pdf/definition/rpkm-calculation).

Where:
RPKM = numReads / (geneLength/1000 * totalNumReads/1,000,000)

