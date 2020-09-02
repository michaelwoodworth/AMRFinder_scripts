# AMRFinder scripts
AMRFinder result analytic workflow

## Motivation
[AMRFinder plus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) is a tool that takes a set of gene sequence data and returns an output table of high confidence matches to reference sequences or hidden markov models (HMM).  This workflow documents input parameters for upstream steps to produce AMRFinder results for reproducibility as well as custom python scripts written to process output tables for downstream analysis and visualization.

## Overall workflow:
1. Prepare reads for assembly
2. Assembly
3. Predict genes
4. Annotate genes with AMRFinder
5. Estimate in-situ gene coverage with CoverageMagic workflow
6. Post-processing

	00. filter for complete genes ([00_amrfinder_filter.py](https://github.com/michaelwoodworth/AMRFinder_scripts/blob/master/01_amrfinder_binary_matrix.py))
	01. create binary presence/absence matrix ([01_amrfinder_binary_matrix.py](https://github.com/michaelwoodworth/AMRFinder_scripts/blob/master/01_amrfinder_binary_matrix.py))
	02. estimate relative abundance of identified AMR genes ([02_amrfinder_validate_and_summarize_RPKM.py](https://github.com/michaelwoodworth/AMRFinder_scripts/blob/master/02_amrfinder_validate_and_summarize_RPKM.py))

7. Analysis in R
8. Produce heatmaps in R

## Workflow detail

### 1. Prepare reads for assembly
Since I work with human-associated metagenomes, in addition to trimming with trimmomatic, I usually attempt to remove reads aligning to reference human genomes.  This can be efficiently done in one step with the Huttenhower lab tool [kneaddata](https://huttenhower.sph.harvard.edu/kneaddata/).

### 2. Assembly
I primarily assemble genomes and metagenomes with [SPAdes](https://cab.spbu.ru/software/spades/).

### 3. Predict genes
Gene prediction at the metagenome, metagenome-assembled genome (MAG), and genome level are performed using [Prodigal](https://github.com/hyattpd/Prodigal)

### 4. Annotate genes with AMRFinder
Genes predicted in step 3 are then analyzed with [AMRFinder plus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/).

### 5. Estimate gene RPKM with CoverageMagic workflow
Genes predicted in step 3 from a metagenome/genome as well as metagenome/genome reads file and filtered tabular magicblast output file are used to estimate AMRFilter-detected gene RPKM using the [CoverageMagic in situ gene coverage workflow](https://github.com/rotheconrad/00_in-situ_GeneCoverage/tree/6812ebd32c5127ce8b72ba8e520799b75f45c895).

### 6. Post-processing

00. 00_amrfinder_filter.py - this python script was written because AMRFinder produces some hits that are incomplete genes that may reduce confidence of your results.  This may be fine in an exploratory analysis, but for my purposes, I prefer to filter only hits that AMRFinder classifies as "ALLELE", "EXACT", "BLASTX", "HMM", which is the default usage.

- *add_partial_end*: 
This script allows users to also include AMRFinder hits that were partial but located at the end of a contig sequence, which could be consistent with a sequencing/assembly issue of a gene that may be complete in host.  This option is flagged with the -m add_partial_end option.

- *just_amr*:
The -j/--just_AMR flag allows filtering of just AMR results.

```console
usage: 00_amrfinder_filter.py [-h] -i  -o  [-m] [-j]

Filter AMRFinder Plus results for high confidence matches.

This script filters AMRFinder output tables for matches, with
default criteria focused on high quality & complete matches.
e.g. >90% identity, >90% match length.

Script options also allow filtering for just AMR determinants.

optional arguments:
  -h, --help      show this help message and exit
  -i , --input    Please specify AMRFinder input tsv file name & path.
  -o , --output   Please specify AMRFinder filtered filename & path.
  -m , --method   Please specify filtered AMRFinder output tsv file name &
                  path. Select from: complete -or- add_partial_end
  -j, --just_AMR  Flag to remove non-AMR AMRFinder results (e.g. virulence,
                  stress)

```

01. 01_amrfinder_binary_matrix.py - this python script searches for .tsv files in an input directory and produces a binary presence/absence matrix for all genes across all samples coded as 0 for absent and 1 as present.

```console
usage: 01_amrfinder_binary_matrix.py [-h] -i INPUT -o OUTPUT [-v]

Create summary matrix of AMRFinder Plus results for plots & analysis.

This script sumamrizes filtered tables from 00_amrfinder_filter.

optional arguments:
  -h, --help      show this help message and exit
  -i , --input    Please specify input directory path.
  -o , --output   Please specify output filename & path.
  -v, --verbose   Increase output messaging detail, print results.
```

02. 02_amrfinder_validate_and_summarize_RPKM.py - this python script performs two main tasks.

- First, a validation step is completed to confirm the presence of all AMRFinder-detected genes in the [in situ gene coverage workflow](https://github.com/rotheconrad/00_in-situ_GeneCoverage/tree/6812ebd32c5127ce8b72ba8e520799b75f45c895) ${unique_id}gene_RPKM.tsv output file.

- Second, [RPKM](https://sites.google.com/site/wiki4metagenomics/pdf/definition/rpkm-calculation) values are tallied for duplicate gene names within a sample (which is printed to STDOUT if the verbose option is selected) and then used to construct a matrix of RPKM values with rows ov unique genes by columns of metagenomes.

