# AMRFinder scripts
AMRFinder result analytic workflow

## Motivation
[AMRFinder plus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/) is a tool that takes a set of gene sequence data and returns an output table of high confidence matches to reference sequences or hidden markov models (HMM) for antimicrobial resistance and virulence determinants.  This workflow documents input parameters for upstream steps to produce AMRFinder results for reproducibility as well as custom python scripts written to process output tables for downstream analysis and visualization.

## Overall workflow:
1. Prepare reads for assembly
2. Assembly
3. Predict genes
4. Annotate genes with AMRFinder
5. Estimate in-situ gene coverage with CoverageMagic workflow
6. Post-processing

	1. filter for complete genes ([00_amrfinder_filter.py](https://github.com/michaelwoodworth/AMRFinder_scripts/blob/master/01_amrfinder_binary_matrix.py))
	2. create binary presence/absence matrix ([01_amrfinder_binary_matrix.py](https://github.com/michaelwoodworth/AMRFinder_scripts/blob/master/01_amrfinder_binary_matrix.py))
	3. estimate relative abundance of identified AMR genes ([02_amrfinder_validate_and_summarize_RPKM.py](https://github.com/michaelwoodworth/AMRFinder_scripts/blob/master/02_amrfinder_validate_and_summarize_RPKM.py))

7. Analysis in R
8. Produce heatmaps in R

## Workflow detail

### 1. Prepare reads for assembly
Since I work with human-associated metagenomes, in addition to trimming with trimmomatic, I usually attempt to remove reads aligning to reference human genomes.  This can be efficiently done in one step with the Huttenhower lab tool [kneaddata](https://huttenhower.sph.harvard.edu/kneaddata/).

```console
# test
kneaddata --input $R1 --input $R2 --run-bmtagger --remove-intermediate-output -db $db_path --output $outdir

# for loop
for ID in `cat acc_list.txt`; do R1=${indir}/${ID}_1.fastq; R2=${indir}/${ID}_2.fastq; kneaddata --input $R1 --input $R2 --run-bmtagger --remove-intermediate-output -t 10 --reference-db $db_path -o $outdir; echo $ID complete.; done

```

### 2. Assembly
I primarily assemble genomes and metagenomes with [SPAdes](https://cab.spbu.ru/software/spades/).

### 3. Predict genes
Gene prediction at the metagenome, metagenome-assembled genome (MAG), and genome level are performed using [Prodigal](https://github.com/hyattpd/Prodigal) or PROKKA.

```console
########## Prodigal example

# test
prodigal -a ${ID}.faa -d ${ID}.fna -f gff -i ${scaffold} -o ${ID}.gff -p meta

# for loop
for ID in `cat acc_list.txt`; do scaffold=${indir}/${ID}/scaffolds.fasta; prodigal -a ${ID}.faa -d ${ID}.fna -f gff -i ${scaffold} -o ${ID}.gff -p meta; echo $ID complete.; done

```

### 4. Annotate genes with AMRFinder
Genes predicted in step 3 are then analyzed with [AMRFinder plus](https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/).

### 5. Estimate gene RPKM with CoverageMagic workflow
Genes predicted in step 3 from a metagenome/genome as well as metagenome/genome reads file and filtered tabular magicblast output file are used to estimate AMRFilter-detected gene RPKM using the [CoverageMagic in situ gene coverage workflow](https://github.com/rotheconrad/00_in-situ_GeneCoverage/tree/6812ebd32c5127ce8b72ba8e520799b75f45c895).

### 6. Post-processing

**1. 00_amrfinder_filter.py** - this python script was written because AMRFinder produces some hits that are incomplete genes that may reduce confidence of your results.  This may be fine in an exploratory analysis, but for my purposes, I prefer to filter only hits that AMRFinder classifies as "ALLELE", "EXACT", "BLASTX", "HMM", which is the default usage.

- *add_partial_end*: 
This script allows users to also include AMRFinder hits that were partial but located at the end of a contig sequence, which could be consistent with a sequencing/assembly issue of a gene that may be complete in host.  This option is flagged with the -m add_partial_end option.

- *just_amr*:
The -j/--just_AMR flag filters and writes a tsv file with just AMR results.

- *virulence_stress*:
The -v/--virulence_stress flag filters and writes a tsv file with non-AMR results.

```console
usage: 00_amrfinder_filter.py [-h] -i  -o  [-m] [-j] [-v]

Filter AMRFinder Plus results for high confidence matches.

This script filters AMRFinder output tables for matches, with
default criteria focused on high quality & complete matches.
e.g. >90% identity, >90% match length.

Script options also allow filtering for just AMR determinants,
or conversely, only non-AMR results (e.g. virulence/stress).

optional arguments:
  -h, --help            show this help message and exit
  -i , --input          Please specify AMRFinder input tsv file name & path.
  -o , --output         Please specify AMRFinder filtered prefix & path for
                        output tsv.
  -m , --method         Please specify filtered AMRFinder output tsv file name
                        & path. Select from: complete -or- add_partial_end
  -j, --just_AMR        Flag to write tsv with just AMR results
  -v, --virulence_stress
                        Flag to write tsv without AMR results (e.g. filter
                        only virulence, stress)

```

**2. 01_amrfinder_binary_matrix.py** - this python script searches for .tsv files in an input directory and produces a binary presence/absence matrix for all genes across all samples coded as 0 for absent and 1 as present.

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

**3. 02_amrfinder_validate_and_summarize_RPKM.py** - this python script performs two main tasks.

- First, a validation step is completed to confirm the presence of all AMRFinder-detected genes in the [in situ gene coverage workflow](https://github.com/rotheconrad/00_in-situ_GeneCoverage/tree/6812ebd32c5127ce8b72ba8e520799b75f45c895) ${unique_id}gene_RPKM.tsv output file.

- Second, [RPKM](https://sites.google.com/site/wiki4metagenomics/pdf/definition/rpkm-calculation) values are tallied for duplicate gene names within a sample (which is printed to STDOUT if the verbose option is selected) and then used to construct a matrix of RPKM values with rows ov unique genes by columns of metagenomes.

```console
usage: 02_amrfinder_validate_and_summarize_RPKM.py [-h] -a  -m  -o  [-v] [-V]

Validate and summarize RPKM of AMRFinder-detected genes 
for plots & analysis.

This script takes the following inputs:
- directory containing filtered AMRFinder tsv files
  (output from step 00)
- directory containing ${uniqueID}_gene_RPKM.tsv files
  (https://github.com/rotheconrad/00_in-situ_GeneCoverage)

With intermediate validation steps (option -V):

- all genes input in AMRFinder tables are tested against all genes 
in the coverage_magic tsv files. If there are any genes that are 
not in the submitted coverage_magic tsv files, these are optionally 
output as: genes_to_validate.tsv

- all duplicated gene RPKM values are summed by sample.  
Input contigs/scaffolds hosting the detected genes are listed with
the summed (deduplicated) RPKM values and gene sequence name and
output as: deduplicated_RPKM.tsv

Genes that have RPKM values are returned with following output:
- specified output file & path containting a single tsv file with 
length / effort normalized AMR gene abundance (RPKM)

optional arguments:
  -h, --help            show this help message and exit
  -a , --amrfinder_tsv_path 
                        Please specify directory path containing filtered
                        AMRFinder tsv files.
  -m , --coverage_magic_path 
                        Please specify directory path containing coverage
                        magic path.
  -o , --output         Please specify output file path (& optional prefix).
  -v, --verbose         Toggle volume of printed output.
  -V, --validate        Write genes_to_validate.tsv and deduplicated.tsv.
```