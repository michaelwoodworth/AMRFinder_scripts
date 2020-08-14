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