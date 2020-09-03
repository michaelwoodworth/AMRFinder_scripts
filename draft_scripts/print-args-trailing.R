# Pheatmap

# This Rscript takes input from post-processed AMRFinder output files and
# matrices to produce heatmaps.

args <- commandArgs(trailingOnly = TRUE)

# check for installation of pacman
	# first install pacman package manager
	if(!require("pacman")) install.packages("pacman",repos = "http://cran.us.r-project.org")

	# then use pacman to check and install dependencies
	pacman::p_load(ggplot2, viridis, pheatmap)


cat(args, sep = '\n')