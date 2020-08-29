#!/usr/bin/env python

'''Estimate RPKM of detected AMR genes 
for plots & analysis.

This script takes the following inputs:
- directory containing filtered AMRFinder tsv files
	(output from step 00)
- directory containing *_gene_RPKM.tsv files
	(https://github.com/rotheconrad/00_in-situ_GeneCoverage)

and returns the following output:
- single tsv file with length / effort normalized AMR gene abundance (RPKM)
'''

import argparse
from collections import defaultdict
import pandas as pd
import glob, os

# def generate_rela_matrix(MG_IDs, verbose):
# 	print('Generating relative abundance matrix...')

# 	matrix = defaultdict(list)
# 	unique_gene_list = []

# 	for Genes in MG_IDs.values():
# 		for gene in Genes:
# 			if gene not in unique_gene_list:
# 				unique_gene_list.append(gene)

# 	for MG in MG_IDs:
# 		if verbose:
# 			print('')
# 			print('MG ID is:', MG, '& genes are:')

# 		for gene in unique_gene_list:
# 			if gene in MG_IDs[MG]:
# 				matrix[MG].append(1)
# 				if verbose:
# 					print(gene)
# 			else:
# 				matrix[MG].append(0)

# 	binary_matrix = pd.DataFrame(matrix, index = unique_gene_list)
# 	binary_matrix.sort_index(inplace=True)
# 	binary_matrix.sort_index(axis=1, inplace=True)	

# 	if verbose:
# 		print('')
# 		print(binary_matrix)

# 	print('')
# 	print(len(unique_gene_list), 'unique genes &', len(MG_IDs),  'samples were included in matrix.')

# 	return binary_matrix

def parse_amrfinder_tsvs(input_directory, verbose):
	# generate dictionary of sequence ids & accession ids
	# to fetch sequences from NCBI entrez

	# initialize dictionary & list objects
	os.chdir(input_directory)			# change to input directory
	file_list = glob.glob("*.tsv")		# generate list of *.tsv files
	AMR_dict = defaultdict(defaultdict(list).copy)		# initialize AMR_dict

	print('Parsing AMRFinder tables...')


	for file in file_list:

		with open(file, 'r') as F:
			for line in F:
				X = line.rstrip().split('\t')
				sequence_name = X[2] 	# sequence name
				accession = X[14]		# NCBI accession for closest sequence
				if sequence_name not in AMR_dict.keys() and accession != 'NA':
					AMR_dict[sequence_name]=accession

	return AMR_dict


def main():
	# configure argparse arguments & pass to functions.
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class = argparse.RawDescriptionHelpFormatter
		)
	parser.add_argument(
		'-a', '--amrfinder_tsv_path',
		help = 'Please specify directory path containing filtered AMRFinder tsv files.',
		metavar = '',
		type=str,
		required=True
		)
	parser.add_argument(
		'-m', '--coverage_magic_path',
		help = 'Please specify directory path containing coverage magic path.',
		metavar = '',
		type=str,
		required=True
		)
	parser.add_argument(
		'-o', '--output',
		help = 'Please specify output filename & path.',
		metavar = '',
		type=str
		)
	parser.add_argument(
		'-v', '--verbose',
		help = 'Toggle volume of printed output.',
		action='store_true'
		)
	args=vars(parser.parse_args())

	AMR_dict= parse_amrfinder_tsvs(args['amrfinder_tsv_path'], args['verbose'])

	if args['verbose']:
		print('')
		print('AMR dictionary:')
		print(AMR_dict)
		# print('')
		# print('Unique gene list:')
		# print(unique_gene_list)
		# print(len(unique_gene_list), 'unique genes detected.')
		print('')
		print(len(AMR_dict),'AMR gene accessions were identified in this set.')

	print('Complete. Output written to:')
	print(args['output'])

if __name__ == "__main__":
	main()