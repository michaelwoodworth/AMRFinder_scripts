#!/usr/bin/env python

'''Estimate relative abundance of detected AMR genes 
for plots & analysis.

This script takes the following inputs:
- directory containing filtered AMRFinder tsv files
	(output from step 00)
- path to reads file for metagenomes

and returns the following output:
- tsv with length / depth normalized AMR gene abundance (RPKM)
'''

import argparse
from collections import defaultdict
import pandas as pd
import glob, os
import 

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

def download_AMR_sequences(AMR_dict, verbose):
	# download nucleotide fasta sequences for AMR genes
	# to fetch sequences from NCBI entrez

	# initialize dictionary & list objects
	AMR_seq_dict = defaultdict(defaultdict(list).copy)		# initialize AMR_dict

	print('Downloading protein FASTA sequences for AMR accessions...')
	print('This step can take time due to NCBI limits...')

	for accession in AMR_dict.values():


	return AMR_dict	



def main():
	# configure argparse arguments & pass to functions.
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class = argparse.RawDescriptionHelpFormatter
		)
	parser.add_argument(
		'-t', '--tsv_path',
		help = 'Please specify directory path containing filtered AMRFinder tsv files.',
		metavar = '',
		type=str,
		required=True
		)
	parser.add_argument(
		'-r', '--read_path',
		help = 'Please specify directory path containing reads.',
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
		'-i', '--percent_identity',
		help = 'Specify percent identity of AMR blast search.',
		default=95,
		required=False
		)
	parser.add_argument(
		'-l', '--length',
		help = 'Specify minimum percent read length for AMR blast search.',
		default=.7,
		required=False
		)	
	parser.add_argument(
		'-v', '--verbose',
		help = 'Toggle volume of printed output.',
		action='store_true'
		)
	args=vars(parser.parse_args())

	AMR_dict= parse_amrfinder_tsvs(args['input'], args['verbose'])
	AMR_seq_dict= download_AMR_sequences(AMR_dict, args['verbose'])

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