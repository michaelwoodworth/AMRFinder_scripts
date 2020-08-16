#!/usr/bin/env python

'''Create summary matrix of AMRFinder Plus results for plots & analysis.

This script sumamrizes filtered tables from 00_amrfinder_filter.
'''

import argparse
from collections import defaultdict
import pandas as pd
import glob, os

def generate_binary_matrix(MG_IDs, verbose):
	print('Generating binary matrix...')

	matrix = defaultdict(list)
	unique_gene_list = []

	for Genes in MG_IDs.values():
		for gene in Genes:
			if gene not in unique_gene_list:
				unique_gene_list.append(gene)

	for MG in MG_IDs:
		if verbose:
			print('')
			print('MG ID is:', MG, '& genes are:')

		for gene in unique_gene_list:
			if gene in MG_IDs[MG]:
				matrix[MG].append(1)
				if verbose:
					print(gene)
			else:
				matrix[MG].append(0)

	binary_matrix = pd.DataFrame(matrix, index = unique_gene_list)
	binary_matrix.sort_index(inplace=True)
	binary_matrix.sort_index(axis=1, inplace=True)	

	if verbose:
		print('')
		print(binary_matrix)

	print('')
	print(len(unique_gene_list), 'unique genes &', len(MG_IDs),  'samples were included in matrix.')

	return binary_matrix

def parse_amrfinder_table(in_dir):
	#parse function that takes amrfinder table
	#	1 - pulls .tsv files from input directory &
	#	2 - creates dictionary of metagenome IDs as key,
	#		sequence name list as values.
	#returns dictionary and file list

	print('Parsing AMRFinder table...')

	MG_IDs = defaultdict(int)

	print('Pulling filenames...')

	os.chdir(in_dir)
	file_list = glob.glob("*.tsv")

	for f in file_list:
		mg_id = f.rstrip().split('.')[0]

		with open(f, 'r') as table:
	# Read AMRFinder table, populate initialized variables.
			sequences = ()

			for line in table:
				sequences = sequences + (line.rstrip().split('\t')[2],)

			MG_IDs.update({mg_id : sequences})

	return MG_IDs, file_list



def main():

	# Configure argparse
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class=argparse.RawDescriptionHelpFormatter
		)
	parser.add_argument(
		'-i', '--input',
		help='Please specify input directory path.',
		metavar = '',
		required=True
		)
	parser.add_argument(
		'-o', '--output',
		help='Please specify output filename & path.',
		metavar = '',
		required=True
		)
	parser.add_argument(
		'-v', '--verbose',
		help='Increase output messaging detail, print results.',
		required=False,
		action='store_true'
		)
	args=vars(parser.parse_args())

	print('Running script...')
	MG_IDs, file_list = parse_amrfinder_table(args['input'])
	binary_matrix = generate_binary_matrix(MG_IDs, args['verbose'])

	# write output tsv file
	binary_matrix.to_csv(args['output'], sep='\t')

	print('Complete. Output written to:')
	print(args['output'])

if __name__ == "__main__":
	main()