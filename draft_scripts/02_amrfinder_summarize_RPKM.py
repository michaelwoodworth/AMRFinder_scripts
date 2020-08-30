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
	# generate dictionary of genes & AMR sequence names

	# initialize dictionary & list objects
	os.chdir(input_directory)			# change to input directory
	file_list = glob.glob("*.tsv")		# generate list of *.tsv files
	AMR_dict = defaultdict(defaultdict(list).copy)		# initialize AMR_dict

	print('Parsing AMRFinder tables...')


	for file in file_list:

		with open(file, 'r') as F:
			for line in F:
				X = line.rstrip().split('\t')
				gene = X[0]	# gene ID / scaffold
				sequence_name = X[2] 	# sequence name
				if gene not in AMR_dict.keys():
					AMR_dict[gene]=sequence_name

	return AMR_dict


def parse_coverage_tsvs(input_directory, verbose, AMR_dict):
	# generate dictionary of prodigal genes & RPKM

	# initialize dictionary & list objects
	os.chdir(input_directory)					# change to input directory
	file_list = glob.glob("*_gene_RPKM.tsv")	# generate *.tsv list
	RPKM_dict = defaultdict(defaultdict(list).copy)		# initialize RPKM_dict

	print('Parsing CoverageMagic RPKM tables...')


	for file in file_list:

		with open(file, 'r') as F:
			for line in F:
				X = line.rstrip().split('\t')
				gene = X[0]	# prodigal gene name / scaffold
				RPKM = X[1] 	# sequence name
				if gene in AMR_dict.keys():
					RPKM_dict[gene]=RPKM

	return RPKM_dict

def validate(RPKM_dict, AMR_dict, verbose):
	# test differences in genes with RPKM values vs those without

	# initialize dictionary & list objects
	validate_dict = defaultdict(defaultdict(list).copy)		# initialize validate_dict
	valid_gene = 0
	test = 0
	RPKM_genes = RPKM_dict.keys()
	# genes_to_validate = 0

	print('Validating results...')


	for gene, sequence_name in AMR_dict.items():
		test += 1

		if gene not in RPKM_dict.keys():
			# genes_to_validate+=1
			print('      gene to validate:', gene)
			print('      sequence_name to validate', sequence_name)
			validate_dict[gene]=sequence_name
		else:
			valid_gene += 1


		# if gene in RPKM_dict.keys():
		# 	valid_gene += 1
		# else:
		# 	# genes_to_validate+=1
		# 	print('      gene to validate:', gene)
		# 	print('      sequence_name to validate', sequence_name)
		# 	validate_dict[gene]=sequence_name

	if verbose:
		print('')
		print('   Number of genes checked:', test)
		print('   Number of valid genes:', valid_gene)
		print('   Number of genes needing validation:', len(validate_dict))

	return validate_dict	

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
		help = 'Please specify output file path & prefix.',
		metavar = '',
		type=str,
		required=True
		)
	parser.add_argument(
		'-v', '--verbose',
		help = 'Toggle volume of printed output.',
		action='store_true'
		)
	args=vars(parser.parse_args())

	if args['verbose']:
		print('')
		print('AMRFinder tsv directory:', args['amrfinder_tsv_path'])
		print('Coverage Magic tsv directory:', args['coverage_magic_path'])
		print('')

	AMR_dict= parse_amrfinder_tsvs(args['amrfinder_tsv_path'], args['verbose'])
	RPKM_dict= parse_coverage_tsvs(args['coverage_magic_path'], args['verbose'], AMR_dict)
	validate_dict= validate(AMR_dict, RPKM_dict, args['verbose'])	

	print('Writing AMR_dict file...')
	with open(f"{args['output']}_AMR_dict.tsv", 'w') as file:
		for key, value in AMR_dict.items():
			newLine = f'{key}\t{value}\n'
			file.write(newLine)

	print('Writing RPKM_dict file...')
	with open(f"{args['output']}_RPKM_dict.tsv", 'w') as file:
		for key, value in RPKM_dict.items():
			newLine = f'{key}\t{value}\n'
			file.write(newLine)

	if args['verbose']:
		#print('')
		#print('AMR dictionary:')
		#print(AMR_dict)
		# print('')
		# print('RPKM dictionary:')
		# print(RPKM_dict)
		# print('')
		# print('Unique gene list:')
		# print(unique_gene_list)
		# print(len(unique_gene_list), 'unique genes detected.')
		print('')
		print(len(AMR_dict),'Genes with AMRFinder hits were identified.')
		print(len(RPKM_dict),'Genes with CoverageMagic RPKM values were identified.')

	#print('Complete. Output written to:')
	#print(args['output'])

if __name__ == "__main__":
	main()