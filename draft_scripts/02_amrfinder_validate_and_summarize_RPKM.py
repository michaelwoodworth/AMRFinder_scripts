#!/usr/bin/env python

'''Estimate RPKM of detected AMR genes 
for plots & analysis.

This script takes the following inputs:
- directory containing filtered AMRFinder tsv files
	(output from step 00)
- directory containing *_gene_RPKM.tsv files
	(https://github.com/rotheconrad/00_in-situ_GeneCoverage)

As an intermediate validation step, all genes input in AMRFinder 
tables are tested against all genes in the coverage_magic tsv files.
If there are any genes that are not in the submitted coverage_magic
tsv files, these are output as:
- genes_to_validate.tsv

Genes that have RPKM values are returned with following output:
- specified output file & path containting a single tsv file with 
length / effort normalized AMR gene abundance (RPKM)

'''

import argparse
from collections import defaultdict
import pandas as pd
import glob, os


def parse_amrfinder_tsvs(input_directory, verbose):
	# generate dictionary of genes & AMR sequence names

	# initialize dictionary & list objects
	os.chdir(input_directory)			# change to input directory
	file_list = glob.glob("*.tsv")		# generate list of *.tsv files
	AMR_dict = defaultdict(defaultdict(list).copy)		# initialize AMR_dict
	MG_IDs = []

	print('Parsing AMRFinder tables...')


	for file in file_list:

		mg_id = file.rstrip().split('.')[0]
		if mg_id not in MG_IDs:
			MG_IDs.append(mg_id)

		with open(file, 'r') as F:
			for line in F:
				X = line.rstrip().split('\t')
				gene = X[0]				# gene ID / scaffold
				sequence_name = X[2] 	# sequence name
				if gene not in AMR_dict.keys():
					AMR_dict[gene]=sequence_name

	return AMR_dict, MG_IDs


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
				gene = X[0]		# prodigal gene name / scaffold
				RPKM = X[1] 	# RPKM value
				if gene in AMR_dict.keys():
					RPKM_dict[gene]=RPKM

	return RPKM_dict

def validate(RPKM_dict, AMR_dict, AMR_input_directory, verbose):
	# test differences in genes with RPKM values vs those without

	# initialize dictionary & list objects
	validate_dict = defaultdict(defaultdict(list).copy)		# initialize validate_dict
	validate_detail_dict = defaultdict(defaultdict(list).copy)		# initialize validate_dict
	total_genes = 0
	valid_gene = 0
	method_hmm = 0

	mg_to_validate = []
	genes_to_validate = []
	file_list = []
	RPKM_genes = RPKM_dict.keys()

	print('Validating results...')


	for gene, sequence_name in AMR_dict.items():
		total_genes += 1

		if gene not in RPKM_dict.keys():
			# print(f"      gene to validate: {gene}")
			# print(f"      sequence_name to validate: {sequence_name}")
			validate_dict[gene]=sequence_name
		else:
			valid_gene += 1

	print('')
	print(f"   Number of genes checked: {total_genes}")
	print(f"   Number of valid genes: {valid_gene} ({round((valid_gene / total_genes * 100), 2)} %)")
	print(f"   Number of genes needing validation: {total_genes - valid_gene} ({round(((total_genes-valid_gene) / total_genes * 100), 2)} %)")
	print('')

	for gene in validate_dict.keys():
		mg_gene_id = gene.rstrip().split('_')[0]
		if mg_gene_id not in mg_to_validate:
			mg_to_validate.append(mg_gene_id)

	if verbose:
		print('')
		print('   Pulling list of metagenomes to validate...')
#		print(mg_to_validate)
		print('')

	os.chdir(AMR_input_directory)		# change to AMRFinder tsv input directory
	
	for mg in mg_to_validate:
		if glob.glob(f"{mg}*.tsv") not in file_list:
			file_list.append(glob.glob(f"{mg}*.tsv"))

	if verbose:
		print('   AMRFinder files to validate:')
		print(file_list)
		print('')

	for mg_file in file_list:
		with open(f"{AMR_input_directory}/{mg_file[0]}", 'r') as file:
			for line in file:
				X = line.rstrip().split('\t')
				protein_id = X[0]		# protein id (with contig)
				gene_symbol = X[1]		# gene symbol
				sequence_name = X[2] 	# sequence name
				scope = X[3]			# core / plus
				element_type = X[4]		# AMR, STRESS, etc.
				element_subtype = X[5]	# AMR, Metal, Acid, etc.
				amr_class = X[6]		# e.g. glycopeptide, aminoglycoside, etc.
				amr_subclass = X[7]		# e.g. vancomycin, streptomycin, etc.
				method_type = X[8]			# e.g. PARTIALP, HMM, EXACTP, etc.
				target_l = X[9]			# target length
				ref_l = X[10]			# ref length
				cov_pct = X[11]			# coverage percent [breadth]
				pct_ident = X[12]		# percent identity to reference
				align_l = X[13]			# alignment length
				accession = X[14]		# NCBI accession for closest sequence
				ref_name = X[15]		# name of closest sequence
				HMM_id = X[16]			# id of closest HMM
				HMM_desc = X[17]		# description of closest HMM

				if protein_id in validate_dict.keys():
					validate_detail_dict[protein_id]=line



	return validate_dict, validate_detail_dict

def generate_RPKM_matrix(AMR_dict, RPKM_dict, MG_IDs, verbose):
	print('Generating RPKM matrix...')

	# AMR_dict - 	key: gene/scaffold name 	value: sequence name (unique)
	# RPKM_dict - 	key: gene/scaffold name  	value: RPKM value

	matrix = defaultdict(list)
	unique_gene_list = []
	AMR_name_scaffold_RPKM = defaultdict(defaultdict(list).copy)		# initialize dict

	for gene, sequence_name in AMR_dict.items():
		if sequence_name not in unique_gene_list:
			unique_gene_list.append(sequence_name)

			AMR_name_scaffold_RPKM[sequence_name]=f"{gene}\t{RPKM_dict[gene]}"

	for mg_id in MG_IDs:
		if verbose:
			print('')
			print(f"   Evaluating {mg_id}...")

		for gene in unique_gene_list:
			

			if gene in AMR_name_scaffold_RPKM.keys():
				matrix[mg_id].append(RPKM_dict[mg_id])
				if verbose:
					print(gene)
			else:
				matrix[mg_id].append(0)

	RPKM_matrix = pd.DataFrame(matrix, index = unique_gene_list)
	RPKM_matrix.sort_index(inplace=True)
	RPKM_matrix.sort_index(axis=1, inplace=True)	

	if verbose:
		print('')
		print(RPKM_matrix)

	print('')
	print(f"{len(unique_gene_list)} unique genes & {len(MG_IDs)} samples were included in matrix.")

	return RPKM_matrix

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

	AMR_dict, MG_IDs= parse_amrfinder_tsvs(args['amrfinder_tsv_path'], args['verbose'])
	RPKM_dict= parse_coverage_tsvs(args['coverage_magic_path'], args['verbose'], AMR_dict)
	validate_dict, validate_detail_dict= validate(RPKM_dict, AMR_dict, args['amrfinder_tsv_path'], args['verbose'])	
	RPKM_matrix= generate_RPKM_matrix(AMR_dict, RPKM_dict, MG_IDs, args['verbose'])

	print('Writing validate_detail_dict file...')
	with open(f"{args['output']}/genes_to_validate.tsv", 'w') as file:
		for key, value in validate_detail_dict.items():
			newLine = f'{key}\t{value}'
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