#!/usr/bin/env python

'''Validate and summarize RPKM of AMRFinder-detected genes 
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

'''

import argparse, csv, glob, os
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
	validate_detail_dict = defaultdict(list)		# initialize validate_dict
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
			# print(f"	  gene to validate: {gene}")
			# print(f"	  sequence_name to validate: {sequence_name}")
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


				# if protein_id in validate_dict.keys():
				# 	validate_detail_dict[protein_id]=line

				if protein_id in validate_dict.keys():
					validate_detail_dict[protein_id]=[protein_id.rstrip().split('_')[0], 
					gene_symbol, sequence_name]

	return validate_dict, validate_detail_dict

def generate_RPKM_matrix(AMR_dict, RPKM_dict, MG_IDs, verbose):
	print('Generating RPKM matrix...')

	# initialize dict, lists
	# AMR_dict - 	key: gene/scaffold name 	value: sequence name (unique)
	# RPKM_dict - 	key: gene/scaffold name  	value: RPKM value

	AMR_scaffold_RPKM = defaultdict(list)
	scaf_mgid_name_rpkm = defaultdict(int)
	matrix = defaultdict(list)
	mg_genes = defaultdict(list)
	dedup_dict = defaultdict(list)

	unique_gene_list = []
	MG_IDs = []

	# generate list of unique AMR genes & MG IDs
	for scaffold, gene_name in AMR_dict.items():
		if gene_name not in unique_gene_list:
			unique_gene_list.append(gene_name)
		mg_id=scaffold.rstrip().split('_')[0]
		if mg_id not in MG_IDs:
			MG_IDs.append(mg_id)
		#print(f' gene_name: {gene_name},   MG_ID: {mg_id}')

	# output counts if verbose
	if verbose:
		print(f"{len(MG_IDs)} MG IDs | {len(unique_gene_list)} unique genes detected")

	# merge AMR_dict & RPKM_dictionaries as scaffold : gene_name, RPKM
	for d in (AMR_dict, RPKM_dict):
		for key, value in d.items():
			AMR_scaffold_RPKM[key].append(value)

	# split into structured dictionary - scaffold : gene name, RPKM,  MG ID
	for scaffold, name_RPKM in AMR_scaffold_RPKM.items():
		if len(name_RPKM) > 1:
			scaf_mgid_name_rpkm[scaffold] = { 'name' : name_RPKM[0],
											  'RPKM' : name_RPKM[1],
											  'mg_id' : scaffold.rstrip().split('_')[0]
											}
		else:
			scaf_mgid_name_rpkm[scaffold] = { 'name' : name_RPKM[0],
											  'RPKM' : 0,
											  'mg_id' : scaffold.rstrip().split('_')[0]
											}

	# tally & deduplicate AMR gene RPKM values by sample
	for MG in MG_IDs:
		dedup_list = defaultdict(list)
		dedup_count = 0

		if verbose:
			print('')
			print(f"Evaluating {MG}...")

		for scaffold, values in scaf_mgid_name_rpkm.items():
			if MG == values['mg_id'] and values['name'] not in dedup_list.keys():
				dedup_list[values['name']]={'mg_id' : values['mg_id'],
											'RPKM': values['RPKM'], 
											'count' : 1, 
											'scaffolds' : [scaffold],
											'gene_name' : values['name']}
			elif MG == values['mg_id'] and values['name'] in dedup_list.keys():
				dedup_list[values['name']]['RPKM']  = float(values['RPKM']) + float(dedup_list[values['name']]['RPKM'])
				dedup_list[values['name']]['count'] = int(dedup_list[values['name']]['count']) +1
				dedup_list[values['name']]['scaffolds'].append(scaffold)
				dedup_count += 1

		if verbose:
			#print(dedup_list)
			print(f"Number of deduplicated genes: {dedup_count}")	
			for gene_name, value in dedup_list.items():
				if value['count'] > 1:
					print('')
					print(f"Duplicate gene is: '{gene_name}' count: {dedup_list[gene_name]['count']}")
					print(f"	 tallied RPKM: {dedup_list[gene_name]['RPKM']}")
					for scaffold, o_values in scaf_mgid_name_rpkm.items():
						if o_values['name'] == gene_name and o_values['mg_id'] == MG:
							print(F"		- original name '{o_values['name']}' | original RPKM {o_values['RPKM']}")

		# add to dedup_dict to store all values across samples
		for gene_name, value in dedup_list.items():
			mg_genes[MG].append(gene_name)
			dedup_dict[f"{MG}_{gene_name}"]=[value['mg_id'], value['RPKM'], 
			value['count'], value['scaffolds'], value['gene_name']]

		#add MG gene RPKM values to matrix
		print('')
		print(f"Updating matrix with {MG} values...")
		for gene in unique_gene_list:
			if gene in mg_genes[MG]:
				matrix[MG].append(dedup_list[gene]['RPKM'])
				if verbose:
					print(f"	 {gene} (RPKM {dedup_list[gene]['RPKM']}) added to matrix")
			else:
				matrix[MG].append(0)
				if verbose:
					print(f"	 {gene} (RPKM 0, ND) added to matrix")		


	# if verbose:
	#	 print('')
	#	 print("MG gene list:")
	#	 print(f"   length: {len(mg_genes)}")
	#	 print(mg_genes)
	#	 print('')
	#	 print("dedup_dict:")
	#	 print(f"   length: {len(dedup_dict)}")
	#	 print(dedup_dict)
	#	 print('')
	#	 print("matrix:")
	#	 print(f"   length: {len(matrix)}")
	#	 print(matrix)	


	print('')
	print('Converting completed matrix to dataframe...')

	RPKM_matrix = pd.DataFrame(matrix, index = unique_gene_list)
	RPKM_matrix.sort_index(inplace=True)
	RPKM_matrix.sort_index(axis=1, inplace=True)	

	print(f"{len(MG_IDs)} MGs evaluated with {len(unique_gene_list)} unique genes.")
	
	return RPKM_matrix, unique_gene_list, dedup_dict

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
		help = 'Please specify output file path (& optional prefix).',
		metavar = '',
		type=str,
		required=True
		)
	parser.add_argument(
		'-v', '--verbose',
		help = 'Toggle volume of printed output.',
		action='store_true'
		)
	parser.add_argument(
		'-V', '--validate',
		help = 'Write genes_to_validate.tsv and deduplicated.tsv.',
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
	RPKM_matrix, unique_gene_list, dedup_dict= generate_RPKM_matrix(AMR_dict, RPKM_dict, MG_IDs, args['verbose'])

# option to write validation tsv files
	if args['validate']:
		print('')
		print(f"Writing validation files...")
		print(f"Validation output path: {args['output']}")
		def list_to_tsv(path, out_file_name, dictionary,
			col_header):
			outfile=f"{path}/{out_file_name}"
			df = pd.DataFrame.from_dict(dictionary,
				orient="index", columns=col_header)
			df.to_csv(outfile,  sep='\t')
			print(f"   ...{out_file_name} complete!")

		print('')
		print(f"   - validate_detail_dict...")

		col_header=['mg_id', 
					'gene_symbol', 
					'sequence_name']
		list_to_tsv(args['output'], 'genes_to_validate.tsv',
			validate_detail_dict, col_header)

		print('')
		print(f"   - dedup_dict...")

		col_header=['mg_id', 
			'RPKM', 'count', 'scaffolds', 'gene_name']
		list_to_tsv(args['output'], 'deduplicated_RPKM.tsv',
			dedup_dict, col_header)


		# validate_detail_dict.to_csv(f"{args['output']}/genes_to_validate.tsv", sep='\t')
		# dedup_dict.to_csv(f"{args['output']}/deduplicated.tsv", sep='\t')		

# write output tsv file
	RPKM_matrix.to_csv(f"{args['output']}/RPKM_matrix.tsv", sep='\t')

	print('')
	print('Looks like everything completed!')
	print(f"Output files written to: {args['output']}")

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