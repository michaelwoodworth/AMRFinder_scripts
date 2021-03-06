#!/usr/bin/env python

'''Filter AMRFinder Plus results for high confidence matches.

This script filters AMRFinder output tables for matches, with
default criteria focused on high quality & complete matches.
e.g. >90% identity, >90% match length.

Script options also allow filtering for just AMR determinants,
or conversely, only non-AMR results (e.g. virulence/stress).
'''

import argparse

def dict_writer(protein_id, sequence_name, HMM_id, 
	method_type, d, line, duplicates):

	if sequence_name in d.values():
		duplicates += 1	

	elif HMM_id in d.values():
		duplicates += 1

	else:
		d[protein_id] = line

	return d, duplicates


def method_filter(infile, outfile, method, just_AMR, no_AMR):
	d = {}				# initialize dictionary for writing matches
	d_amr = {}			# initialize dictionary for writing just AMR matches
	d_no_amr = {}		# initialize dictionary for writing just non-AMR matches
	total = 0			# counter for number of AMRFinder lines/matches
	method_match = 0	# counter for number of genes matching method criteria
	method_fail = 0		# counter for number of genes failing method criteria filter
	duplicates = 0		# counter for number of duplicate genes
	just_AMR_count = 0
	no_AMR_count = 0

	if method == "complete":
		methods = ["ALLELE", "EXACT", "BLAST", "HMM"]
	elif method == "add_partial_end":
		methods = ["ALLELE", "EXACT", "BLAST", "HMM", "PARTIAL_CONTIG_END"]

	with open(infile, 'r') as file:
		for line in file:
			total += 1
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
#			print(line)

			# print(method_type, element_type)
			if method_type.startswith(tuple(methods)):

				# add option to filter just AMR results
				if just_AMR:
					if element_type == 'AMR':
						d_amr, duplicates = dict_writer(protein_id, sequence_name,
							HMM_id, method_type, d_amr, line, duplicates)
						just_AMR_count += 1
				if no_AMR:
					if element_type != 'AMR':
						d_no_amr, duplicates = dict_writer(protein_id, sequence_name,
							HMM_id, method_type, d_no_amr, line, duplicates)
						no_AMR_count += 1
				d, duplicates = dict_writer(protein_id, sequence_name,
					 HMM_id, method_type, d, line, duplicates)

				method_match += 1

			else:
				method_fail += 1

	print('Writing output file(s)...')
	with open(f"{outfile}.tsv", 'w') as file:
		for key, value in d.items():
			file.write(value)

	if no_AMR:
		with open(f"{outfile}_noAMR.tsv", 'w') as noAMR_file:
			for key, value in d_no_amr.items():
				noAMR_file.write(value)

	if just_AMR:
		with open(f"{outfile}_justAMR.tsv", 'w') as file:
			for key, value in d_amr.items():
				file.write(value)
			
	print('Done.')
	print('')

	print('Filter statistics:')
	print(' - total hits in unfiltered AMRFinder table: ', total)
	print(' - lines passing filter:', method_match, "(", round((method_match/total * 100), 1), "% ).")
	print(' - lines failing filter:', method_fail, "(", round((method_fail/total * 100), 1), "% ).")
	if just_AMR:
		print(' - just AMR results:', just_AMR_count)
	if no_AMR:
		print(' - non-AMR results:', no_AMR_count)
	print(' - duplicate HMM ids:', duplicates, "(", round((duplicates/total * 100), 1), "% ).")


def main():
	# configure argparse arguments & pass to method_filter
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class = argparse.RawDescriptionHelpFormatter
		)
	parser.add_argument(
		'-i', '--input',
		help = 'Please specify AMRFinder input tsv file name & path.',
		metavar = '',
		type=str,
		required=True
		)
	parser.add_argument(
		'-o', '--output',
		help = 'Please specify AMRFinder filtered prefix & path for output tsv.',
		metavar = '',
		type=str,
		required=True
		)
	parser.add_argument(
		'-m', '--method',
		help = 'Please specify filtered AMRFinder output tsv file name & path.\
		Select from: complete -or- add_partial_end',
		metavar = '',
		type=str,
		default = 'complete',
		choices=['complete', 'add_partial_end']
		)
	parser.add_argument(
		'-j', '--just_AMR',
		help = 'Flag to write tsv with just AMR results',
		required=False,
		action='store_true'
		)
	parser.add_argument(
		'-v', '--virulence_stress',
		help = 'Flag to write tsv without AMR results (e.g. filter only virulence, stress)',
		required=False,
		action='store_true'
		)
	args=vars(parser.parse_args())

	print('')
	print("File:", args['input'])
	print("Method:", args['method'])

	print('Filtering results...')
	method_filter(args['input'], args['output'], 
		args['method'], args['just_AMR'], args['virulence_stress'])

if __name__ == "__main__":
	main()