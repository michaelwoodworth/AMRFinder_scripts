#!/usr/bin/env python

'''Filter AMRFinder Plus results for high confidence matches.

This script filters AMRFinder output tables for matches, with
default criteria focused on high quality & complete matches.
e.g. >90% identity, >90% match length.

Parameter options allow user to vary match length, identity,
and inclusion of partial matches at the end of a contig.
'''

import argparse

def dict_writer(protein_id, HMM_id, method_type, d, line, duplicates):

	if HMM_id in d:
		duplicates += 1

		for m_type in methods:
			if method_type.startswith(mtype):
				d[protein_id] = line
				break

	else:
		d[protein_id] = line

	#print(d)
	return d, duplicates


def method_filter(infile, outfile, method):
	d = {}				# initialize dictionary for writing matches
	total = 0			# counter for number of AMRFinder lines/matches
	method_match = 0	# counter for number of genes matching method criteria
	method_fail = 0		# counter for number of genes failing method criteria filter
	duplicates = 0		# counter for number of duplicate genes

	print("file:", infile)
	print("method:", method)

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

			if method_type.startswith(tuple(methods)):
				#print(line)

				d, duplicates = dict_writer(protein_id, HMM_id, method_type, d, line, duplicates)
				method_match += 1

			else:
				method_fail += 1

	with open(outfile, 'w') as file:
		for key, value in d.items():
			file.write(value)


	print('Number of total hits in unfiltered AMRFinder table: ', total)
	print('Number of lines passing filter:', method_match,"(",(method_match/total),"% ).")
	print('Number of lines failing filter:', method_fail,"(",(method_fail/total),"% ).")
	print('Number of duplicate HMM ids:', duplicates,"(",(duplicates/total),"% ).")


def main():
	# configure argparse arguments & pass to method_filter
	parser = argparse.ArgumentParser(
		description=__doc__,
		formatter_class = argparse.RawDescriptionHelpFormatter
		)
	parser.add_argument(
		'-i', '--in_file',
		help = 'Please specify AMRFinder output tsv file name & path.',
		metavar = '',
		type=str,
		required=True
		)
	parser.add_argument(
		'-o', '--out_file',
		help = 'Please specify AMRFinder filtered filename & path.',
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
		choices=['complete', 'add_partial_end']
		)
	args=vars(parser.parse_args())

	print('Filtering results...')
	method_filter(args['in_file'], args['out_file'], args['method'])
	print('\n')

if __name__ == "__main__":
	main()