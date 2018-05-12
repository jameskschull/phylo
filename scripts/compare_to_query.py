# Given a file containing query insertions/deletions (hg38 reference), compares
# to provided query and writes results to file.

# There are two cases:
# 1. infile contains mm10 or canFam3 deletions. In this case, we search canFam3 or mm10 
# (respectively), writing to the outfile any loci that are not deleted in the other query.
# This supports the HM/HD hypothesis.

# Infile contains hg38 coordinates of deletion as well as chainID mapping to 2nd query

# 2. infile contains mm10 or canFam3 insertions. In this case, we use the hg38 coordinates
# of each locus to get the chain for canFam3/mm10 respectively. We then find 

# Infile contains hg38 coordinates of insertion as well as the chainID mapping to 2nd query

import sys
import time
import pickle
import edlib
from Bio import SeqIO

MARGIN = 5
SIMILARITY_THRESHOLD = ?

def querydel(infile, query2_chainfile, outfile):

	out = open('outfile', 'w')

	querydels = pickle.load(open(infile, 'r'))
	print "Loaded infile"

	query2_chaindict = pickle.load(open(query2_chainfile, 'r'))
	print "Loaded chains"

	# For each chain mapped to in the querydels file
	for i, key in enumerate(querydels.keys()):

		chain = query2_chaindict[key]
		chain_len = len(chain)
		ref_chain_start = chain[0][5]

		# For each querydel that maps to that chain
		for entry in querydels[key]:

			ref_insertion_start = entry[1]
			ref_insertion_end = entry[2]

			# Inspect chain if insertion is contained
			if ref_chain_start - MARGIN <= ref_insertion_start:

				ref_curr_coord = ref_chain_start

				line_n = 1

				while ref_curr_coord - MARGIN <= insertion_start_coord and line_n < chain_len:

					curr_line = chain[line_n]
					gapless_block_size = int(curr_line[0])
					ref_gap_size = int(curr_line[1])

					#if the whole insertion is contained in a gapless block, keep it
					if ref_curr_coord - MARGIN <= ref_insertion_start and curr_coord + gapless_block_size + MARGIN >= insertion_end_coord:
						print "Insertion alignment found"
						print '\t'.join(entry) + '\n'
						print 'Line in chainfile: {}'.format(curr_line)
						print 'Curr coord: {}'.format(ref_curr_coord)
						print 'End of gapless block: {}'.format(ref_curr_coord + gapless_block_size)
						out.write('\t'.join(entry) + '\n')

					ref_curr_coord += gapless_block_size + ref_gap_size
					line_n += 1

	return 


# 1. Use the query (mouse) coordinates of the insertion to get the mouse sequence
# 2. Use the reference (human) coordinates (where start and end are actually the same) to get the dog chain and find the corresponding dog coordinates
# 3. Use those coordinates to get the dog sequence
# 4. Find edit distance between dog/mouse sequence

def queryins(infile, query2_chainfile, outfile):

	out = open(outfile, 'w')

	insertions = pickle.load(open(infile, 'r'))
	query2_chains = pickle.load(open(query2_chainfile, 'r'))

	query1_whole_genome = SeqIO.to_dict(SeqIO.parse('../wholegenomes/fasta/mm10.fa', 'fasta'))
	query2_whole_genome = SeqIO.to_dict(SeqIO.parse('../wholegenomes/fasta/canFam3.fa', 'fasta'))

	print len(query1_whole_genome)
	print len(query2_whole_genome)

	# For each chain
	for i, key in enumerate(insertions.keys()):

		chain = query2_chains[key]
		if chain is None: 
			print "Chain not found!"
			continue

		q2_chr = chain[7]
		q2_chain_start = int(chain[10])
		ref_chain_start = int(chain[5])
		
		# For each insertion that maps to that chain
		for entry in insertions[key]:

			insertion_size = int(entry[0])

			# GET Q1 SEQUENCE
			q1_chr = entry[4]
			q1_start = int(entry[5])
			q1_end - int(entry[6])
			q1_seq = query1_whole_genome[q1_chr][q1_start:q1_end]

			# GET Q2 SEQUENCE
			ref_position = entry[2] # since this is an insertion, ref start and end are the same

			ref_left = ref_position - ref_chain_start # bp to ref start point
			q2_curr = q2_chain_start

			for line in chain[1:]:

				gapless_block_size = int(line[0])
				ref_block_size = int(line[1])
				query_block_size = int(line[2])

				if ref_left - gapless_block_size < 0:
					q2_start = q2_curr + ref_left
					break

				ref_left -= gapless_block_size
				q2_curr += gapless_block_size

				if ref_left - ref_block_size < 0:
					q2_start = q2_curr
					break

				ref_left -= ref_block_size
				q2_curr += query_block_size

			q2_seq = query2_whole_genome[q2_chr][q2_start:q2_start + insertion_size]

			# compare to threshold, write if they're similar enough
			if edlib.align(q1_seq, q2_seq)["editDistance"] < SIMILARITY_THRESHOLD:
				out.write('\t'.join(entry))

	return


def main():

	infile = sys.argv[1]
	query_chainfile = sys.argv[2]
	outfile = sys.argv[3]
	mode = sys.argv[4]

	if mode == 'insertion':
		queryins(infile, query_chainfile, outfile)
	elif mode == 'deletion':
		querydel(infile, query_chainfile, outfile)

	return


if __name__ == '__main__':
	main()