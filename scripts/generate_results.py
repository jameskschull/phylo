## EDLIB METHOD ##

## Takes indel locations and uses reference-query chains to generate appropriate results.
## Call using appropriate args for the case at hand.

## Mode args: 

# di: reference deletion, search for insertion in query (e.g. hg38-.mm10+.canFam3+.txt)
# dd: reference deletion, search for deletion in query (e.g. hg38-.mm10+.canFam3+.loxAfr3-.txt)
# ii: reference insertion, search for insertion in query (e.g. hg38+.mm10-.canFam3+.txt)
# id: reference insertion, search for deletion in query (e.g. hg38+.mm10-.canFam3+.loxAfr3-txt)

import sys
import time
import pickle
import edlib
from Bio import SeqIO
from collections import defaultdict

INS_THRESHOLD = 0.7
MARGIN = 5
L_MARGIN = 10

def ref_deletion_find_insertion(sites_dict, chain_dict, outfile):

	out = open(outfile, 'w')
	print "Writing evidence to {}".format(outfile)

	query1_whole_genome = SeqIO.to_dict(SeqIO.parse('/cluster/u/jschull/phylo/wholegenomes/fasta/mm10.fa', 'fasta'))
	print "Loaded query1 genome."
	query2_whole_genome = SeqIO.to_dict(SeqIO.parse('/cluster/u/jschull/phylo/wholegenomes/fasta/canFam3.fa', 'fasta'))
	print "Loaded query2 genome."

	# For each chain
	for chain_id in sites_dict.keys():

		chain = chain_dict.get(chain_id, None)
		if chain is None: 
			print "Chain not found!"
			continue

		# For each insertion that maps to that chain
		for site in sites_dict[chain_id]:

			ref_curr_coord = int(chain[0][5]) # ref chain start
			ref_ins_position = int(site[3]) # since this is an insertion, ref start and end are the same
			insertion_size = int(site[1])

			### GET Q1 SEQUENCE ###
			q1_chr, q1_strand, q1_start, q1_end = site[5], site[6], int(site[7]), int(site[8])
			q1_chrom_size = mm_chrom_sizes[q1_chr]

			#######################

			### GET Q2 SEQUENCE ###
			
			q2_curr_coord = int(chain[0][10]) # q2 chain start
			q2_chr, q2_chrom_size, q2_strand, q2_start, q2_end = chain[0][7], int(chain[0][8]), chain[0][9], None, None 

			for line in chain[1:]:

				if len(line) < 3: continue

				gapless_block_size = int(line[0])
				query_gap_size = int(line[1])
				ref_gap_size = int(line[2])

				# Case 1: still before ref insertion start
				if ref_curr_coord + gapless_block_size < ref_ins_position - MARGIN:
					
					ref_curr_coord += gapless_block_size
					q2_curr_coord += gapless_block_size

				# Case 2: viable insertion found
				elif (ref_ins_position - MARGIN <= ref_curr_coord + gapless_block_size <= ref_ins_position + MARGIN
					 and query_gap_size == 0 and ref_gap_size > 0):

					q2_start = q2_curr_coord + gapless_block_size
					q2_end = q2_start + ref_gap_size
					break

				# Case 3: move beyond insertion start
				elif ref_curr_coord + gapless_block_size > ref_ins_position + MARGIN:
					break

				# Account for gaps
				ref_curr_coord += query_gap_size
				q2_curr_coord += ref_gap_size

			#######################

			# Viable insertion wasn't found
			if q2_start is None: continue
			
			q1_forward_start, q1_forward_end, q2_forward_start, q2_forward_end = None, None, None, None

			# Account for strand 
			if q1_strand == '-':
				q1_forward_start = q1_chrom_size - q1_end
				q1_forward_end = q1_chrom_size - q1_start
				q1_seq = query1_whole_genome[q1_chr][q1_forward_start:q1_forward_end].reverse_complement()
			elif q1_strand == '+':
				q1_seq = query1_whole_genome[q1_chr][q1_start:q1_end]

			if q2_strand == '-':
				q2_forward_start = q2_chrom_size - q2_end
				q2_forward_end = q2_chrom_size - q2_start
				q2_seq = query2_whole_genome[q2_chr][q2_forward_start:q2_forward_end].reverse_complement()
			elif q2_strand == '+':
				q2_seq = query2_whole_genome[q2_chr][q2_start:q2_end]

			q1_seq, q2_seq = str(q1_seq.seq), str(q2_seq.seq)

			# Calculate similarity
			len_longer_seq = max(len(q1_seq), len(q2_seq))

			similarity = (len_longer_seq - int(edlib.align(q1_seq, q2_seq)["editDistance"]))/float(len_longer_seq)

			# Compare to threshold, write if they're similar enough
			if similarity > INS_THRESHOLD:
				out.write('\t'.join(site) + '\t' + 'similarity: {}'.format(round(similarity, 2)) + '\t' + 'canStrand: {}'.format(q2_strand) + '\t' + 'mouseBEDcoords: {}-{}'.format(q1_forward_start, q1_forward_end) + '\t' + 'dogBEDcoords: {}-{}'.format(q2_forward_start, q2_forward_end) + '\n')

	return

def ref_deletion_find_deletion(sites_dict, chain_dict, outfile):

	print "Reference deletion, finding deletions in query2."
	print 'Margin: {}'.format(L_MARGIN)
	out = open(outfile, 'w')
	results = set([])

	# For each chain
	for chain_id in sites_dict.keys():

		print chain_id 

		chain = chain_dict.get(chain_id, None)
		if chain is None: 
			print "Chain not found!"
			continue

		# For each insertion that maps to that chain
		for site in sites_dict[chain_id]:

			print "Gene ID: {}".format(site[0])

			ref_curr_coord = int(chain[0][5]) # ref chain start
			ref_ins_position = int(site[3])

			print "Insertion position: {}".format(ref_ins_position)

			#### SEARCH CHAIN ####

			for line in chain[1:]:

				if len(line) < 3: 
					print "Reached end of chain"
					break

				gapless_block_size = int(line[0])
				query_gap_size = int(line[1])
				ref_gap_size = int(line[2])

				# print "Current coordinate: {}. Insertion coordinate: {}. Gapless size: {}. Query gap size: {}".format(ref_curr_coord, ref_ins_position, gapless_block_size, query_gap_size)

				# Case 1: still before ref insertion start
				if ref_curr_coord + gapless_block_size < ref_ins_position - L_MARGIN:	
					ref_curr_coord += gapless_block_size

				# Case 2: deletion found
				elif ref_curr_coord <= ref_ins_position - L_MARGIN and ref_curr_coord + gapless_block_size >= ref_ins_position + L_MARGIN and ref_gap_size == 0:
					print "Found a deletion!"
					print "GAPLESS: {}, RCC: {}, RCC + GAPLESS: {}".format(gapless_block_size, ref_curr_coord, ref_curr_coord + gapless_block_size)
					results.add('\t'.join(site) + '\n')
					break

				# Case 3: move beyond insertion range
				elif ref_curr_coord > ref_ins_position + L_MARGIN:
					print "No deletion found..."
					break

				# Account for gaps
				ref_curr_coord += query_gap_size

	for result in results:
		out.write(result)

	return


def ref_insertion_find_insertion(sites_dict, chain_dict, outfile):

	out = open(outfile, 'w')
	print "Writing evidence to {}".format(outfile)

	ref_whole_genome = SeqIO.to_dict(SeqIO.parse('/cluster/u/jschull/phylo/wholegenomes/fasta/hg38.fa', 'fasta'))
	print "Loaded reference genome."
	query2_whole_genome = SeqIO.to_dict(SeqIO.parse('/cluster/u/jschull/phylo/wholegenomes/fasta/mm10.fa', 'fasta'))
	print "Loaded query2 genome."

	# For each chain
	for chain_id in sites_dict.keys():

		chain = chain_dict.get(chain_id, None)
		if chain is None: 
			print "Chain not found!"
			continue

		# For each insertion that maps to that chain
		for site in sites_dict[chain_id]:

			ref_curr_coord = int(chain[0][5]) # ref chain start
			ref_chr, ref_ins_start, ref_ins_end = site[2], int(site[3]), int(site[4]) # since this is an insertion, ref start and end are the same
			insertion_size = int(site[1])

			ref_seq = ref_whole_genome[ref_chr][ref_ins_start:ref_ins_end]

			### GET Q2 SEQUENCE ###
			
			q2_curr_coord = int(chain[0][10]) # q2 chain start
			q2_chr, q2_chrom_size, q2_strand, q2_start, q2_end = chain[0][7], int(chain[0][8]), chain[0][9], None, None 

			### FIND Q2 COORDINATES ###
			for line in chain[1:]:

				if len(line) < 3:
					print "Reached end of chain, found nothing."
					continue

				gapless_block_size = int(line[0])
				query_gap_size = int(line[1])
				ref_gap_size = int(line[2])

				# print 'RCC: {}. Insertion start: {}. Gapless size: {}.'.format(ref_curr_coord, ref_ins_start, gapless_block_size)
				## Deal with gapless block ##

				# Insertion start contained in gapless block
				if ref_curr_coord - MARGIN <= ref_ins_start and ref_curr_coord + gapless_block_size + MARGIN >= ref_ins_start:

					q2_start = q2_curr_coord + (ref_curr_coord + gapless_block_size - ref_ins_start)
					q2_end = q2_start + insertion_size
					break

				ref_curr_coord += gapless_block_size
				q2_curr_coord += gapless_block_size

				## Deal with gaps ##

				# print 'RCC: {}. Insertion start: {}. Query gap size: {}.'.format(ref_curr_coord, ref_ins_start, query_gap_size)

				# Insertion start contained in query gap
				if ref_curr_coord  - MARGIN <= ref_ins_start and ref_curr_coord + query_gap_size + MARGIN >= ref_ins_start:

					q2_start = q2_curr_coord 
					q2_end = q2_start + insertion_size
					break

				ref_curr_coord += query_gap_size
				q2_curr_coord += ref_gap_size

			#######################
			
			# Didn't find insertion range
			if q2_start == None or q2_end == None: continue

			q1_forward_start, q1_forward_end, q2_forward_start, q2_forward_end = None, None, None, None

			# Account for strand 
			if q2_strand == '-':
				q2_forward_start = q2_chrom_size - q2_end
				q2_forward_end = q2_chrom_size - q2_start
				q2_seq = query2_whole_genome[q2_chr][q2_forward_start:q2_forward_end].reverse_complement()
			elif q2_strand == '+':
				q2_seq = query2_whole_genome[q2_chr][q2_start:q2_end]

			ref_seq, q2_seq = str(ref_seq.seq), str(q2_seq.seq)

			# Calculate similarity
			len_longer_seq = max(len(ref_seq), len(q2_seq))

			similarity = (len_longer_seq - int(edlib.align(ref_seq, q2_seq)["editDistance"]))/float(len_longer_seq)

			# Compare to threshold, write if they're similar enough
			if similarity > INS_THRESHOLD:
				out.write('\t'.join(site) + '\t' + 'similarity: {}'.format(round(similarity, 2)) + '\n')

	return

def ref_insertion_find_deletion(sites_dict, chain_dict, outfile):

	out = open(outfile, 'w')
	results = set([])

	# For each chain
	for chain_id in sites_dict.keys():

		chain = chain_dict.get(chain_id, None)
		if chain is None: 
			print "Chain not found!"
			continue

		# For each insertion that maps to that chain
		for site in sites_dict[chain_id]:

			ref_curr_coord = int(chain[0][5]) # ref chain start
			ref_ins_start, ref_ins_end = int(site[3]), int(site[4])

			#### SEARCH CHAIN ####

			for line in chain[1:]:

				if len(line) < 3: continue

				gapless_block_size = int(line[0])
				query_gap_size = int(line[1])
				ref_gap_size = int(line[2])

				# Case 1: still before ref insertion start
				if ref_curr_coord + gapless_block_size < ref_ins_start - MARGIN:	
					
					ref_curr_coord += gapless_block_size

				# Case 2: deletion found
				elif (ref_curr_coord + gapless_block_size >= ref_ins_start - MARGIN
				      and ref_curr_coord + gapless_block_size <= ref_ins_start + MARGIN
				      and ref_gap_size == 0
				      and query_gap_size > 0
				      and ref_curr_coord + gapless_block_size + query_gap_size >= ref_ins_end - MARGIN
				      and ref_curr_coord + gapless_block_size + query_gap_size <= ref_ins_end + MARGIN):
					
					results.add('\t'.join(site) + '\n')
					break

				# Case 3: move beyond insertion range
				elif ref_curr_coord > ref_ins_position + MARGIN:
					
					break

				# Account for gaps
				ref_curr_coord += query_gap_size

	for result in results:
		out.write(result)

	return


def get_sites_dict(sitesfile, start, end, mode):

	print "Loading sites."
	sites_dict = defaultdict(list)

	if mode == 'dd':
		col = 18
	elif mode == 'ii':
		col = 12


	with open(sitesfile, 'r') as f:
		
		for lineNum, line in enumerate(f.readlines(), 1):

			if lineNum < start: 
				continue

			if lineNum > end: 
				break

			line = line.split()

			if len(line) == 2: continue

			sites_dict[line[col]].append(line)

			print line[col]

	print "Sites loaded. \n"
	return sites_dict

def get_chain_dict(chainfile, sites_dict):

	print "Loading chains."
	begin_chainload = time.time()

	chainIDs = sites_dict.keys()

	# max/min ID to load
	maxID = str(max([int(chainID) for chainID in chainIDs]))
	minID = str(min([int(chainID) for chainID in chainIDs]))

	foundFirstChain = False
	loadedAllChains = False

	chain_dict = defaultdict(list)

	with open(chainfile, 'r') as f:

		# Initialize 'current chain'
		curr_chain_id = -1
		curr_chain = []
		num_chains = 0

		for line in f.readlines():

			line = [word.strip() for word in line.split()]

			# ignore comments and blank line at end of each chain
			if len(line) == 0 or line[0].startswith('#'): continue 

			################ Deal with line ################

			# Only start building dict once reached min chainID
			if foundFirstChain == False: 

				# Set to true once we've reached our first relevant chain
				if line[0] == 'chain' and line[12] == minID:
					foundFirstChain = True
					curr_chain_id = line[12]
				else:
					continue

			# In relevant section of chain files
			else:

				if line[0] == 'chain':

					# Add loaded chain to dictionary
					chain_dict[curr_chain_id] = curr_chain
					curr_chain = []
					# print "Adding chain to dictionary."

					# If we've reached designated limit, open new file
					if curr_chain_id == maxID:
						loadedAllChains = True
						break			

					curr_chain_id = line[12]

			curr_chain.append(line)

			################ Move to next line ################

		# Edge case: maxID is last ID in file
		if not loadedAllChains: chain_dict[curr_chain_id] = curr_chain

	print "Loaded chains in {} minutes.\n".format((time.time() - begin_chainload)/60)

	return chain_dict


def main():

	sitesfile = sys.argv[1]
	chainfile = sys.argv[2]
	outfile = sys.argv[3]
	mode = sys.argv[4]

	start = int(sys.argv[5])
	if len(sys.argv)==7:
		end = int(sys.argv[6])
	else:
		end = float('inf')

	sites_dict = get_sites_dict(sitesfile, start, end, mode)
	chain_dict = get_chain_dict(chainfile, sites_dict)

	if mode == 'di':
		ref_deletion_find_insertion(sites_dict, chain_dict, outfile)
	elif mode == 'dd':
		ref_deletion_find_deletion(sites_dict, chain_dict, outfile)
	elif mode == 'ii':
		ref_insertion_find_insertion(sites_dict, chain_dict, outfile)
	elif mode == 'id':
		ref_insertion_find_deletion(sites_dict, chain_dict, outfile)

if __name__ == '__main__':
	main()