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

# python compare_to_query.py ../sorted/mm10/hg38plus_mm10minus.bed.canFam3.pickled ../chains/canFam3/hg38.canFam3.all.chain.pickled ../sorted/evidence/hg38plus_mm10minus_canFam3plus.bed 'deletion'

# python compare_to_query.py ../sorted/loxAfr3/hg38.canFam3.dels.mm10.ins.txt.loxAfr3ID ../chains/loxAfr3/hg38.loxAfr3.all.chain 

import sys
import time
import pickle
import edlib
from Bio import SeqIO
from collections import defaultdict

MARGIN = 5
UPPER_THRESHOLD = 0.8 
LOWER_THRESHOLD = 0.7

# for filtering outgroup
def search_for_q2_deletions(sites_dict, chain_dict, outfile):

	print "Searching for insertions in query2."
	print "Writing results to {}.".format(outfile)

	out = open(outfile, 'w')

	chain_ids = sites_dict.keys()
	found_evidence = set([])

	# For each chain mapped to in the querydels file
	for chain_id in chain_ids:

		chain = chain_dict.get(chain_id, None)

		if chain == None: 
			print 'Chain not found!'
			continue
		else:
			print "Chain {} found.".format(chain_id)

		chain_len = len(chain)
		ref_chain_start = int(chain[0][5])

		# For each querydel that maps to that chain
		for site in sites_dict[chain_id]:

			ref_insertion_start = int(site[3])
			ref_insertion_end = int(site[4])
			ref_curr_coord = ref_chain_start

			site_missing = True

			# Search chain
			for line in chain[1:]:

				print "Current reference coordinate: {}".format(ref_curr_coord)

				gapless_block_size = int(line[0])

				print "Gapless block size: {}".format(gapless_block_size)
				print "Current allowable range: {}-{} + \n".format(ref_curr_coord - MARGIN, ref_curr_coord + gapless_block_size + MARGIN)

				# print "Ref coord: {}, insertion start: {}, gapless_block_size: {}.".format(ref_curr_coord, ref_insertion_start, gapless_block_size)

				# If we have gone beyond the insertion, break
				if ref_curr_coord - MARGIN > ref_insertion_start: 
					print "Gone past insertion start, breaking."
					break

				# If whole insertion is contained in a gapless block, site is not evidence
				if ref_curr_coord - MARGIN <= ref_insertion_start and ref_curr_coord + gapless_block_size + MARGIN >= ref_insertion_end:
					print "Insertion found in query 2."
					site_missing = False
					break

				# If on last line of chain, end loop
				if len(line) == 1: 
					print "Last line of chain"
					continue

				query_gap_size = int(line[1])
				ref_curr_coord += gapless_block_size + query_gap_size

			if site_missing == True:
				print "Insertion missing in query2: evidence found."
				found_evidence.add('\t'.join(site) + '\n')

	for evidence in found_evidence:
		out.write(evidence)

	return

# For hg + q1 - q2 + case and filtering by outgroup
def compare_to_q2_chain(sites_dict, chain_dict, outfile, mode):

	print "Searching for evidence, looking for {}s in second query.".format(mode)
	print "Writing results to {}.".format(outfile)

	out = open(outfile, 'w')

	chain_ids = sites_dict.keys()
	found_evidence = set([])

	# For each chain mapped to in the querydels file
	for chain_id in chain_ids:

		chain = chain_dict.get(chain_id, None)

		if chain == None: 
			print 'Chain not found!'
			continue
		else:
			print "Chain {} found.".format(chain_id)

		# Reference coordinate of chain start
		ref_chain_start = int(chain[0][5])

		# For each site that maps to that chain
		for site in sites_dict[chain_id]:

			rangeOver = False

			ref_curr_coord = ref_chain_start
			ref_ins_start, ref_ins_end = int(site[3]), int(site[4]) # reference coords of insertion
			gapless_bp, gap_bp = 0, 0 # number of bp within insertion range that are gapless/gap

			################ SEARCH CHAIN ################

			for line in chain[1:]:

				# Edge case: last line of chain has only gapless block
				if len(line) < 3:
					break

				gapless_block_size = int(line[0])
				query_gap_size = int(line[1])

				# If moved past insertion, stop search
				if rangeOver == True:
					break

				########## PROCESS GAPLESS BLOCK AND QUERY GAP ##########

				for i, block_size in enumerate([gapless_block_size, query_gap_size]):

					# Case 1: Insertion range starts within the block
					if ref_curr_coord + block_size > ref_ins_start and ref_curr_coord <= ref_ins_start:

						# i): only part of insertion range is contained within block
						if ref_curr_coord + block_size < ref_ins_end:
							if i == 0: 
								gapless_bp += ref_curr_coord + block_size - ref_ins_start
							elif i == 1:
								gap_bp += ref_curr_coord + block_size - ref_ins_start

						# ii): whole insertion range is contained within block
						elif ref_curr_coord + block_size >= ref_ins_end:
							rangeOver = True
							if i == 0:
								gapless_bp += ref_ins_end - ref_ins_start
							elif i == 1:
								gap_bp += ref_ins_end - ref_ins_start
							

					# Case 2: We are already in the range of the insertion
					elif ref_curr_coord + block_size > ref_ins_start and ref_curr_coord > ref_ins_start:

						# i) only part of insertion range is contained within block
						if ref_curr_coord + block_size < ref_ins_end:
							if i == 0:
								gapless_bp += block_size
							elif i == 1:
								gap_bp += block_size

						# ii) rest of insertion range is contained within block:
						elif ref_curr_coord + block_size >= ref_ins_end:
							rangeOver = True
							if i == 0:
								gapless_bp += ref_ins_end - ref_curr_coord
							elif i == 1:
								gap_bp += ref_ins_end - ref_curr_coord

					if rangeOver == False:
						ref_curr_coord += block_size
					else: 
						break

				#########################################################

			if gapless_bp == 0 and gap_bp == 0: continue

			print "Gapless bp: {}, gap bp: {}. Total insertion size: {}.".format(gapless_bp, gap_bp, ref_ins_end-ref_ins_start)
			gapless_percentage = gapless_bp/float(gapless_bp + gap_bp)

			if gapless_percentage > UPPER_THRESHOLD and mode=='insertion':
				print "Insertion found! Gapless percentage: {}".format(gapless_percentage)
				found_evidence.add('\t'.join(site) + '\n')

			elif gapless_percentage < LOWER_THRESHOLD and mode=='deletion':
				print "Deletion found! Gapless percentage: {}".format(gapless_percentage)
				found_evidence.add('\t'.join(site) + '\n')
			
	for evidence in found_evidence:
		out.write(evidence)

	return

# 1. Use the query (mouse) coordinates of the insertion to get the mouse sequence
# 2. Use the reference (human) coordinates (where start and end are actually the same) to get the dog chain and find the corresponding dog coordinates
# 3. Use those coordinates to get the dog sequence
# 4. Find edit distance between dog/mouse sequence

def double_insertion(sites_dict, chain_dict, outfile):

	out = open(outfile, 'w')
	print "Writing evidence to {}/.".format(outfile)

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

		q2_chr = chain[0][7]
		q2_strand = chain[0][9]
		q2_chain_start = int(chain[0][10])
		ref_chain_start = int(chain[0][5])
		
		# For each insertion that maps to that chain
		for site in sites_dict[chain_id]:

			print '\t'.join(site)

			insertion_size = int(site[1])

			# GET Q1 SEQUENCE
			q1_chr = site[5]
			q1_strand = site[6]
			q1_start = int(site[7])
			q1_end = int(site[8])
			q1_seq = query1_whole_genome[q1_chr][q1_start:q1_end]

			print "Mouse coordinates (from site): {}: {}-{} (strand: {})".format(q1_chr, q1_start, q1_end, q1_strand)

			# GET Q2 SEQUENCE
			ref_ins_position = int(site[3]) # since this is an insertion, ref start and end are the same
			print "Human insertion start position: {}".format(ref_position)

			ref_left = ref_ins_position - ref_chain_start # bp to ref start point
			q2_curr_coord = q2_chain_start

			for line in chain[1:]:

				#if on last chain line, break
				if len(line) < 3: break

				gapless_block_size = int(line[0])
				query_gap_size = int(line[1])
				ref_gap_size = int(line[2])

				if ref_left - gapless_block_size <= 0:
					q2_start = q2_curr_coord + ref_left
					break

				ref_left -= gapless_block_size
				q2_curr_coord += gapless_block_size

				if ref_left - query_gap_size <= 0:
					q2_start = q2_curr_coord + ref_left
					break

				ref_left -= query_gap_size
				q2_curr_coord += ref_gap_size

			q2_seq = query2_whole_genome[q2_chr][q2_start:q2_start + insertion_size]

			# Account for strand 
			if q2_strand == '-':
				q2_seq = q2_seq.reverse_complement()

			if q1_strand == '-':
				q1_seq = q1_seq.reverse_complement()

			q1_seq, q2_seq = str(q1_seq), str(q2_seq)

			print "Mouse sequence: {}".format(q1_seq)
			print "Dog sequence: {}".format(q2_seq)

			# Calculate similarity
			len_longer_seq = max(len(q1_seq), len(q2_seq))

			similarity = (len_longer_seq - int(edlib.align(q1_seq, q2_seq)["editDistance"]))/float(len_longer_seq)

			# compare to threshold, write if they're similar enough
			if similarity > UPPER_THRESHOLD:
				print "Sites have similarity of {}: evidence found! \n".format(similarity)
				out.write('\t'.join(site) + '\t' + 'similarity: {}'.format(round(similarity, 2)) + '\n')
			else:
				print "Sites have similarity of {}: insufficient. \n".format(similarity)

	return

# returns dictionary with (chainID : list of sites) items
def get_sites_dict(sitesfile, start, end):

	print "Loading sites."
	sites_dict = defaultdict(list)

	with open(sitesfile, 'r') as f:
		
		for lineNum, line in enumerate(f.readlines(), 1):

			if lineNum < start: 
				continue

			if lineNum > end: 
				break

			line = line.split()

			sites_dict[line[len(line) - 1]].append(line)

	print "Sites loaded. \n"
	return sites_dict


# returns dictionary of (chainID : chain) items 
def get_chain_dict(chainfile, sites_dict):

	print "Loading chains."
	begin_chainload = time.time()

	# strings
	chainIDs = sites_dict.keys()
	# print "Num chain IDs: {}".format(len(chainIDs))

	maxID = str(max([int(chainID) for chainID in chainIDs]))
	# print "Max ID to load: {}.".format(maxID)

	minID = str(min([int(chainID) for chainID in chainIDs]))
	# print "Min ID to load: {}.".format(minID)

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

		# print "Loaded these keys: {}.".format(sorted(chain_dict.keys()))
	return chain_dict

def main():

	sitesfile = sys.argv[1] # file containing indels
	chainfile = sys.argv[2] # name of query2 (to compare to)
	outfile = sys.argv[3] # file to write evidence to

	# 0 = deletion in query 1 (look for insertion in q2), 
	# 1 = insertion in query 1 (look for insertion in q2)
	# 2 = __ in query 1 (look for deletion in q2)
	mode = int(sys.argv[4]) 

	# first line to look at in sitesfile
	start = int(sys.argv[5]) 

	# last line to look at in sitesfile
	if len(sys.argv)==7:
		end = int(sys.argv[6])
	else:
		end = float('inf')

	sites_dict = get_sites_dict(sitesfile, start, end)
	chain_dict = get_chain_dict(chainfile, sites_dict)

	# hg + q1 - q2 +
	if mode == 0:
		compare_to_q2_chain(sites_dict, chain_dict, outfile, 'insertion')

	# hg - q1 + q2 +
	elif mode == 1:
		double_insertion(sites_dict, chain_dict, outfile)

	# hg + q1 - q2 + outgroup -
	elif mode == 2:
		compare_to_q2_chain(sites_dict, chain_dict, outfile, 'deletion')

	return


if __name__ == '__main__':
	main()