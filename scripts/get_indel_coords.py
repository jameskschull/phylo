# Given a file of select reference sites to search, a query chainfile,
# an argument indicating which sequence to search for insertions, and a start and end line number
# (read only these lines) this writes to an outfile with the coordinates of all insertions or 
# deletions (with reference to query), depending on the argument provided.

# EXAMPLE #
# If reference is hg38, and query is mm10, pass as args: target_sites.pickled, hg38.mm10.all.chain.pickled, 
# hg38.mm10.deletions.bed, "deletion"

# The outfile will contain all loci where there has been an apparent SINE/LINE insertion in human (i.e. a deletion
# in the query).

# Args format:
# python get_indel_coords.py SELECTFILE CHAINFILE IN/DEL START (END)

# Outfile format:
# INSERTION SIZE, REF CHROM, REF START, REF END, QUERY CHROM, QUERY START, QUERY END

import sys
import subprocess 
import time
import pickle
from collections import defaultdict

def get_indels(sites_dict, chain_dict, outfile, target, target_col, non_target_col):

	# initialize variables for para
	insertion_start_ref = -1
	insertion_start_query = -1
	insertion_end_ref = -1
	insertion_end_query = -1

	print "Getting query {}s.".format(target)
	
	print "Writing results to {}.".format(outfile)
	out = open(outfile, 'w')

	chain_ids = sites_dict.keys()
	num_chains = len(chain_ids)
	num_sites_found = 0

	# Only keep unique sites
	found_sites = set([])

	# Search each chain for indels
	for c, chain_id in enumerate(chain_ids, 1):

		print "Searching chain {} out of {}.".format(c, num_chains)

		chain = chain_dict.get(chain_id, None)

		# print "CHAIN ID: {}".format(chain_id)

		if chain == None: 
			print 'Chain not found!'
			continue

		ref_start = int(chain[0][5])
		query_start = int(chain[0][10])

		for site in sites_dict[chain_id]:

			ensemblID = site[0]

			site_start = int(site[2])
			site_end = int(site[3])

			# offset from start of chain
			ref_offset = 0 
			query_offset = 0

			for i, line in enumerate(chain[1:], 1):

				if line == '': continue

				gapless_block_size = int(line[0])
				ref_offset += gapless_block_size
				query_offset += gapless_block_size

				# if we have stepped beyond the end of the site, break
				if ref_start + ref_offset > site_end:
					break

				# if we haven't reached the beginning of the site, skip line and update coords
				if ref_start + ref_offset < site_start:
					ref_offset += int(line[1]) if len(line) > 1 else 0
					query_offset += int(line[2]) if len(line) > 1 else 0
					continue

				# either insertion in reference or query
				insertion_size = int(line[target_col]) if len(line) > 1 else 0

				# if there is an insertion of >= 50bp, record it
				if insertion_size >= 50 and line[non_target_col] == '0':

					insertion_start_ref = ref_start + ref_offset
					insertion_start_query = query_start + query_offset

					# insertion start/end is same coord for 'deletion' genome
					if target==0: #deletion
						insertion_end_ref = insertion_start_ref + insertion_size
						insertion_end_query = insertion_start_query
					elif target==1: #insertion
						insertion_end_ref = insertion_start_ref
						insertion_end_query = insertion_start_query + insertion_size

					entry = ensemblID + '\t' + str(insertion_size) + '\t' + chain[0][2] + '\t' +  str(insertion_start_ref) + '\t' + str(insertion_end_ref)  + '\t' + chain[0][7] + '\t' + chain[0][9] + '\t' + str(insertion_start_query) + '\t' + str(insertion_end_query) + '\n'
					found_sites.add(entry)
					num_sites_found += 1
					# print(entry)
					# print line

				ref_offset += int(line[1]) if len(line) > 1 else 0
				query_offset += int(line[2]) if len(line) > 1 else 0

		print "Search completed. {} sites found.".format(num_sites_found)

	# print "These are the unique sites found: '\n"
	for site in found_sites:
		# print site
		out.write(site)
							
	return

# returns dictionary with (chainID : list of sites) items
def get_sites_dict(sitesfile, start, end):

	print "Loading sites."
	sites_dict = defaultdict(list)

	with open(sitesfile, 'r') as f:
		
		for lineNum, line in enumerate(f.readlines(), 1):

			if lineNum < start: continue

			if lineNum > end: break

			line = line.split()

			sites_dict[line[4]].append(line)

	print "Sites loaded. \n"
	return sites_dict


# returns dictionary of (chainID : chain) items 
def get_chain_dict(chainfile, sites_dict):

	print "Loading chains."
	begin_chainload = time.time()

	# strings
	chainIDs = sites_dict.keys()
	# print "Num chain IDs: {}".format(len(chainIDs))
	# print "ChainIDs to load: {}.".format(chainIDs)

	maxID = str(max([int(chainID) for chainID in chainIDs]))
	# print "Max ID to load: {}.".format(maxID)

	minID = str(min([int(chainID) for chainID in chainIDs]))
	# print "Min ID to load: {}.".format(minID)

	foundFirstChain = False
	loadedAllChains = False

	with open(chainfile, 'r') as f:

		chain_dict = defaultdict(list)

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

	sitesfile = sys.argv[1]
	chainfile = sys.argv[2]
	outfile = sys.argv[3]
	target = int(sys.argv[4])
	start = int(sys.argv[5])

	# TARGET:
	# 1 is insertion
	# 0 is deletion

	# Load sites
	if len(sys.argv)==7:
		end = int(sys.argv[6])
	else:
		end = float('inf')

	target_col = 1000
	non_target_col = 2000

	# contains column of gap we're looking for
	if target==0:
		target_col = 1
		non_target_col = 2
	elif target==1:
		target_col = 2
		non_target_col = 1

	sites_dict = get_sites_dict(sitesfile, start, end)

	# Load chains
	chain_dict = get_chain_dict(chainfile, sites_dict)

	get_indels(sites_dict, chain_dict, outfile, target, target_col, non_target_col)

	return

if __name__ == '__main__':
	main()