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

## DEBUGGING ##

# python get_indel_coords.py ../sorted/mm10/introns.mm.bed.pickled ../chains/mm10/hg38.mm10.all.chain.pickled ../sorted/mm10/hg38plus_mm10minus.bed 'deletion'
# python get_indel_coords.py ../sorted/canFam3/introns.can.bed.pickled ../chains/canFam3/hg38.canFam3.all.chain.pickled ../sorted/canFam3/hg38plus_canFam3minus.bed 'deletion'

# OPTIMIZATIONS TO IMPLEMENT:
# - Only search chain once for each sites_dict[chainID]
# - Parallelize by splitting into multiple chains

# Call with pickled_sites, pickled_chain, outfile, target (either 
# 'deletion' or 'insertion').

def get_indel_coords(sites_dict, chain_dict, outfile, target):

	print "Getting query {}s.".format(target)
	
	print "Writing results to {}.".format(outfile)
	out = open(outfile, 'w')

	# contains column of gap we're looking for
	if target=="deletion":

		target_col = 1
		non_target_col = 2
	elif target=="insertion":
		target_col = 2
		non_target_col = 1

	chain_ids = sites_dict.keys()
	num_chains = len(chain_ids)
	num_sites_found = 0

	# Only keep unique sites
	found_sites = set([])

	# Search each chain for indels
	for c, chain_id in enumerate(chain_ids, 1):

		print "Searching chain {} out of {}.".format(c, num_chains)

		chain = chain_dict.get(chain_id, None)

		if chain == None: 
			print 'Chain not found!'
			continue

		ref_start = int(chain[0][5])
		query_start = int(chain[0][10])

		for site in sites_dict[chain_id]:

			print "Searching coordinates of new site: \n"
			print site

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
					if target=="deletion":
						insertion_end_ref = insertion_start_ref + insertion_size
						insertion_end_query = insertion_start_query
					elif target=="insertion":
						insertion_end_ref = insertion_start_ref
						insertion_end_query = insertion_start_query + insertion_size

					entry = ensemblID + '\t' + str(insertion_size) + '\t' + chain[0][2] + '\t' +  str(insertion_start_ref) + '\t' + str(insertion_end_ref)  + '\t' + chain[0][7] + '\t' + chain[0][9] + '\t' + str(insertion_start_query) + '\t' + str(insertion_end_query) + '\n'
					found_sites.add(entry)
					num_sites_found += 1
					print(entry)
					print line

				ref_offset += int(line[1]) if len(line) > 1 else 0
				query_offset += int(line[2]) if len(line) > 1 else 0

		print "Search completed. {} sites found.".format(num_sites_found)

	print "These are the unique sites found: '\n"
	for site in found_sites:
		print site
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

	chainIDs = sites_dict.keys()
	print "Num chain IDs: {}".format(len(chainIDs))
	maxID = max(sorted(chainIDs))
	loadedAllChains = False

	with open(chainfile, 'r') as f:

		chain_dict = defaultdict(list)

		curr_chain_id = -1
		curr_chain = []
		num_chains = 0

		for line in f.readlines():

			line = [word.strip() for word in line.split()]

			if len(line) == 0: continue # ignore blank line at end of each chain

			if line[0] == 'chain':

				# if not first chain
				if curr_chain_id > 0:
					chain_dict[curr_chain_id] = curr_chain
					curr_chain = []

					# If we've reached designated limit, open new file
					if curr_chain_id == maxID:
						loadedAllChains = True
						break			

				curr_chain_id = line[12]

			curr_chain.append(line)

		# Edge case: maxID is last ID in file
		if not loadedAllChains: chain_dict[curr_chain_id] = curr_chain

		print "Loaded chains in {} minutes.\n".format((time.time() - begin_chainload)/60)

		return chain_dict

def main():

	sitesfile = sys.argv[1]
	chainfile = sys.argv[2]
	outfile = sys.argv[3]
	target = sys.argv[4]
	start = int(sys.argv[5])

	# Load sites
	if len(sys.argv)==7:
		end = int(sys.argv[6])
	else:
		end = float('inf')

	name_parts = chainfile.replace('/', '.').split('.')

	# Build outfile name
	outfile = name_parts[5] + '.' + name_parts[6] + '.' + target + '.' + str(start) + '-' + str(end) + '.txt'
	
	sites_dict = get_sites_dict(sitesfile, start, end)

	# Load chains
	chain_dict = get_chain_dict(chainfile, sites_dict)

	
	get_indel_coords(sites_dict, chain_dict, outfile, target)

	return

if __name__ == '__main__':
	main()