# Given a file of select reference sites to search, a query chain, an outfile,
# and an argument indicating which sequence to search for insertions, this
# writes to an outfile with the coordinates of all insertions or deletions (with reference to query),
# depending on the argument provided.

# EXAMPLE #
# If reference is hg38, and query is mm10, pass as args: target_sites.pickled, hg38.mm10.all.chain.pickled, 
# hg38.mm10.deletions.bed, "deletion"

# The outfile will contain all loci where there has been an apparent SINE/LINE insertion in human (i.e. a deletion
# in the query).

import sys
import subprocess 
import time
import pickle

# Call with pickled_insertions, pickled_chain, outfile, target (either 
# 'deletion' or 'insertion').
def get_indel_coords(pickled_sites, pickled_chain, outfile, target):
	
	out = open(outfile, 'w')

	t1 = time.time()
	sites_dict = pickle.load(open(pickled_sites, 'r'))
	t2 = time.time()
	print ("Loaded introns in {} minutes".format((t2-t1)/60))
	print sites_dict.keys()

	chain_dict = pickle.load(open(pickled_chain, 'r'))
	t3 = time.time()
	print ("Loaded chains in {} minutes".format((t3-t2)/60))
	print chain_dict.keys()

	# print ("{} chainIDs to search".format(num_ids))

	# contains column of gap we're looking for
	if target=="deletion":
		target_col = 1
		non_target_col = 2
	elif target=="insertion":
		target_col = 2
		non_target_col = 1

	# Search each chain for indels
	for c, chain_id in enumerate(sites_dict.keys()):

		print ("Searching for chainID {}".format(chain_id))

		start_time = time.time()
		chain = chain_dict.get(chain_id, None)

		if chain == None: 
			print 'Chain not found!'
			continue

		ref_start = int(chain[0][5])
		query_start = int(chain[0][10])

		for site in sites_dict[chain_id]:

			site_start = int(site[2])
			site_end = int(site[3])

			ref_offset = 0 # offset from start of chain
			query_offset = 0

			for i, line in enumerate(chain[1:], 1):

				# if we have stepped beyond the end of the site, break
				if ref_start + ref_offset > site_end:
					break

				# if we haven't reached the beginning of the site, skip line
				if ref_start + ref_offset < site_start:
					continue

				gapless_block_size = int(line[0])
				ref_offset += gapless_block_size
				query_offset += gapless_block_size

				# either insertion in reference or query
				insertion_size = int(line[target_col]) 

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

					print ("Writing")
					out.write(str(insertion_size) + '\t' + chain[0][2] + '\t' +  str(insertion_start_ref) + '\t' + str(insertion_end_ref)  + '\t' + chain[0][7] + '\t' + str(insertion_start_query) + '\t' + str(insertion_end_query) + '\n')

				ref_offset += int(line[1])
				query_offset += int(line[2])
			
		end_time = time.time()

		# if c == 0: 
			# print "Time for first chainID: {} minutes".format((end_time - start_time)/60)
			# print "Rough projected time: {} hours".format((end_time - start_time)*num_ids/60/60)
							
	return


def main():

	pickled_sites = sys.argv[1]
	pickled_chain = sys.argv[2]
	outfile = sys.argv[3]
	target = sys.argv[4]

	get_indel_coords(pickled_sites, pickled_chain, outfile, target)

	return

if __name__ == '__main__':
	main()