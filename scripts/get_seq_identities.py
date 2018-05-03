# Given bed file with coordinates of hg38 insertions, compare with other query 
# sequence to see if the insertion is present/absent there

# USAGE: python get_seq_identities.py HG_INSERTIONS_FILE QUERY_CHAIN_FILE OUTFILE

# Step 1: overlapSelect to get ensemblID 
# Step 2: Join with chainID for opposite species
# Step 3: Find matching coordinates in hg38 chains
# Step 4: See if there's a gap; if it's present, write to file
# Step 5: From other file, compare regions

import sys
import subprocess
import random
import time

# BUGGY - finding insertions in human only #

def main():

	hg_insertions = sys.argv[1]
	query_chainfile = sys.argv[2]
	outfile = sys.argv[3]

	with open(hg_insertions, 'r') as f:

		out = open(outfile, 'w') 
		num_insertions = int(subprocess.check_output(['wc', '-l', '{}'.format(hg_insertions)]).split()[0])

		# For each insertion
		for c, insertion in enumerate([insertion.split() for insertion in random.sample(f.readlines(), 100)]):
			
			start_time = time.time()

			print "Checking insertion {}".format(c)
			chain_id = insertion[5]
			insertion_start_coord = int(insertion[2])
			insertion_end_coord = int(insertion[2])

			# Compare to chain
			chain = subprocess.check_output(['chainFilter', '-id={}'.format(chain_id), '{}'.format(query_chainfile)])
			
			found_chain = time.time()
			print "Time to find chain {}".format(found_chain-start_time)

			lines = [line.split() for line in chain.splitlines()]

			chain_start_coord = int(lines[0][5])
			# print lines[0]

			# Inspect if chain includes whole insertion area 
			if chain_start_coord <= insertion_start_coord:

				curr_coord = chain_start_coord # start from beginning of chain

				line_n = 1

				while curr_coord -5 <= insertion_start_coord:

					curr_line = lines[line_n]
					# print lines[line_n]
					# print "Current coord: {}".format(curr_coord)
					# print "Block end coord: {}".format(curr_coord + int(curr_line[0]))

					# if whole insertion is contained in gapless block, keep
					if curr_coord - 5 <= insertion_start_coord and curr_coord + int(curr_line[0]) + 5 >= insertion_end_coord:
						print "Insertion alignment found"
						print '\t'.join(insertion) + '\n'
						print 'Line in chainfile: {}'.format(curr_line)
						print 'Curr coord: {}'.format(curr_coord)
						print 'End of gapless block: {}'.format(curr_coord + int(curr_line[0]))
						out.write('\t'.join(insertion) + '\n')

					curr_coord = curr_coord + int(curr_line[0]) + int(curr_line[1])

					line_n += 1


			end_time = time.time()
			if c == 0: print "Projected time: {} hours".format((end_time-start_time)*num_insertions/60/60)
			if c == 20: break

if __name__ == '__main__':
	main()
			





