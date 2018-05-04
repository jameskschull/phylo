# Given chain IDs of hg-query chains with overlapping hg TEs, 
# finds coordinates of all gaps in the human sequence, implying insertion
# in hg, and no insertion in query

# This assumes that outgroup has no insertion - if they do, the results
# may indicate a deletion in the query 

# 

import sys
import subprocess 
import time

def get_target_chain_file(target_sites, chainfile, outfile):

	with open(target_sites, 'r') as f:
		with open(chainfile, 'r') as c:
			out = open(outfile, 'w')

			target_ids = set([line.split()[6] for line in f.readlines()])

			# Search each chain for indels
			for c, target_id in enumerate(target_ids):

				print "Finding chain {}".format(c)
				chain = subprocess.check_output(['chainFilter', '-id={}'.format(target_id), '{}'.format(chainfile)])
				out.write(chain)

	return

# Call with target_sites, chainfile, outfile
def get_indel_coords(target_sites, chainfile, outfile):

	with open(target_sites, 'r') as f:
		with open(chainfile, 'r') as c:
			out = open(outfile, 'w')

			target_ids = set([line.split()[6] for line in f.readlines()])
			num_ids = len(target_ids)
			print "{} IDs to search".format(num_ids)

			# Search each chain for indels
			for c, target_id in enumerate(target_ids):

				start_time = time.time()

				print "Searching chain {}".format(c)

				chain = subprocess.check_output(['chainFilter', '-id={}'.format(target_id), '{}'.format(chainfile)])
				lines = [line.split() for line in chain.splitlines()]

				offset = 0 # offset from start of chain
				for line in lines[1:len(lines)-2]:
					
					offset += int(line[0])
					insertion_size = int(line[1])

					# if there is an hg38 insertion of >= 50bp, record it
					if insertion_size >= 50:

						insertion_start = int(lines[0][5]) + offset
						insertion_end = int(lines[0][5]) + offset + insertion_size
						out.write(lines[0][2] + '\t' + str(insertion_start) + '\t' + str(insertion_end) + '\t' + str(insertion_size) + '\n')

					offset += int(line[1])
				
				end_time = time.time()

				if c == 0: print "Projected time: {} hours".format((end_time - start_time)*num_ids/60/60)

				
	return


def main():

	target_sites = sys.argv[1]
	chainfile = sys.argv[2]
	outfile = sys.argv[3]

	# get_target_chain_file(target_sites, chainfile, outfile)
	get_indel_coords(target_sites, chainfile, outfile)

	return

if __name__ == '__main__':
	main()