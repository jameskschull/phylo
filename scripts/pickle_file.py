# Converts a chain/insertions file to a dictionary and serializes it using pickle

import pickle
from collections import defaultdict
import sys
import time
				
def main():

	file = sys.argv[1]
	filetype = sys.argv[2]

	if filetype == 'chain':
		num_chains_per_file = int(sys.argv[3]) # last chain to pickle (exclusive)

	startTime = time.time()
	with open(file, 'r') as f:

		if filetype == 'chain':
			print "Loaded chain file in {} seconds.".format((time.time() - startTime))
			print "Creating outfiles with {} chains per file.".format(num_chains_per_file)

		if filetype == 'chain':

			file_no = 1

			curr_chain_id = -1
			curr_chain = []
			num_chains = 0

			# initialize first file
			out = open("pickled/"+file+".{}".format(file_no), 'w')
			return_dict = defaultdict(list)

			for line in f.readlines():

				# Ignore empty lines
				if len(line.strip())==0: continue

				line = line.split()

				# Ignore commented lines
				if line[0][0] != '#':

					if line[0] == 'chain':

						# if not first chain
						if curr_chain_id > 0:
							# print "Adding chain."
							return_dict[curr_chain_id] = curr_chain
							curr_chain = []
							num_chains += 1

							# If we've reached designated limit, open new file
							if num_chains == num_chains_per_file:
								print "Chain limit reached, opening new file."
								pickle.dump(return_dict, out)
								out.close()

								num_chains = 0 # reset count
								file_no += 1
								return_dict = defaultdict(list)
								out = open(file+".pickled{}".format(file_no), 'w')
								print "File number {} opened".format(file_no)		

						curr_chain_id = line[12]
						curr_chain.append(line)
					else:
						curr_chain.append(line)

			return_dict[curr_chain_id] = curr_chain
			pickle.dump(return_dict, out)

		elif filetype == 'bed':

			out = open(file+".pickled", 'w')

			return_dict = defaultdict(list)

			for line in f.readlines():

				line = line.split()

				return_dict[line[4]].append(line)

			pickle.dump(return_dict, out)

		# q1 deletions labelled with q2 chain IDs
		elif filetype == 'q1dels':

			out = open(file+".pickled", 'w')

			return_dict = defaultdict(list)

			for line in f.readlines():

				line = line.split()

				return_dict[line[7]].append(line)

			pickle.dump(return_dict, out)


	return 

if __name__ == '__main__':
	main()




