# Converts a chain/insertions file to a dictionary and serializes it using pickle

import pickle
from collections import defaultdict
import sys

file = sys.argv[1]
filetype = sys.argv[2]

if filetype == 'chain':
	num_chains_per_file = sys.argv[3] # last chain to pickle (exclusive)

with open(file, 'r') as f:

	if filetype == 'chain':

		file_no = 1

		curr_chain_id = -1
		curr_chain = []
		num_chains = 0

		# initialize first file
		out = open(file+".pickled{}".format(file_no), 'w')
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
						# print "Adding chain"
						return_dict[curr_chain_id] = curr_chain
						curr_chain = []
						num_chains += 1
						
						# If we've reached designated limit, open new file
						if num_chains == num_chains_per_file:
							pickle.dump(return_dict, out)
							file_no += 1
							out.close()
							out = open(file+".pickled{}".format(file_no), 'w')
							return_dict = defaultdict(list)

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
				




