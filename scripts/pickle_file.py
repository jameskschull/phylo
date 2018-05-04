# Converts a chain/insertions file to a dictionary and serializes it using pickle

import pickle
from collections import defaultdict
import sys

file = sys.argv[1]
filetype = sys.argv[2]

with open(file, 'r') as f:

	out = open(file+".pickled", 'w')

	return_dict = defaultdict(list)

	if filetype == 'chain':

		curr_chain_id = -1
		curr_chain = []

		for line in f.readlines():

			line = line.split()

			if len(line) > 1:
				if line[0] == 'chain':
					# if not first chain
					if curr_chain_id > 0:
						# print "Adding chain"
						return_dict[curr_chain_id] = curr_chain
						curr_chain = []
					# not first chain
					curr_chain_id = line[12]
					curr_chain.append(line)

		return_dict[curr_chain_id] = curr_chain

	elif filetype == 'insertion':

		for line in f.readlines():

			line = line.split()

			return_dict[line[4]] = line

	pickle.dump(return_dict, out)
				




