# phylo

Approach: look at >50bp intronic insertions.

HG38 INTRONS: hg38, Genes, GENCODE v24, BED format, introns only.
--> introns.sorted.bed
951 423 lines

sed -i 's/chr23/chrX/' introns.sorted.bed
sed -i 's/chr24/chrY/' introns.unsort

INTRONS LABELLED WITH MM10/CANFAM3 chainIDs:

Label introns with EnsemblIDs and sort by ID:

./label_and_sort_byEnsembl.sh ../mappings/ensembl_ids.bed ../sorted/introns.sorted.bed

N unique introns in original file: 
333364 ../sorted/introns.sorted.bed.unique.bed
N unique introns after labelling with overlapSelect: 
439313 temp_file.1.bed
After cleaning: 
439313 temp_file.2.bed
Final unique count after sorting: 
439313 ../sorted/introns.sorted.bed.byID

INTRONS LABELLED WITH CHAIN IDs:

join <(sort file1.txt) <(sort file2.txt)

For canFam3:
join -1 4 -2 1 introns.sorted.bed.byID ../chains/canFam3/canFam3.trimmed.chain_ids > canFam3/introns.can.bed

324848 canFam3/introns.can.bed

ensemblID chr hgStart hgEnd chainID

For mm10:
join -1 4 -2 1 introns.sorted.bed.byID mm10.trimmed.chain_ids > introns.mm.bed
309 681

For monDom5:
join -1 4 -2 1 introns.sorted.bed.byID monDom5.trimmed.chain_ids > introns.monDom.bed

For loxAfr3:
join -1 4 -2 1 introns.sorted.bed.byID ../chains/loxAfr3/loxAfr3.trimmed.chain_ids > loxAfr3/introns.loxAfr3.bed


PICKLING:

python pickle_file.py ../chains/mm10/hg38.mm10.all.chain 'chain'
python pickle_file.py ../chains/canFam3/hg38.canFam3.all.chain 'chain'
python pickle_file.py ../chains/monDom5/hg38.monDom5.all.chain 'chain'

3 678 476 chains in mm10


DELS:

python get_indel_coords.py ../sorted/introns.can.bed.pickled ../chains/canFam3/hg38.canFam3.all.chain.pickled ../sorted/hg38insertions.vsCanFam3.bed

python get_indel_coords.py ../sorted/introns.mm.bed.pickled ../chains/mm10/hg38.mm10.all.chain.pickled ../sorted/hg38.mm10.dels.bed

python get_indel_coords.py ../sorted/introns.mm.bed.pickled ../chains/mm10/hg38.mm10.all.chain.pickled ../sorted/hg38insertions.vsmm10.bed


INDELS --> EVIDENCE

Label with q2 IDs
join -1 1 -2 1 <(sort -k1.5 -d ../../jobs/indels/mm10dels_job/hg38.mm10.dels.txt) ../../chains/canFam3/canFam3.trimmed.chain_ids > hg38.mm10.dels.canFamID.txt.unsort

join -1 1 -2 1 <(sort -k1.5 -d ../../jobs/indels/canFam3dels_job/hg38.canFam3.dels.txt) ../../chains/mm10/mm10.trimmed.chain_ids > hg38.canFam3.dels.mm10ID.txt

join -1 1 -2 1 <(sort -k1.5 -d ../results/hg38.mm10.dels.canFam3.ins.txt) ../chains/loxAfr3/loxAfr3.trimmed.chain_ids > ../sorted/loxAfr3/hg38.mm10.dels.canFam3.ins.txt.loxAfr3ID


-------------------------
 # Inspect chain to see if insertion is contained in gapless block
			# if ref_chain_start - MARGIN <= ref_insertion_start:

			# 	ref_curr_coord = ref_chain_start

			# 	line_n = 1

			# 	while ref_curr_coord - MARGIN <= ref_insertion_start and line_n < chain_len:

			# 		curr_line = chain[line_n]
			# 		gapless_block_size = int(curr_line[0])
			# 		ref_gap_size = int(curr_line[1])

			# 		#if the whole insertion is contained in a gapless block, keep it
			# 		if ref_curr_coord - MARGIN <= ref_insertion_start and curr_coord + gapless_block_size + MARGIN >= insertion_end_coord:
			# 			print "Insertion alignment found"
			# 			print '\t'.join(site) + '\n'
			# 			print 'Line in chainfile: {}'.format(curr_line)
			# 			print 'Curr coord: {}'.format(ref_curr_coord)
			# 			print 'End of gapless block: {}'.format(ref_curr_coord + gapless_block_size)
			# 			out.write('\t'.join(site) + '\n')

			# 		ref_curr_coord += gapless_block_size + ref_gap_size
			# 		line_n += 1

	# return 

-------------------------
Old compare (only looking at gapless blocks)

def querydel(sites_dict, chain_dict, outfile):

	print "Searching for evidence."
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

		ref_chain_start = int(chain[0][5])

		# For each querydel that maps to that chain
		for site in sites_dict[chain_id]:

			# reference (hg38) coordinates of the insertion
			ref_insertion_start = int(site[3]) 
			ref_insertion_end = int(site[4])

			ref_curr_coord = ref_chain_start

			# Search chain
			for line in chain[1:]:

				gapless_block_size = int(line[0])

				# print "Ref coord: {}, insertion start: {}, gapless_block_size: {}.".format(ref_curr_coord, ref_insertion_start, gapless_block_size)

				# If we have gone beyond the insertion, break
				if ref_curr_coord > ref_insertion_start: 
					# print "Gone past insertion start, breaking."
					break

				# If whole insertion is contained in a gapless block, keep it and end search
				if ref_curr_coord - MARGIN <= ref_insertion_start and ref_curr_coord + gapless_block_size + MARGIN >= ref_insertion_end:
					print "Evidence found."
					found_evidence.add('\t'.join(site) + '\n')

					break

				# If on last line of chain, end loop
				if len(line) == 1: 
					print "Last line of chain"
					continue

				query_gap_size = int(line[1])
				ref_curr_coord += gapless_block_size + query_gap_size

	for evidence in found_evidence:
		out.write(evidence)
		# print evidence

	return

	----------------
	if type_of_search in [0, 2]:

		jobfile.write('python /cluster/u/jschull/phylo/scripts/{} {} {} {} {} {} {}'.format(script, sitesfile, chainfile, outfile_name + str(line_start) + '-' + str(line_end) + '.txt', str(type_of_search), line_start, line_end) + '\n')
		line_start = line_end + 1
		line_end = line_start + lines_per_job - 1

	elif type_of_search == 1:

		script_name = job_dir + '/scripts/job'+str(file_no)+'.sh'
		with open(script_name, 'w') as f:
			f.write('#!/bin/sh' + '\n')
			f.write('export PYTHON_EGG_CACHE=results' + '\n')
			f.write('export PYTHONPATH=$PYTHONPATH:/cluster/u/yatisht/aligner/comparative/try/lib64/python2.7/site-packages/' + '\n')
			f.write('python /cluster/u/jschull/phylo/scripts/{} {} {} {} {} {} {}'.format(script, sitesfile, chainfile, outfile_name + str(line_start) + '-' + str(line_end) + '.txt', str(type_of_search), line_start, line_end) + '\n')

		jobfile.write(script_name + '\n')
		# call('chmod +x ' + script_name)

		line_start = line_end + 1
		line_end = line_start + lines_per_job - 1

		file_no += 1
