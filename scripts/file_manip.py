# Removes from BED file all HAP/Un/etc chroms
import sys

def rename_chroms(index, infile, outfile):

	with open(infile, 'r') as f:

		out = open(outfile, 'w')

		lines = [line.split(",") for line in f.readlines()]

		# X chr = 23
		# Y chr = 24
		for line in lines:
			chrn = line[0][index:] 

			if len(line[0]) <= 5:

				if chrn == 'X':
					out.write('chr23' + '\t' + '\t'.join(line[1:4]))
				elif chrn == 'Y':
					out.write('chr24' + '\t' + '\t'.join(line[1:4]))
				else:
					out.write('chr' + chrn + '\t' + '\t'.join(line[1:4]))

def clean_target_file(infile, outfile):

	with open(infile, 'r') as f:
		out = open(outfile, 'w')

		# out.write('intronChr' + '\t' + 'TEstart' + '\t' + 'TEend' + '\t' + 'INTstart' + '\t' + 'INTend' + '\t' + 'EnsemblID' + '\n')

		lines = [line.split() for line in f.readlines()]
		# print lines[0]

		for i, line in enumerate(lines):
			new_line = line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[6] + '\n'
			out.write(new_line)
			# if i==0: print new_line

def main():

	argv = sys.argv[1:]
	filename = argv[1]
	outfilename = argv[2]

	if argv[0] == 'clean':
		clean_target_file(filename, outfilename)

	# for TES
	#rename_chroms(3, '../unprocessed/hg38_TEs.bed', '../processed/hg38_TEs_proc.bed')
	
	# for ensembl, introns, TEs
	# rename_chroms(0, '../unproc/ensembl_ids.preproc.bed', '../unproc/ensembl_ids.proc.bed')
	# rename_chroms(3, '../unprocessed/introns_unproc.bed', '../sorted/introns.bed')
	# rename_chroms(3, '../sorted/hg38_TEs_uncut.bed', '../sorted/hg38s_TEs.bed')

	# Clean target loci file

	# clean_target_file("../sorted/target_sites.uncut.bed", "../sorted/target_sites.bed")

	# Clean insertions files
	#clean_target_file("../unproc/hg38.canFam3.insertions.bed.labelled", '../unproc/hg38.canFam3.insertions.bed.trimmed')

if __name__ == '__main__':
	main()