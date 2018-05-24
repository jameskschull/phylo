# Takes a sitesfile and and number of lines per job, and creates a parasol
# joblist for get_indel_coords.py

# EXAMPLE USAGE: python create_job_list.py '../sorted/mm10/introns.mm.bed' '../chains/mm10/hg38.mm10.all.chain' '../jobs/indels/mm10dels_job/mm10dels' 'deletion' 15000

import sys

sitesfile_name = sys.argv[1]
chainfile_name = sys.argv[2]
jobfile_name = sys.argv[3]
type_of_search = sys.argv[4]
lines_per_job = int(sys.argv[5])

sitesfile = '/cluster/u/jschull/phylo/' + '/'.join(sitesfile_name.split('/')[1:])
chainfile = '/cluster/u/jschull/phylo/' + '/'.join(chainfile_name.split('/')[1:])

jobfile = open(jobfile_name, 'w')

num_lines = sum(1 for line in open(sitesfile_name, 'r'))

print "{} lines in sitesfile.".format(num_lines)

line_start = 1
line_end = line_start + lines_per_job - 1

while line_start < num_lines:

	jobfile.write('python /cluster/u/jschull/phylo/scripts/get_indel_coords.py {} {} \'{}\' {} {}'.format(sitesfile, chainfile, type_of_search, line_start, line_end) + '\n')
	line_start = line_end + 1
	line_end = line_start + lines_per_job - 1

print "Jobs written to {}.".format(jobfile_name)