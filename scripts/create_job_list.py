# Takes a sitesfile and and number of lines per job, and creates a parasol
# joblist for get_indel_coords.py

# EXAMPLE USAGE: python create_job_list.py 'indels' '../sorted/mm10/introns.mm10.bed' '../chains/mm10/hg38.mm10.all.chain' 'results/hg38.mm10.dels' '../jobs/indels/mm10dels_job/mm10dels' 0 500

# For compare_to_query:
# python create_job_list.py 'compare' '../sorted/mm10/hg38.mm10.dels.canFamID.txt' '../chains/canFam3/hg38.canFam3.all.chain' 'results/hg38.mm10del.canFam3ins' '../jobs/evidence/mm10del_canFam3ins_job/mm10dels_canFam3ins' 0 100
# python create_job_list.py 'compare' '../sorted/canFam3/hg38.canFam3.dels.mm10ID.txt' '../chains/mm10/hg38.mm10.all.chain' 'results/hg38.canFam3del.mm10ins' '../jobs/evidence/canFam3del_mm10ins_job/canFam3dels_mm10ins' 0 100
# python create_job_list.py 'compare' '../sorted/mm10/hg38.mm10.ins.canFamID.txt' '../chains/canFam3/hg38.canFam3.all.chain' 'results/hg38.mm10ins.canFam3ins' '/cluster/u/jschull/phylo/jobs/evidence/mm10ins_canFam3ins_job/mm10ins_canFam3ins' 1 100

#python create_job_list.py 'compare' '../sorted/loxAfr3/hg38.canFam3.dels.mm10.ins.txt.loxAfr3ID' '../chains/loxAfr3/hg38.loxAfr3.all.chain' 'results/hg38.canFam3del.mm10ins.loxAfr3del' '/cluster/u/jschull/phylo/jobs/evidence/canFam3del_mm10ins_loxAfr3del_job/canFam3del_mm10ins_loxAfr3del' 2 50
#python create_job_list.py 'compare' '../sorted/loxAfr3/hg38.mm10.dels.canFam3.ins.txt.loxAfr3ID' '../chains/loxAfr3/hg38.loxAfr3.all.chain' 'results/hg38.mm10del.canFam3ins.loxAfr3del' '/cluster/u/jschull/phylo/jobs/evidence/mm10del_canFam3ins_loxAfr3del_job/mm10del_canFam3ins_loxAfr3del' 2 50


# mm10del_canFam3ins_job'
import sys
from subprocess import call

script_type = sys.argv[1]
sitesfile_name = sys.argv[2]
chainfile_name = sys.argv[3]
outfile_name = sys.argv[4]
jobfile_name = sys.argv[5]
type_of_search = int(sys.argv[6])
lines_per_job = int(sys.argv[7])

if script_type == 'indels':
	script = 'get_indel_coords.py'
elif script_type == 'compare':
	script = 'compare_to_query.py'

sitesfile = '/cluster/u/jschull/phylo/' + '/'.join(sitesfile_name.split('/')[1:])
chainfile = '/cluster/u/jschull/phylo/' + '/'.join(chainfile_name.split('/')[1:])

jobfile = open(jobfile_name, 'w')

job_dir = jobfile_name.split('/')
job_dir = '/'.join(job_dir[:len(job_dir) - 1])

num_lines = sum(1 for line in open(sitesfile_name, 'r'))

print "{} lines in sitesfile.".format(num_lines)

line_start = 1
line_end = line_start + lines_per_job - 1

file_no = 1
while line_start < num_lines:

	if type_of_search in [0, 2]:

		jobfile.write('python /cluster/u/jschull/phylo/scripts/{} {} {} {} {} {} {}'.format(script, sitesfile, chainfile, outfile_name + str(line_start) + '-' + str(line_end) + '.txt', str(type_of_search), line_start, line_end) + '\n')
		line_start = line_end + 1
		line_end = line_start + lines_per_job - 1

	elif type_of_search == 1:

		script_name = job_dir + '/scripts/job'+str(file_no)+'.sh'
		print script_name
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


print "Jobs written to {}.".format(jobfile_name)

