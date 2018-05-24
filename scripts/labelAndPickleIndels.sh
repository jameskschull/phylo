#!/bin/sh

# Takes a file of query (q1) insertions/deletions, as well as the name of another query (q2)
# species, and labels the file with the appropriate chain IDs. Then pickles file. 

# e.g for mm10 dels: ./labelAndPickleIndels.sh ../sorted/mm10/hg38plus_mm10minus.bed canFam3

infile=$1
q2=$2

q2_chain_ids='../chains/'$q2'/'$q2'.trimmed.chain_ids'
printf 'Infile: '$infile'\n'
printf 'Q2 chain ID file: '$q2_chain_ids'\n \n'

# Sort infile by EnsemblID
sort -k1.5,1 -d $infile > $infile.byID 
printf '1. Sorted infile by EnsemblID. \n \n'

# Map ensemblID to q2 chainID
join -1 1 -2 1 $infile.byID $q2_chain_ids > $infile.$q2
printf '2. Joined infile with chain IDs, new file named: '$infile.$q2'.\n \n'

rm $infile.byID

printf '3. Pickling file. \n \n'
# Pickle labelled file
python pickle_file.py $infile.$q2 'q1dels'

