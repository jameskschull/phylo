#!/bin/sh

# Takes as args: [ensemblIDs] [infile]

ensemblIDs=$1
orig=$2

uniq $orig > $orig.unique.bed
echo "N unique introns in original file: "
wc -l $orig.unique.bed

overlapSelect -mergeOutput $ensemblIDs $orig.unique.bed temp_file_with_dups.1.bed
echo "N unique introns after labelling with overlapSelect: "
uniq temp_file_with_dups.1.bed > temp_file.1.bed
wc -l temp_file.1.bed
#head temp_file.1.bed

python file_manip.py 'clean' temp_file.1.bed temp_file.2.bed
echo "After cleaning: "
wc -l temp_file.2.bed
#head temp_file.2.bed

sort -k4.5 -d temp_file.2.bed > temp_file.3.bed
echo "Final unique count after sorting: "
uniq temp_file.3.bed > $orig.byID
wc -l $orig.byID

rm $orig.unique.bed
rm temp_file.1.bed
rm temp_file.2.bed
rm temp_file.3.bed
