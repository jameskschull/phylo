#!/bin/sh

## 

origFile=$1
otherSpecies=$2

sort -k1.4,1n -k2,2n $origFile > $origFile.sorted
echo "Sorted"
overlapSelect -mergeOutput -inFmt=bed ../sorted/ensembl_ids.bed $origFile.sorted $origFile.labelled
echo "Labelled with EnsemblIDs"
python file_manip.py
echo "Trimmed"
sort -k5.5 -d $origFile.trimmed > $origFile.byID
echo "Sorted by ID"
join -1 5 -2 1 $origFile.byID ../sorted/$otherSpecies.trimmed.chain_ids > $origFile.$otherSpecies.chains
echo "Labelled with chain IDs"

rm $origFile.sorted
rm $origFile.labelled
rm $origFile.trimmed
rm $origFile.byID
