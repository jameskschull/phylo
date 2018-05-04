# phylo

Approach: look at >50bp intronic insertions.

HG38 INTRONS: hg38, Genes, GENCODE v24, BED format, introns only.
--> introns.sorted.bed
951 423 lines

INTRONS LABELLED WITH MM10/CANFAM3 chainIDs:

Label introns with EnsemblIDs and sort by ID:

./label_and_sort_byEnsembl.sh ../mappings/ensembl_ids.bed ../sorted/introns.sorted.bed

N unique introns in original file: 
333364 ../sorted/introns.sorted.bed.unique.bed
N unique introns after labelling with overlapSelect: 
426213 temp_file.1.bed
After cleaning: 
426213 temp_file.2.bed
Final unique count after sorting: 
426213 ../sorted/introns.sorted.bed.byID

INTRONS LABELLED WITH CHAIN IDs:

For canFam3:
join -1 4 -2 1 introns.sorted.bed.byID canFam3.trimmed.chain_ids > introns.can.bed
315 509

ensemblID chr hgStart hgEnd chainID

For mm10:
join -1 4 -2 1 introns.sorted.bed.byID mm10.trimmed.chain_ids > introns.mm.bed
309 681

