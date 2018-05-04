# phylo

Approach: look at >50bp intronic insertions.

HG38 INTRONS: hg38, Genes, GENCODE v24, BED format, introns only.
--> introns.sorted.bed
951 423 lines

INTRONS LABELLED WITH MM10/CANFAM3 chainIDs:

Label introns with EnsemblIDs:

overlapSelect -mergeOutput ensembl_ids.bed introns.sorted.bed introns.withIDs.bed
canFam3.trimmed.chainIDs
mm10.trimmed.chainIDs

