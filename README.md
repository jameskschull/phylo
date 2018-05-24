# phylo

Approach: look at >50bp intronic insertions.

HG38 INTRONS: hg38, Genes, GENCODE v24, BED format, introns only.
--> introns.sorted.bed
951 423 lines

sed -i 's/chr23/chrX/' introns.sorted.bed
sed -i 's/chr24/chrY/' introns.sorted.bed

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