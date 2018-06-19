# Intronic Insertions for Phylogenomics

This README includes basic instructions for how to use the code contained in the scripts folder. While the pipeline is ultimately intended for extension to other species, there are currently a number of locations in the code that are human/dog/mouse specific. 

The aim is to build this pipeline such that given n species, evidence for all permutations of insertion/deletion amongst those species (excepting outgroup, since we only consider evidence that is deleted in outgroup) will be generated in an automated fashion.

The following sections describe the stages in the pipeline. 

## Step 1: Label introns with chainIDs

The first step is to label the introns file with the chainIDs of the species in which we wish to find insertions or deletions. These chainIDs, and the chain files themselves, all refer to chains mapping between the chosen reference species and each respective query species (e.g. hg38.mm10.all.chain, in which case hg39 is the reference and mm10 is the query). 

The file 
```
/cluster/u/jschull/phylo/sorted/introns.byID.bed
```
contains all reference (hg38) introns, labelled and sorted by their Ensembl IDs.

To label this file, we use a simple join. For example, to label it with canFam3 chain IDs, we use:

```
join -1 4 -2 1 introns.byID.bed ../chains/canFam3/canFam3.trimmed.chain_ids > canFam3/introns.canFam3.bed
```

## Step 2: Find query1 deletions/insertions

This step of the pipeline generates files containing sites of the form (ref +, query -) or (ref -, query +). We refer to the first case as a *deletion*, and the second as an *insertion*. 

###### Input 

To generate these files, we use the script **get_indels.py**, which takes arguments of the following format:
```
python get_indels.py \[SITESFILE\] \[CHAINFILE\] \[OUTFILE\] \[TYPE\] \[START\] \[END\]
```
- SITESFILE is a file containing the BED coordinates of sites within which we would like to search (e.g. introns.canFam3.bed)
- CHAINFILE is a chain file for the query we are searching (e.g. hg38.canFam3.all.chain)
- OUTFILE is the name of the file (including path) that we would like to write results to.
- TYPE indicates whether we are searching for deletions or insertions. TYPE is either 0 or 1: 0 indicates deletions, and 1 indicates insertions.
- START is the first line in SITESFILE that we should look at. This is a parallelization parameter, allowing us to look at a select number of lines in the SITESFILE.
- END is an optional argument indicating the last line in SITESFILE that we should look at. If left absent, all lines after START will be looked at.

###### Output

This will output a file of the format:

ENSEMBL_ID, INSERTION SIZE, REF CHROM, REF START, REF END, QUERY CHROM, QUERY STRAND, QUERY START, QUERY END


## Step 3: Generate evidence

The final step is to use the files containing insertions/deletions to compare to all the *other* queries, and generate the final evidence. This step uses the scripts **create_job_list.py** to quickly generate jobfiles for submission to parasol, and **generate_evidence.py** to actually compare the sites to the relevant queries.

The script **create_job_list.py** takes arguments of the following format:
```
python create_job_list.py \[SCRIPT-TYPE\] \[SITESFILE\] \[CHAINFILE\] \[OUTFILE\] \[JOBFILE\] \[TYPE-OF-SEARCH\] \[LINES-PER-JOB\]
```




