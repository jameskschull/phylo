# Intronic Insertions for Phylogenomics

This README includes basic instructions for how to use the code contained in the scripts folder. While the pipeline is ultimately intended for extension to other species, there are currently a number of locations in the code that are human/dog/mouse specific. 

The aim is to build this pipeline such that given n species, evidence for all permutations of insertion/deletion amongst those species (excepting outgroup, since we only consider evidence that is deleted in outgroup) will be generated in an automated fashion.

The following sections describe the stages in the pipeline. 

## Directories



## Pipeline

### Step 1: Label introns with chainIDs

The first step is to label the introns file with the chainIDs of the species in which we wish to find insertions or deletions. These chainIDs, and the chain files themselves, all refer to chains mapping between the chosen reference species and each respective query species (e.g. hg38.mm10.all.chain, in which case hg38 is the reference and mm10 is the query). 

The file 
```
/cluster/u/jschull/phylo/sorted/introns.byID.bed
```
contains all reference (hg38) introns, labelled and sorted by their Ensembl IDs.

To label this file, we use a simple join. For example, to label it with canFam3 chain IDs, we use:

```
join -1 4 -2 1 introns.byID.bed ../chains/canFam3/canFam3.trimmed.chain_ids > canFam3/introns.canFam3.bed
```

### Step 2: Find query1 deletions/insertions

This step of the pipeline generates files containing sites of the form (ref +, query -) or (ref -, query +). We refer to the first case as a *deletion*, and the second as an *insertion*. To generate these files, we use the script **get_indels.py**; however, in order to parallelize the process with parasol, we use an intermediate script **create_job_list.py** to create a list of jobs that call get_indels.py.

##### Creating a job list

The script **create_job_list.py** takes arguments of the following format:
```
python create_job_list.py [SCRIPT-TYPE] [SITESFILE] [CHAINFILE] [OUTFILE] [JOBFILE] [TYPE-OF-SEARCH] [LINES-PER-JOB]
```
- SCRIPT-TYPE is either 'indels' or 'evidence': the former indicates that the jobs will be calling **get_indels.py**, whereas the latter indicates that the jobs will be calling **generate_evidence.py**. 
- SITESFILE is the name of the file containing the candidate sites that we wish to compare to the relevant query. Note that this file needs to be labelled with the chain IDs of the query being compared to. 
- CHAINFILE is the name of the chain file for the query we are searching.
- OUTFILE is the name of the file that we would like to write results to.
- JOBFILE is the name of the file to which each individual job will be written. 
- TYPE-OF-SEARCH differs depending on the script type. For **get_indels.py** jobs, TYPE-OF-SEARCH will be either 0 or 1. For **generate_evidence.py** jobs, TYPE-OF-SEARCH will be one of the four search types described in the section on **generate_evidence.py**.
- LINES-PER-JOB is the number of lines in SITESFILE that you would like each job to look at. 

##### Input 

The script **get_indels.py** takes arguments of the following format:
```
python get_indels.py [SITESFILE] [CHAINFILE] [OUTFILE] [TYPE] [START] [END]
```
- SITESFILE is a file containing the BED coordinates of sites within which we would like to search (e.g. introns.canFam3.bed)
- CHAINFILE is a chain file for the query we are searching (e.g. hg38.canFam3.all.chain)
- OUTFILE is the name of the file (including path) that we would like to write results to.
- TYPE indicates whether we are searching for deletions or insertions. TYPE is either 0 or 1: 0 indicates deletions, and 1 indicates insertions.
- START is the first line in SITESFILE that we should look at. This is a parallelization parameter, allowing us to look at a select number of lines in the SITESFILE.
- END is an optional argument indicating the last line in SITESFILE that we should look at. If left absent, all lines after START will be looked at.

Note, however, that these will be filled in by **create_job_list.py**.

Running these jobs will create files of the following format (which can be concatenated together to generate the final indels file):

ENSEMBL_ID, INSERTION SIZE, REF CHROM, REF START, REF END, QUERY CHROM, QUERY STRAND, QUERY START, QUERY END


### Step 3: Generate evidence

The final step is to use the files containing insertions/deletions to compare to all the *other* queries, and generate the final evidence. This step uses the scripts **create_job_list.py** to quickly generate jobfiles for submission to parasol, and **generate_evidence.py** to actually compare the sites to the relevant queries.

##### Generating evidence

The script **generate_evidence.py** takes one of four 'mode' arguments (these are 'TYPE-OF-SEARCH' in **create_job_list.py**). These indicate the method that should be used to compare from queryX to queryY, and are differentiated as follows:

1. 'di': Reference deletion, search for insertion in query (e.g. hg38-.mm10+.canFam3+.txt).
2. 'dd': Reference deletion, search for deletion in query (e.g. hg38-.mm10+.canFam3+.loxAfr3-.txt).
3. 'ii': Reference insertion, search for insertion in query (e.g. hg38+.mm10-.canFam3+.txt).
4. 'id': Reference insertion, search for deletion in query (e.g. hg38+.mm10-.canFam3+.loxAfr3-txt).

The other arguments are also provided by **create_job_list.py**.

The current methods for each mode are as follows:

##### DI

We use the coordinates of the queryX insertion (in reference and queryX) to extract the queryX and queryY sequences. We compare them using edlib, and consider them a shared insertion if they have above INS_THRESHOLD similarity, currently 0.7.

##### DD

We navigate through the queryY chain file; we require that the reference coordinate of the queryX insertion--plus a small MARGIN (currently 5bp on either side)--is contained in a gapless block.

##### II

To remain consistent, we also use edlib comparison for this method, extracting the reference and queryY sequences by using the coordinates of the insertion. We consider them a shared insertion if they have above INS_THRESHOLD similarity.

##### ID

We search the queryY chain to find a gap in the query sequence whose reference start and end coordinate are each within a MARGIN range (on either side) of the reference *insertion* start and end coordinate. We require also that this is a one-sided gap; if these conditions are met, we consider the site evidence. 


