# Metagenomic-virus host prediction

## Linux workflow for predicting hosts of viruses found using metagenomics.

# Background
- Viruses found in metagenomic datasets are often not linked to their hosts due to high biodiversity in samples (many potential pairings = difficult to infer host)
- Despite this, host nucleic acids are often sequenced alongside their viruses; we must therefore find the needle in the haystack
- Given a cohort of samples, some containing a virus lineage, and some not, it may be possible to predict hosts
- This workflow can be used to identify viruses & cellular taxa in samples, and then pinpoint which possible host taxa occur at high prevalence alongisde viruses
- Follow-up with statistical testing can establish whether these taxa are found at higher prevalence than in virus negative samples
- Once computational prediction of a host is done, targeted validation is possible
- The workflow assumes a eukaryotic host, and otherwise remains agnostic to host taxonomy. Though untested, it should also work for prokaryotes with a minor change

# Dependencies:
Prerequisite programmes in $PATH:
- bwa
- blast
- bbmap
- pathoscope
- entrez-direct
- samtools

Alternatively, a conda environment.yml is provided in this repository for fast setup.
Get conda here:
https://docs.conda.io/en/latest/miniconda.html#linux-installers

Create, activate, and check environment:
```
conda hostPredict create -f environment.yml
conda activate hostPredict
conda list
```

# Required input files:
- Metagenomic sequence reads from a suitable cohort. These are assumed to be paired-end fastq.gz, with extensions "_1.fastq.gz" & "_2.fastq.gz"
- Virus reference genomes in a file named "viruses.fas". Remove spaces & special characters from headers, e.g.:
```
>Kirkoviridae_7_sampled_pig_stool
>Kirkoviridae_8_sampled_bovine_stool
>Kirkoviridae_9_sampled_pig_stool
>Redondoviridae_100_sampled_human_oral
>Redondoviridae_101_sampled_human_periodontal_tissue
>Redondoviridae_102_sampled_human_gut
```
- The SILVA SSU and LSU NR99 databases, merged. SILVA databases are available here: https://www.arb-silva.de/index.php?id=272. For version 138.1:
```
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138.1/Exports/SILVA_138.1_LSURef_NR99_tax_silva.fasta.gz
zcat SILVA_138.1_* | sed 's/ /_/g' | sed 's/;/_/g' | sed 's/\//_/g' | sed 's/://g' | sed '/^[^>]/ y/uU/tT/' | gzip > SILVA_138.1_LSU_SSU_Ref_NR99_tax_silva.fasta.gz
rm SILVA_138.1_?SURef_NR99_tax_silva.fasta.gz
```
- GenBank nt v5 database. To download, use the script "update_blastdb.pl" provided with BLAST+ (N.B. database is >100 GB)
- Ensure taxonomy is properly integrated into the database (see below)
- Due to long download time, this should be run on a compute node (most login nodes will terminate the process before it can finish)
```
cd /desired_database_location
perl /BLAST_installation_directory/bin/update_blastdb.pl --blastdb_version 5 nt --decompress
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
tar -xzf taxdb.tar.gz
```
- To ensure the taxonomy file you downloaded can be seen by BLAST, you need to edit the ~/.bashrc file (so this will be run on all future ssh client sessions)
```
nano ~/.bashrc
# Paste following two lines:
BLASTDB=/desired_database_location
export BLASTDB

# Exit the file (ctrl + X), save (hit Y), retain file name (hit Enter). 
# Restart ssh client session
```








# Usage:

- Steps 1 and 2 are provided as slurm batch scripts. Adjust workload manager information and variables for your system
- Step 3 is a short bash script to be run after steps 1 and 2 are finished, and a list of virus positive samples is known (see below)

## Step 1: Identify virus positive samples (GenerateVirusMatrix.sh)
- Outputs a data matrix and processed sam files
- First two matrix columns contain sample names & read counts for later normalisation
- Remaining data columns (one per reference genome) contain raw virus read counts after quality control
- Counts are merged per sample (i.e. forward and reverse read files = 1 data row with total values)

## Step 2: Metagenomic analysis of eukaryotic taxa in samples
(name.sh)
- Outputs...
- 

## Step 3: Identify potential host taxa enriched in virus positive samples 
(NAME.sh)

- Manually inspect virus output matrix - apply any cutoffs for calling a sample positive for a virus lineage
- Generate a text file, **positives.txt**, containing a list of samples positive for your virus lineage of interest, e.g.: 
```
SRR6713816
SRR6713817
SRR6713818
````
- **From the output directory** of steps 1 and 2, run step 3, provided as a bash script
- Outputs 

