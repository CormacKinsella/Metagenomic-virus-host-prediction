# Metagenomic-virus host prediction

## Linux workflow for predicting hosts of viruses found using metagenomics.

# Background
- Viruses found in metagenomic datasets are often not linked to their hosts due to high biodiversity in samples (many potential pairings = difficult to infer host)
- Despite this, host nucleic acids are often sequenced alongside their viruses; we must therefore find the needle in the haystack
- Given a cohort of samples, some containing a virus lineage, and some not, it may be possible to predict hosts
- This workflow can be used to identify viruses & cellular taxa in samples, and then pinpoint which possible host taxa occur at high prevalence alongside viruses
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
- Ensure taxonomy is properly integrated into the database (taxdb download and unpack, see below)
- Due to long download time, this should be run on a compute node (most login nodes will terminate the process before it can finish)
```
cd /desired_database_location
perl /BLAST_installation_directory/bin/update_blastdb.pl --blastdb_version 5 nt --decompress
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz
tar -xzf taxdb.tar.gz
```
# Usage:

- Steps 1 and 2 are provided as slurm batch scripts. Adjust workload manager information and variables for your system
- Step 3 is a short bash script to be run after steps 1 and 2 are finished, and a list of virus positive samples is known (see below)

## Step 1: Identify virus positive samples (GenerateVirusMatrix.sh)
- Outputs a data matrix and processed SAM files
- First two matrix columns contain sample names & read counts for later normalisation
- Remaining data columns (one per reference genome) contain raw virus read counts after quality control
- Counts are merged per sample (i.e. forward and reverse read files = 1 data row with total values)

## Step 2: Metagenomic analysis of eukaryotic taxa in samples (GenerateTaxonLists.sh)
- Outputs a list of taxa found in each sample, plus additional files for generating further outputs (e.g. taxon lists at alternative taxonomic ranks, read counts of a particular taxon, etc.)
- Script carries out taxonomy lookups for each unique NCBI taxid found (N.B. BLASTn can also output a scientific name using the outfmt option "ssciname", but this will be the lowest available taxonomy rank, normally but not always the binomial. Using lookups adds a further, rather slow step, but ensures a consistent taxonomic rank (if available), and allows any rank to be selected (e.g. genus, family, etc.). Note that BLASTn is still set to output ssciname, however this is for subsequent parsing of fmt6 files to output read counts for specific genera)
- No matrix is generated (too many taxa), so for taxa of interest (see step 3), columns can be manually generated from the fmt6 outputs, e.g.:
```
for i in *merge*fmt6; do grep "Entamoeba" $i | awk '$3 == 100' | awk '$4 > 99' | wc -l; done 
```
 
## Step 3: Find potential hosts/taxa of interest (IdentifyOverRepresented.sh)
Requires you to have identified a list of virus positive samples (**virusPositives.txt**, from step 1 output, e.g.:
```
SRR6713816
SRR6713817
SRR6713818
````
- The script is run from the output directory of step 2, where you have put **virusPositives.txt**
```
bash IdentifyOverRepresented.sh
```
- Outputs a folder containing the list of virus positive samples and a summary of which taxa are commonly found in them
- It is important to take taxa of interest identified here, and confirm if any statistical link is found between them and the viruses (i.e. are they truly over-represented in virus positive samples, or just prevalent across all samples)
