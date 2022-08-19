# Metagenomic-virus host prediction

## Linux workflow for predicting hosts of viruses found using metagenomics.

# Background
- Viruses found in metagenomic datasets are often not linked to their hosts due to high biodiversity in samples (many potential pairings = difficult to infer host)
- Despite this, host nucleic acids are often sequenced alongside their viruses; we must therefore find the needle in the haystack
- Given a cohort of samples, some containing a virus lineage, and some not, it may be possible to predict hosts
- This workflow can be used to identify viruses & cellular taxa in samples, and then pinpoint which possible host taxa occur at high prevalence alongisde viruses
- Follow-up with statistical testing can establish whether these taxa are found at higher prevalence than in virus negative samples
- Once computational prediction of a host is done, targeted validation is possible
- The workflow assumes a eukaryotic host, and otherwise remains agnostic to host taxonomy. Though untested, it should also work for prokaryotes with minor changes

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
- Virus reference genomes in a file named "viruses.fas". Remove spaces or special characters from headers, e.g.:
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
- GenBank nt v5 database. To download, use the script provided with BLAST+ (N.B. database is >100 GB)
```
update_blastdb.pl --blastdb_version 5 nt --decompress
```

# Usage:

## Step 1: Identify virus positive samples

## Step 2: Metagenomic analysis of cellular taxa in samples

## Step 3: Identify potential host taxa enriched in virus positive samples


# Notes:



