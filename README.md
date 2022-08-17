# Metagenomic-virus host prediction

## Linux workflow for identifying over-represented eukaryotic taxa in virus positive samples.

# Highlights
- Viruses found in metagenomic datasets are often not linked to their hosts due to high biodiversity in samples (many potential pairings = difficult to infer host)
- Despite this, hosts nucleic acids are often sequenced alongside their viruses
- Given a cohort of samples, some containing the host, and some not - it may be possible to predict likely hosts
- A set of virus positive samples should contain an enrichment in the eukaryotic/prokaryotic host versus virus negative samples
- Follow-up tests can then be targeted, given a specific putative host
- This workflow carries out the necessary sub-steps
- It remains agnostic to host taxonomy (though you may target eukaryotic vs. prokaryotic hosts for example) 

# Usage:

Prerequisite programmes in $PATH:
- xyz

Alternatively, a conda environment.yml is provided in this repository for fast setup.
Get conda here:
https://docs.conda.io/en/latest/miniconda.html#linux-installers

Create, activate, and check environment:
```
conda hostPredict create -f environment.yml
conda activate hostPredict
conda list
```

## Step 1: Identify virus positive samples

## Step 2: Metagenomic analysis of cellular taxa in samples

## Step 3: Identify potential host taxa enriched in virus positive samples

# Required input files:
- SILVA SSU and LSU NR99 databases, merged (https://www.arb-silva.de/download/arb-files/). 
- Local GenBank nt database for BLAST searches
- fasta or fastq files derived from cohort metagenomic sequencing


# Notes:



