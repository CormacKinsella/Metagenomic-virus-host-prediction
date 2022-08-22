#!/bin/bash
#SBATCH -N 1
#SBATCH -t 120:00:00
#SBATCH --mem=60G
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=USER@example.com

#########################################################################
# SET THESE PARAMETERS (N.B. $TMPDIR variable refers to compute node scratch space, adjust this to your system)

# Activate conda environment: conda activate hostPredict

# Data. Assumes paired-end fastq.gz, with extensions "_1.fastq.gz" & "_2.fastq.gz"
cp /directory/*fastq.gz $TMPDIR

# SILVA SSU and LSU database
cp /directory/SILVA_138.1_LSU_SSU_Ref_NR99_tax_silva.fasta.gz $TMPDIR
silva=$TMPDIR/SILVA_138.1_LSU_SSU_Ref_NR99_tax_silva.fasta.gz

# GenBank v5 nt database
cp /directory/nt* $TMPDIR
cp /directory/taxdb.b?? $TMPDIR
BLASTDB=$TMPDIR
export BLASTDB

# Set taxonomy Rank (the taxonomic rank used for over-representation search)
# E.g. class; order; family; genus
taxRank=genus

# Out directory
outdir=/directory/outputs

# Set number of threads to use
thread_num=8

# Set wd to compute node
cd $TMPDIR

#########################################################################

# Start script

# Index reference
bwa index $silva

# Set up a while loop for analysis (ensures partial output even if job times out)
# Identify samples to analyse
ls -1 *_1.fastq.gz | sed 's/_1.fastq.*//' > sample_list.txt

while read line
do

#  Map reads to rRNA db
for i in $line*fastq.gz; do bwa mem -t $thread_num $silva $i | samtools view -S -F4 - > ${i}.rRNA.sam; done

# Reassign multi-mapped reads to best hit in sam file
for i in *sam; do touch updated_$i; done
ls -1 *sam | grep -v updated_ > process_files
while read line; do pathoscope ID -alignFile $line -fileType sam -outDir . ; done < process_files; rm pathoid-sam-report.tsv process_files
for i in updated_*; do mv $i ${i#updated_}; done

# Exclude prokaryotic hits from sam files (remove "-v" to exclude eukaryotes)
for i in *sam; do grep -v "Bacteria\|Archaea" $i > ${i%sam}subset.sam; done

# Reformat sams
for i in *subset.sam; do reformat.sh in=$i out=${i%sam}fasta; done

# Run BLAST against nt database
for i in *.subset.fasta; do blastn -query $i -db nt -outfmt "6 qseqid sacc pident length qlen mismatch gapopen sstart send evalue bitscore staxid ssciname stitle" -num_threads $thread_num | sort -k1,1 -k11,11nr -k10,10n | sort -u -k1,1 --merge > ${i}.fmt6; done

# Merge sample files (done here to ensure no upstream loss of identical read names in F and R files)
for i in *_1.fastq.gz*fmt6; do cat $i ${i%1.fastq.gz*}2.fastq.gz.rRNA.subset.fasta.fmt6 >> ${i%1.fastq.gz*}merged.fastq.gz.rRNA.subset.fasta.fmt6; done

# Filter only high quality alignments (100% identity for 100 bp), take out uniq taxids
for i in *merged*fmt6; do awk '$3 == 100' $i | awk '$4 > 99' | cut -f12 | sort | uniq > ${i}.hq.txids; done

# Taxonomy lookup on txids, extract selected rank
for i in *txids; do touch ${i}.taxList; done
for i in *txids; do while read line; do efetch -db taxonomy -id $line -format xml | xtract -pattern Taxon -block "*/Taxon" -unless Rank -equals "no rank" -tab "\n" -element Rank,ScientificName >> ${i}.taxList; done < $i; done
for i in *taxList; do sort $i | uniq > ${i}.sort; done; for i in *sort; do mv $i ${i%.sort}; done
for i in *taxList; do grep "$taxRank" $i | cut -f2 > ${i}.$taxRank.uniq; done

# Final clear-up
gzip *rRNA.sam
mv *taxList* *merge*fmt6 *rRNA.sam.gz $outdir
rm *txids *subset*
done < sample_list.txt
