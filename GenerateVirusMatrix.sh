#!/bin/bash
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH --mem=60G
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=USER@example.com

#########################################################################
# SET THESE PARAMETERS (N.B. $TMPDIR variable refers to compute node scratch space, adjust this to your system)

# Activate conda environment: conda activate hostPredict

# Data. Assumes paired-end fastq.gz, with extensions "_1.fastq.gz" & "_2.fastq.gz"
cp /directory/*fastq.gz $TMPDIR

# Mapping reference transfer
cp /directory/viruses.fas $TMPDIR
reference_file=$TMPDIR/viruses.fas

# Out directory
outdir=/directory/outputs

# Set number of threads to use
thread_num=8

# Set wd to compute node
cd $TMPDIR

#########################################################################

# Start script

# Index reference
bwa index $reference_file
makeblastdb -in $reference_file -dbtype nucl

# Store sample names and read counts for later normalisation
ls -1 *_1.fastq.gz | sed 's/_1.fastq.gz//' > names
while read line; do count_forward=`echo $(zcat ${line}_1.fastq.gz | wc -l)/4|bc`; count_reverse=`echo $(zcat ${line}_2.fastq.gz | wc -l)/4|bc`; echo $(($count_forward + $count_reverse)) >> counts; done < names
paste names counts > temp_data; rm names counts
echo -e 'sample_name\tread_count' > header
cat header temp_data > sample_info; rm header temp_data

# Map reads to reference
for i in *fastq.gz; do bwa mem -t $thread_num $reference_file $i | samtools view -S -F4 - > ${i}.vir.sam; done

# Reassign multi-mapped reads to best hit in sam file
for i in *sam; do touch updated_$i; done
ls -1 *sam | grep -v updated_ > process_files
while read line; do pathoscope ID -alignFile $line -fileType sam -outDir . ; done < process_files; rm pathoid-sam-report.tsv process_files
for i in updated_*; do mv $i ${i#updated_}; done

# QC of reads
for i in *sam; do reformat.sh in=$i out=${i}.fasta; done
for i in *sam.fasta; do blastn -query $i -db $reference_file -outfmt "6 qseqid sacc pident length qlen mismatch gapopen sstart send evalue bitscore staxid stitle" -num_threads $thread_num -word_size 11 -gapopen 5 -gapextend 2 -penalty -3 -reward 2 -dust yes | sort -k1,1 -k11,11nr -k10,10n | sort -u -k1,1 --merge > ${i}.fmt6; done
for i in *fmt6; do awk '$4 > 40' $i | cut -f1 | sort | uniq > $i.hq; done
for i in *sam; do grep -Ff ${i}.fasta.fmt6.hq $i > ${i%sam}hq.sam; done

# Merge paired-end outputs (done late, since identical forward and reverse read names will otherwise be lost)
for i in *_1.fastq.gz.vir.hq.sam; do cat $i ${i%1.fastq.gz*}2.fastq.gz.vir.hq.sam >> ${i%1.fastq.gz*}merged.fastq.gz.vir.hq.sam; done

# Construct data matrix
grep ">" $reference_file | sed 's/>//' > header_lines
while read line; do for i in *merged*hq.sam; do cut -f1-3 $i | grep $line | wc -l >> $line.col; done; done < header_lines
temp_outs=`sed 's/$/.col/' header_lines | tr '\n' '\t'`
paste $temp_outs > temp_matrix; rm *.col
tr '\n' '\t' < header_lines > header; rm header_lines; echo "" >> header
cat header temp_matrix > named_temp_matrix
paste sample_info named_temp_matrix > virus_matrix.txt; rm sample_info *temp_matrix header

# Clearup
mv *merged*hq.sam virus_matrix.txt $outdir
