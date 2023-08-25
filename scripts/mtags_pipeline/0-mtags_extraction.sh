#!/bin/sh

# mTags extraction from clean fastq.gz files

# Load modules and apps

module load usearch/9.2.64
module load blast/2.7.1
module load seqtk/1.3

# Path to files

DATA_DIR=data/clean_fastq
FASTA_OUT_DIR=data/blast/mtags
BLAST_OUT_DIR=data/blast/blast_files
SAMPLE=$(ls ${DATA_DIR} | awk -F'.' '{print $1}' | sort -u | awk "NR == ${SLURM_ARRAY_TASK_ID}")

DB_BLAST="/home/aobiol/db/blastdb/Mitags_V4_ref_v15_BLAST"

# Parameters

NUM_THREADS=4

# Script

date
echo '# Begin mTags extraction'

echo '# Unzipping'
gzip -dc ${DATA_DIR}/${SAMPLE}*.gz \
  > ${FASTA_OUT_DIR}/${SAMPLE}.fastq

seqs_original=$(echo $(cat ${FASTA_OUT_DIR}/${SAMPLE}.fastq | wc -l) / 4 | bc) # count total seqs in the clean fastq file

date
echo '# Filtering'
usearch \
  -fastq_filter ${FASTA_OUT_DIR}/${SAMPLE}.fastq \
  -fastq_minlen 70 \
  -fastaout ${FASTA_OUT_DIR}/${SAMPLE}.temp.spaces.fasta

seqs_filtered_spaces=$(grep -c '^>' ${FASTA_OUT_DIR}/${SAMPLE}.temp.spaces.fasta)

echo '# Remove spaces from fasta'
cat ${FASTA_OUT_DIR}/${SAMPLE}.temp.spaces.fasta | \
  perl -pe 's/ /\//' \
  > ${FASTA_OUT_DIR}/${SAMPLE}.temp.fasta

seqs_filtered=$(grep -c '^>' ${FASTA_OUT_DIR}/${SAMPLE}.temp.fasta)

if [ "$seqs_filtered" != "$seqs_filtered_spaces" ]; then
	echo 'WARNING: Number of sequences in fasta with and without spaces do not match'
fi

echo '# Remove fastq and fasta with spaces files'
rm ${FASTA_OUT_DIR}/${SAMPLE}.fastq
rm ${FASTA_OUT_DIR}/${SAMPLE}.temp.spaces.fasta

date
echo '# Mapping'
blastn \
 -max_target_seqs 1 \
 -db ${DB_BLAST} \
 -outfmt '6 qseqid sseqid pident length qcovhsp qstart qend sstart send evalue' \
 -perc_identity 90 \
 -qcov_hsp_perc 70 \
 -query ${FASTA_OUT_DIR}/${SAMPLE}.temp.fasta \
 -out ${BLAST_OUT_DIR}/${SAMPLE}.blast \
 -num_threads ${NUM_THREADS}

date
echo '# Blast finished'

mtags_map=$(awk '{print $1}' ${BLAST_OUT_DIR}/${SAMPLE}.blast | sort -u | wc -l)

echo '# Take hits'

awk '{print $1}' ${BLAST_OUT_DIR}/${SAMPLE}.blast | sort -u > ${BLAST_OUT_DIR}/${SAMPLE}.hits

date
echo '# Extracting sequences with seqtk'

seqtk \
subseq ${FASTA_OUT_DIR}/${SAMPLE}.temp.fasta \
${BLAST_OUT_DIR}/${SAMPLE}.hits > ${FASTA_OUT_DIR}/${SAMPLE}.nolabel.fna

mtags_fna=$(grep -c '^>' ${FASTA_OUT_DIR}/${SAMPLE}.nolabel.fna)

if [ "$mtags_map" != "$mtags_fna" ]; then
	echo 'WARNING: Total extracted sequences do not match with those in map'
fi

echo '# Remove wrapping and add barcodelabel'
cat ${FASTA_OUT_DIR}/${SAMPLE}.nolabel.fna | \
  awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' | \
  tr '\t' '\n' | \
  sed "s/>.*/&;barcodelabel=${SAMPLE};/" \
  > ${FASTA_OUT_DIR}/${SAMPLE}.mtags.fna

mtags_fna_label=$(grep -c '^>' ${FASTA_OUT_DIR}/${SAMPLE}.mtags.fna)

if [ "$mtags_fna" != "$mtags_fna_label" ]; then
	echo 'WARNING: Fastas with and without label contain different sequences'
fi

#Removing intermediate files

rm ${FASTA_OUT_DIR}/${SAMPLE}.temp.fasta
rm ${FASTA_OUT_DIR}/${SAMPLE}.nolabel.fna

#Summary

echo "SUMMARY SAMPLE ${SAMPLE}; Initial seqs: $seqs_original; Filtered seqs: $seqs_filtered; mTags retrieved: $mtags_fna_label"
