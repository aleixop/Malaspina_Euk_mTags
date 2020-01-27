#!/bin/sh

# mTags classification

# Load modules

module load usearch/9.2.64
module load python/3.6.5

# Files and paths

FASTA=<fasta_with_extracted_mtags>
MAP=<output_tophits_map>
DB=<database>
OUT_TABLE=<out_dir_for_table>
OUT_MAP=<out_dir_for_map>
OUT_NAME=<name_of_final_file>
MITAGS_CONSENSUS_TAXONOMY=mtags_consensus_tax_v2.py # python script to filter map
MAKE_OTU_TABLE=make_otu_table.py # python script to create OTU table

# Top hits map

usearch \
-usearch_local ${FASTA} \
-db ${DB} \
-uc ${MAP} \
-id 0.97 \
-mincols 70 \
-strand both \
-top_hits_only \
-maxaccepts 0 \
-maxrejects 0 \
-threads 24

# Make consensus taxonomy

${MITAGS_CONSENSUS_TAXONOMY} \
  --tax_separator '_' \
  --tax_sense 'asc' \
  --pair_separator '/' \
  --output_file ${OUT_MAP}/${OUT_NAME}.uc \
  ${MAP}

# Make otu table

${MAKE_OTU_TABLE} \
  --sample_identifier 'barcodelabel=' \
  --output_file ${OUT_TABLE}/${OUT_NAME}_otuTable.txt \
  ${OUT_MAP}/${OUT_NAME}.uc
