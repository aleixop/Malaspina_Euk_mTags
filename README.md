# A metagenomic assessment of microbial eukaryotic diversity in the global ocean

This repository contains code and data included in:

--------
**Obiol A., Giner C.R., Sánchez P., Duarte C.M., Acinas S.G., Massana R.** (2020). *A metagenomic assessment of microbial eukaryotic diversity in the global ocean*.

--------

## Scripts

### mTags extraction and classification pipeline

  - `scripts/mtags_pipeline/0-mtags_extraction.sh`: Bash script for mTags extraction.
  - `scripts/mtags_pipeline/1-mtags_classification.sh`: Bash script for mTags classification. 
  - `scripts/mtags_pipeline/mtags_consensus_tax_v2.py`: Python function to make consensus taxonomy from the above script's output.
  - `scripts/mtags_pipeline/make_otu_table.py`: auxiliary Python function to build an OTU table out from the above script's output.

### Data processing and analysis (R scripts)

  - `scripts/data_analysis/1.1-mtags_piconano_taxlevels.R`: assessment of mTags taxonomic resolution.
  - `scripts/data_analysis/1.2-mtags_piconano_diversity.R`: diversity as seen by mTags in pico and nano fraction.
  - `scripts/data_analysis/2-comparison_amplicons.R`: comparison between mTags and V4/V9 amplicons.
  - `scripts/data_analysis/3-hmm_vs_nohmm.R`: comparison of mTags extraction methods (BLAST vs HMM).
  - `scripts/data_analysis/4-150_vs_100bp.R`: comparison of taxonomic resolution between 100 bp and 150 bp long mTags. 
  - `scripts/data_analysis/stat_smooth_func.R`: auxiliary function for scatter plots in Fig. S3. Original source [here](https://gist.github.com/kdauria/524eade46135f6348140).

## Tables

  - `data/tables/MPN_VP_metaG_OTUtable.txt`: mTags table.
  - `data/tables/MPN_VP_ampliconsV4_ASVtable.txt`: amplicons table for V4 region.
  - `data/tables/MPN_VP_ampliconsV9_ASVtable.txt`: amplicons table for V9 region.
  - `data/tables/MPN_hmm_comparison_table.txt`: mTags table with both BLAST and HMM extraction methods in order to compare them.
  - `data/tables/MPN_amplicons_mtags_comparison_table.txt`: amplicons and mTags table with shared samples between the 3 datasets.

## mTags sequences

  - `data/mtags/MPN_VP_mTags.fasta.gz`: original fasta file with all retrevied mTags.
  - `data/mtags/MPN_VP_mTags_101bp.fasta.gz`: fasta file with all retrieved mTags trimmed to 101bp. 

## Contigs sequences

  - `data/contigs/MPN_VP_contigs.xlsx`: Excel file with all retrieved contigs and their associated information.
  - `data/contigs/MPN_VP_contigs.fasta.gz`: fasta file with all retrieved contigs.
  - `data/contigs/diplonemea_tree_alignment.phy`: alignment of Diplonemea contigs used to build Fig. 5.

## Metadata

  - `data/metadata/MPN_VP_metadata.txt`: metadata for each sample code.

## *eukaryotesV4* database

*eukaryotesV4*, the database built for this study, is available [in the following repository](https://github.com/aleixop/eukaryotesV4).

In order to extract 18S-V4 mTags, a 92% clustered version of *eukaryotesV4* was used. It can be found in `data/db/eukaryotesV4_92clust.fasta.gz`.

## Raw data

  - Amplicons raw data can be found at European Nucleotide Archive under accession number [PRJEB23771](https://www.ebi.ac.uk/ena/data/view/PRJEB23771).
