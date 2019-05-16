#!/bin/bash
#SBATCH -o sb-costea2017-metaphlan2-%j.out
set -euo pipefail

# NOTE: make sure that bowtie2 is in the path

# The ERA run accession, used to create the fastq.gz file names
accession=$1
# Number of cpus to use
nproc=$2

# Set paths based on the DATA_PATH set in the .env file
export $(cat .env | xargs)
out_path=$DATA_PATH/costea2017/metaphlan2

metaphlan2.py \
    --input_type bowtie2out \
    --nproc $nproc \
    --min_cu_len 0 \
    --stat avg_g \
    $out_path/${accession}_bowtie2.bz2 \
    $out_path/${accession}_profiled_metagenome_custom.tsv

# Submit all accessions to Slurm
#
# nproc=2
# for i in {1971003..1971031}
# do
#     accession=ERR$i
#     # bash costea2017-metaphlan2-custom.sh $accession $nproc
#     sbatch -c $nproc costea2017-metaphlan2-custom.sh $accession $nproc
# done

# Afterwards, merge all accessions into a single table with
# /path/to/merge_metaphlan_tables.py *_profiled_metagenome_custom.tsv > costea2017_metaphlan2_custom_profiles.tsv
