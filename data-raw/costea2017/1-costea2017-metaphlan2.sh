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
reads_path=$DATA_PATH/costea2017/reads
out_path=$DATA_PATH/costea2017/metaphlan2

metaphlan2.py \
    --bowtie2out $out_path/${accession}_bowtie2.bz2 \
    --nproc $nproc --input_type fastq \
    $reads_path/${accession}_1.fastq.gz,$reads_path/${accession}_2.fastq.gz \
    $out_path/${accession}_profiled_metagenome.tsv

# Submit all accessions to Slurm
#
# nproc=8
# for i in {1971003..1971031}
# do
#     accession=ERR$i
#     sbatch -c $nproc costea2017-metaphlan2.sh $accession $nproc
# done

# Afterwards, merge all accessions into a single table with
# /path/to/merge_metaphlan_tables.py *_profiled_metagenome.tsv > costea2017_metaphlan2_profiles.tsv
