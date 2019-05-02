#!/bin/bash
#SBATCH -o sb_costea2017-motus2-map_snv_%j.out
#SBATCH -c 4
set -euo pipefail

# Replace N with number of processors / cores you want to use
# If submitting to a cluster with sbatch, also adjust in "SBATCH -c N" above
nproc=4

# Set paths based on the DATA_PATH set in the .env file
export $(cat .env | xargs)
reads_path=$DATA_PATH/costea2017/reads
out_path=$DATA_PATH/costea2017/motus2/bam

# The ERA run accession, used to create the fastq.gz file names; read in as command line arg
accession=$1

eval "$(/home/mrmclare/applications/miniconda3/bin/conda shell.zsh hook)"
conda activate motus

# TODO: Create out path if it doesn't exist; edit conda activation command to use .env file

# Align reads
motus map_snv \
    -f $reads_path/${accession}_1.fastq.gz \
    -r $reads_path/${accession}_2.fastq.gz \
    -t $nproc \
    > $out_path/${accession}.bam
