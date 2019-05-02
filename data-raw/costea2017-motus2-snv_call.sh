#!/bin/bash
#SBATCH -o sb_costea2017-motus2-snv_call_%j.out
#SBATCH -c 16
set -euo pipefail

# Replace N with number of processors / cores you want to use
# If submitting to a cluster with sbatch, also adjust in "SBATCH -c N" above
nproc=16

# Set paths based on the DATA_PATH set in the .env file
export $(cat .env | xargs)
bam_path=$DATA_PATH/costea2017/motus2/bam
out_path=$DATA_PATH/costea2017/motus2/snvs

eval "$(/home/mrmclare/applications/miniconda3/bin/conda shell.zsh hook)"
conda activate motus

# TODO: check that out path doesn't exist; edit conda activation command to use .env file

# Call SNVs
motus snv_call \
    -d $bam_path \
    -o $out_path \
    -t $nproc
