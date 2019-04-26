#!/bin/bash
#SBATCH -o sb_costea2017-motus2_%j.out
#SBATCH -c 4

# Replace N with number of processors / cores you want to use
# If submitting to a cluster with sbatch, also adjust in "SBATCH -c N" above
nproc=4

# Fill in the paths
reads_path=/home/mrmclare/data/mgs_bias/costea2017/reads
out_path=/home/mrmclare/data/mgs_bias/costea2017/motus2

# The ERA run accession, used to create the fastq.gz file names; read in as command line arg
accession=$1

eval "$(/home/mrmclare/applications/miniconda3/bin/conda shell.zsh hook)"
conda activate motus

motus profile \
    -f $reads_path/${accession}_1.fastq.gz \
    -r $reads_path/${accession}_2.fastq.gz \
    -o $out_path/${accession}-counts.motus \
    -M $out_path/${accession}.mgc \
    -n $accession \
    -c \
    -g 1 \
    -y insert.raw_counts \
    -t $nproc

# Afterwards, merge all with
#> FLIST=$(echo ERR*-counts.motus | tr ' ' ,)
#> motus merge -i $FLIST -o costea2017-phase3-counts.motus
