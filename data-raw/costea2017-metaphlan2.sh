#!/bin/bash
#SBATCH -o sb_metaphlan2_%j.out
#SBATCH -c N

# NOTE: make sure that bowtie2 is in the path

# Replace N with number of processors / cores you want to use
# If submitting to a cluster with sbatch, also adjust in "SBATCH -c N" above
nproc=N

# Fill in the paths
reads_path=/folder/containing/reads
out_path=/folder/for/metaphlan2/output
metaphlan2_path=/path/to/metaphlan2.py

# The ERA run accession, used to create the fastq.gz file names; read in as command line arg
accession=$1

$metaphlan2_path \
    --bowtie2out $out_path/${accession}_bowtie2.bz2 \
    --nproc $nproc --input_type fastq \
    $reads_path/${accession}_1.fastq.gz,$reads_path/${accession}_2.fastq.gz \
    $out_path/${accession}_profiled_metagenome.tsv
