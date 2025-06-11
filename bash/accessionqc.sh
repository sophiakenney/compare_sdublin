#!/bin/bash
# Raw Read Accession and Quality Control

# --- write job start time
echo "Job started at $(date)"

# --- activate conda env
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate bioinfo

# --- download raw reads from sra

WD=/your/path/here


cd ${WD}/0.0-raw 

# make sure the accession list is in the 0.0-raw directory

for i in `cat acclist.txt`; do printf ${i}"\t"; fasterq-dump ${i} --split-files; done

# compress files 
gzip SRR*


# ---- run FastQC and MultiQC on raw reads ----

mkdir -p ${WD}/0.1-rawqc

# --- run fastqc
fastqc *.gz -o "${WD}/0.1-rawqc" -t 10

# --- multiqc

export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8

multiqc ${WD}/0.1-rawqc*_fastqc.zip --interactive

# --- run seqkit
seqkit stats *.gz -T > ${WD}/0.1-rawqc/rawstats.txt

# --- write job end time
echo "Job ended at $(date)"

