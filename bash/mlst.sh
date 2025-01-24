#!/bin/bash
# MLST

# --- write job start time
echo "Job started at $(date)"

#activate conda env
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate mlst

OUT=/your/path/here/3.2-sero/mlst

cd /your/path/here/3.0-unicycler_entero/assemblies

mlst --scheme senterica_achtman_2 *.fasta > ${OUT}/mlst.tsv

# --- write job end time
echo "Job ended at $(date)"