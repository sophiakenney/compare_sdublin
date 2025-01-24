#!/bin/bash
# Pangenome analysis: Roary

# --- write job start time
echo "Job started at $(date)"

#activate conda env
module load anaconda3
source activate base
conda activate roary_env

#set directory paths
GFF=/your/path/here/4.0-prokka/gff
OUT=/your/path/here/4.1-roary

cd ${GFF}

roary -e --mafft -p 10 *.gff -f ${OUT} -v

# --- write job end time
echo "Job ended at $(date)"