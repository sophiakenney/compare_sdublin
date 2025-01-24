#!/bin/bash
# IQTree

#--- write job start time
echo "Job started at $(date)"

#--- activate env
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate iqtree

#--- set paths 
OUT=/your/path/here/5.1-iqtree
ALN=/your/path/here/4.1-roary

cd ${OUT}

# --- run without pathogen
iqtree -s ${ALN}/core_gene_alignment.aln \
       -m GTR+I+G4 \
       --seed 1234 \
       --verbose \
       --ufboot 1000 \
       -bnni \
       -nt AUTO \
       --prefix iq1


#--- write job end time
echo "Job ended at $(date)"