#!/bin/bash
#ModelTest for Phylogeny


#--- write job start time
echo "Job started at $(date)"

#--- set paths
ROARY=/your/path/here/4.1-roary
OUT=/your/path/here/5.0-modeltest

#--- run program
/syour/path/here/modeltest/bin/modeltest-ng -i ${ROARY}/core_gene_alignment.aln \
-o ${OUT} \
-t ml \
-d nt \
-T raxml \
-p 8

#--- write job end time
echo "Job ended at $(date)"
