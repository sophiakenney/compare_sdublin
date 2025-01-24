#!/bin/bash
#Abricate for plasmids

# --- write start time
echo "Job started at $(date)"

# --- activate conda env
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate abricate

#set directory paths
ASSEMB=/your/path/here/3.0-unicycler_entero/assemblies
PF=/your/path/here/4.2-abricate

#change to working directory
cd ${ASSEMB}


abricate --db plasmidfinder *_assembly.fasta > ${PF}/enteroplasmid-all.tab

#summarize into presence absence matrix
abricate --summary ${PF}/enteroplasmid-all.tab > ${PF}/enteroplasmid-all_summary.txt

echo "PF ended at $(date)"

# --- write end time
echo "Job started at $(date)"