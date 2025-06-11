#!/bin/bash
# AMRFinderPlus

# --- write job start time
echo "Job started at $(date)"

#activate conda env
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate amrfinder

#set paths
PROT=/your/path/here/4.0-prokka/faa
ANNO=/your/path/here/4.0-prokka/gff
ASSEMB=/your/path/here/3.0-unicycler_entero/assemblies
OUT=/your/path/here/4.3-amrfinder

# Run

cd ${ASSEMB} #change to directory with assemblies

for f in *_assembly.fasta;
do

echo "running" ${f%_assembly.fasta}

amrfinder -p ${PROT}/${f%_assembly.fasta}.faa \
-n ${f} \
-g ${ANNO}/${f%_assembly.fasta}.gff \
-a prokka \
-O Salmonella --plus \
-o ${OUT}/${f%_assembly.fasta} ;

done

echo "Job ended at $(date)"