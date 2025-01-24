#!/bin/bash
# Genome Annotation - Prokka

# --- write job start time
echo "Job started at $(date)"

#activate conda env
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate prokka_env

#set directory paths
ASSEMB=/your/path/here/3.0-unicycler_entero/assemblies
OUT=/your/path/here/4.0-prokka

cd ${ASSEMB}

for i in *_assembly.fasta ; do

mkdir ${OUT}/${i%_*.fasta}

prokka ${i} \
--outdir ${OUT}/${i%_*.fasta} \
--prefix ${i%_*.fasta} \
--force \
--kingdom Bacteria \
--genus Salmonella \
--species enterica

done

# copy all gff files to one directory

mkdir -p ${OUT}/gff

for i in *_assembly.fasta ; do cp ${OUT}/${i%_*.fasta}/${i%_*.fasta}*.gff ${OUT}/gff/ ; done

# copy all faa files to one directory

mkdir -p ${OUT}/faa

for i in *_assembly.fasta ; do cp ${OUT}/${i%_*.fasta}/${i%_*.fasta}*.faa ${OUT}/faa/ ; done

# copy all summary files to one directory

mkdir -p ${OUT}/txt

for i in *_assembly.fasta ; do cp ${OUT}/${i%_*.fasta}/${i%_*.fasta}*.txt ${OUT}/txt/ ; done

# --- write job end time
echo "Job ended at $(date)"