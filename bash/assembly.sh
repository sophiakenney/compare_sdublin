#!/bin/bash
# Genome Assembly

# --- write job start time
echo "Job started at $(date)"

# --- activate conda env
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate unicycler

# --- set paths up for unicycler
FILT="/your/path/here/2.1-k2filteredentero"
UNI=/your/path/here/Unicycler
OUT=/your/path/here/3.0-unicycler_entero

cd ${FILT} #change to directory with filtered reads

for f in *.filt_1P.fq;
do

echo "assemble" ${f%.filt_1P.fq}

${UNI}/unicycler-runner.py -1 ${FILT}/${f} \
-2 ${FILT}/${f%.filt_1P.fq}.filt_2P.fq \
-o ${OUT}/${f%.filt_1P.fq} ;

done

mkdir ${OUT}/assemblies #make directory for all assemblies

for f in *.filt_1P.fq;
do

cd ${OUT}/${f%.filt_1P.fq}
cat assembly.fasta > ${f%.filt_1P.fq}_assembly.fasta
mv ${f%.filt_1P.fq}_assembly.fasta ../assemblies
cd ..;

done

# --- write job start end
echo "Job started at $(date)"