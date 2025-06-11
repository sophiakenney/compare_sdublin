#!/bin/bash
# Run SeqSero2

# --- write job start time
echo "Job started at $(date)"

# --- activate conda env
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate seqsero

# --- set directory paths
ASSEMB=/your/path/here/3.0-unicycler_entero/assemblies
OUT=/your/path/here/3.2-sero/seqsero

# --- change to working directory
cd ${ASSEMB}

# -m "k" for raw reads and genome assembly k-mer
# -t "4" for genome assembly as input data type
# -i path to input file
# -d output directory

for f in *.fasta ;
do

  echo ${f%_assembly.fasta}

  SeqSero2_package.py -m k -t 4 \
  -i ${ASSEMB}/${f} \
  -d ${OUT}/${f%_assembly.fasta}_seqsero \
  -n ${f%_assembly.fasta} ;

done

# create one file
cd ${OUT}

{ head -n 1 SRR###_seqsero/SeqSero_result.tsv; tail -n +2 -q *_seqsero/*_result.tsv; } > enteroseqsero_summary.tsv # you will need to edit the SRR### to be an actual ID

# --- write job end time
echo "Job ended at $(date)"