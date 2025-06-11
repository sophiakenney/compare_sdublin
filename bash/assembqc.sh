#!/bin/bash
# Assembly QC

# ---- Run CheckM ----

# --- write job start time
echo "Job started at $(date)"

# --- activate conda env
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate checkm

# --- set paths
export CHECKM_DATA_PATH=/your/path/here/checkM
ASSEMB=/your/path/here/3.0-unicycler_entero/assemblies
OUT=/your/path/here/3.1-uniqc/checkm

cd ${ASSEMB}

echo "CheckM started at $(date)"

  checkm lineage_wf -f "${OUT}/checkm.tsv" -x fasta ${ASSEMB} ${OUT} --tab_table

echo "CheckM ended at $(date)"

conda deactivate

# ---- Run QUAST ----

echo "QUAST started at $(date)"


python /your/path/to/quast/quast.py \
-o /your/path/to/assembqc/3.1-uniqc/quast \
--min-contig 500 \
--threads 4 \
${ASSEMB}/*_assembly.fasta

echo "QUAST ended at $(date)"

# ---write job end time
echo "Job ended at $(date)"