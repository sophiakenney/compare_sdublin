#!/bin/bash
# Read Classification - Kraken2

# --- write job start time to log file
echo "Job started at $(date)"

# --- Set the directory trimmed files
WD="/your/path/here"

cd ${WD}/1.0-trimmed

mkdir -p ${WD}/2.0-kraken2

# --- activate bioinfo env for numpy
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate bioinfo

# --- Set variable paths
KRAKEN_DB="/your/path/to/kraken2/krakendbstandard"
OUTPUT="/your/path/here/2.0-kraken2" # set to directory where you'd like your output
KRAKEN2_DIR="/your/path/to/kraken2"

# --- Loop over all .fq files
for file in "$WD"/*_1P.fq.gz; do

    fname=$(basename $file)
    filename="${fname%_*}"

    echo "${filename} running"
    echo " "

    echo "kraken2 2"
    ${KRAKEN2_DIR}/kraken2 --threads 16 --db ${KRAKEN_DB} \
    --gzip-compressed \
    --paired --classified-out ${OUTPUT}/${filename}.classified-out.R#.fastq \
    --unclassified-out ${OUTPUT}/${filename}.unclassified-out.R#.fastq \
    --report ${OUTPUT}/"${filename}.k2report" \
    --report-minimizer-data \
    --minimum-hit-groups 3 \
    ${filename}_1P.fq.gz ${filename}_2P.fq.gz > ${OUTPUT}/"${filename}.kraken2" &
    wait $!
    echo "kraken2 2 complete"
done

# --- write job end time to log file
echo "Job ended at $(date)"