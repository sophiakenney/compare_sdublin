#!/bin/bash
# Extract Enterobacteriaceae reads

# --- write job start time to log file
echo "Job started at $(date)"

# --- activate bioinfo env for numpy
source ~/.bashrc
eval "$(conda shell.bash hook)"
conda activate bioinfo

# --- Set the directory containing trimmed reads
TR="/your/path/here/1.0-trimmed" 
cd ${TR}

mkdir -p "your/path/here/2.1-k2filteredentero"

# --- Set variable paths
KRAKEN_DB="/your/path/to/kraken2/krakendbstandard"
KRAKEN2_DIR="/your/path/to/kraken2"
KRAKEN_TOOLS="/your/path/to/KrakenTools"
K2="/your/path/here/2.0-kraken2" # kraken2 output directory
FILT="/your/path/here/2.1-k2filteredentero" # extracted reads output directory

#--- Loop over all .fastq files

for file in *1P.fq.gz; do

    fname="${file%_1P*}"

    echo "${fname} running"
    echo " "

   # Start Read Extraction
   python ${KRAKEN_TOOLS}/extract_kraken_reads.py -t 543 \ # using Enterobacteriaceae taxID
    --include-children \
    -k ${K2}/${fname}.kraken2 \
    -s1 ${TR}/${fname}_1P.fq.gz \
    -s2 ${TR}/${fname}_2P.fq.gz \
    -o ${FILT}/${fname}.filt_1P.fq \
    -o2 ${FILT}/${fname}.filt_2P.fq \
    --fastq-output \
    -r ${K2}/${fname}.k2report &
    wait $!

    echo "${fname} done"
done

# --- Run seqkit stats to account for read loss from raw to filtered 

cd ${FILT}

seqkit stats *.fq -T > k2filtstats.txt


# --- write job end time to log file
echo "Job ended at $(date)"