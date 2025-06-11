#!/bin/bash
# Combine Prokka txt files into one table

# Check if at least one argument is provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 file1.txt file2.txt ... output.txt"
    exit 1
fi

# Get the output file name from the last argument
output_file="${@: -1}"

# Get the list of input files
input_files=("${@:1:$#-1}")

# Initialize the output file and add the header row
echo -e "Filename\torganism\tcontigs\tbases\tCDS\trRNA\trepeat_region\ttRNA\ttmRNA" > "$output_file"

# Loop through all input files
for file in "${input_files[@]}"; do
    # Extract values from the file
    organism=$(grep "organism:" "$file" | cut -d ' ' -f2-)
    contigs=$(grep "contigs:" "$file" | cut -d ' ' -f2)
    bases=$(grep "bases:" "$file" | cut -d ' ' -f2)
    CDS=$(grep "CDS:" "$file" | cut -d ' ' -f2)
    rRNA=$(grep "rRNA:" "$file" | cut -d ' ' -f2)
    repeat_region=$(grep "repeat_region:" "$file" | cut -d ' ' -f2)
    tRNA=$(grep "tRNA:" "$file" | cut -d ' ' -f2)
    tmRNA=$(grep "tmRNA:" "$file" | cut -d ' ' -f2)

    # Append the values as a new row in the output file
    echo -e "$file\t$organism\t$contigs\t$bases\t$CDS\t$rRNA\t$repeat_region\t$tRNA\t$tmRNA" >> "$output_file"
done

echo "Combined summary with filenames as rows saved to $output_file"