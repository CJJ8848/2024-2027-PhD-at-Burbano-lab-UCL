#!/bin/bash

# Check if input file is provided
if [[ -z "$1" ]]; then
    echo "Usage: $0 <input_fasta_file>"
    exit 1
fi

# Input FASTA file
fasta_file="$1"

# Derive the output file name
output_file="${fasta_file%.*}_lengths_sorted.txt"

# Temporary unsorted file
unsorted_file="${fasta_file%.*}_lengths_unsorted.txt"

# Initialize variables
sample=""
length=0

# Create the unsorted file and add a header
echo -e "Sample\tLength" > "$unsorted_file"

# Read the FASTA file line by line
while IFS= read -r line; do
    if [[ $line == \>* ]]; then
        # If it's a header line, save the previous sample and length
        if [[ -n $sample ]]; then
            echo -e "$sample\t$length" >> "$unsorted_file"
        fi
        # Extract the sample name (without the ">")
        sample="${line#>}"
        length=0
    else
        # Add the length of the current line to the sequence length
        length=$((length + ${#line}))
    fi
done < "$fasta_file"

# Save the last sample and length
if [[ -n $sample ]]; then
    echo -e "$sample\t$length" >> "$unsorted_file"
fi

# Sort the file by length (2nd column) in descending order
sort -n -k2 -t$'\t'  "$unsorted_file" > "$output_file"

# Remove the unsorted file
rm "$unsorted_file"

echo "Sample names and lengths saved to $output_file"
