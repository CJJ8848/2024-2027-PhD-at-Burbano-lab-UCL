#!/bin/bash

# Input and output files
input_file="./step1/noNs_extracted_region.fasta"
output_file="noNs_formatted_sequences.fasta"

# Initialize variables
{ rm -f "$output_file"; touch "$output_file"; }

# Read through the file and format sequences
while IFS= read -r line || [ -n "$line" ]; do
    # Check if the line starts with '>', indicating a header line
    if [[ "$line" =~ ^\> ]]; then
        # If a sequence is buffered, write it to the output file
        if [[ -n "$sequence" ]]; then
            echo "$sequence" >> "$output_file"
        fi
        # Write the header to the output and reset the sequence buffer
        echo "$line" >> "$output_file"
        sequence=""
    else
        # Append the line to the sequence, removing any newlines
        sequence+="$line"
    fi
done < "$input_file"

# Write the last buffered sequence to the file
if [[ -n "$sequence" ]]; then
    echo "$sequence" >> "$output_file"
fi

echo "Formatted sequences saved to $output_file."
