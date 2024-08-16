#!/bin/bash

# Define the path to the directory containing the .final.fasta files
path="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_modern85/haplotype_selected"
#path="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/haplotype_selected"
output_file="summarym.txt"

# Initialize the output file
echo -e "sample\tTFA\tHTF\tgroup" > $output_file

# Iterate over each .final.fasta file in the directory
for file in "$path"/*.final.fasta; do
    # Extract the sample name from the file name
    sample=$(basename "$file" .final.fasta)

    # Initialize TFA and HTF as NULL
    TFA="NULL"
    HTF="NULL"

    # Check for the presence of TFA and HTF in the file
    if grep -q ">TFA" "$file"; then
        TFA=$(grep -o ">TFA[^ ]*" "$file" | tr -d '>')
    fi
    if grep -q ">HTF" "$file"; then
        HTF=$(grep -o ">HTF[^ ]*" "$file" | tr -d '>')
    fi

    # Determine the group
    if [ "$TFA" == "NULL" ] && [ "$HTF" == "NULL" ]; then
        group="noboth"
    elif [ "$TFA" == "NULL" ]; then
        group="noTFA"
    elif [ "$HTF" == "NULL" ]; then
        group="noHTF"
    else
        group=""
    fi

    # Append the results to the output file
    echo -e "$sample\t$TFA\t$HTF\t$group" >> $output_file
done
