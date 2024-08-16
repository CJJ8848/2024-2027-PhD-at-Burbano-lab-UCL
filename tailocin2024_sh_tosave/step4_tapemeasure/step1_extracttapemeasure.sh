#!/bin/bash

# Define the input FASTA file and the output file for the extracted region
input_fasta="../../tailocin_extract/haplotype_selected/all_merged_hmfa_samples_tailocin_markexceptlongest.fasta"
output='/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/tapemeasure'
mkdir -p $output
step1=$output/step1
mkdir -p $step1
output_fasta="$step1/extracted_region_all.fasta"

# Clear the output file if it already exists
> "$output_fasta"

# Loop through each sequence name in the .fai index file
while read -r line; do
    # Extract the sequence name (first column in .fai file)
    seq_name=$(echo "$line" | cut -f1)

    # Extract the desired region using samtools faidx and append to the output file
    samtools faidx "$input_fasta" "${seq_name}:10249-10838" >> "$output_fasta"
done < "${input_fasta}.fai"

echo "Extraction complete. Extracted regions are saved in $output_fasta."


#!/bin/bash

# Define the input FASTA file and the output files for the extracted regions and those containing 'N's
input_fasta="../../tailocin_extract/haplotype_selected/all_merged_hmfa_samples_tailocin_markexceptlongest.fasta"
output_fasta="$step1/noNs_extracted_region.fasta"
output_modern_ns_fasta="$step1/modernNs.fasta"
output_historical_ns_fasta="$step1/historicalNs.fasta"

# Clear the output files if they already exist
> "$output_fasta"
> "$output_modern_ns_fasta"
> "$output_historical_ns_fasta"

# Loop through each sequence name in the .fai index file
while read -r line; do
    # Extract the sequence name (first column in .fai file)
    seq_name=$(echo "$line" | cut -f1)

    # Extract the desired region using samtools faidx
    extracted_region=$(samtools faidx "$input_fasta" "${seq_name}:10249-10838")

    # Modify the header to "${seq_name}_tape"
    modified_region=$(echo "$extracted_region" | sed "1s/.*/>${seq_name}_tape/")

    # Check if the extracted region contains 'N'
    if echo "$modified_region" | grep -q 'N'; then
        # Append to the appropriate file based on the sequence name prefix
        if [[ $seq_name == p* ]]; then
            echo "$modified_region" >> "$output_modern_ns_fasta"
        else
            echo "$modified_region" >> "$output_historical_ns_fasta"
        fi
    else
        # Append to the output file for extracted regions
        echo "$modified_region" >> "$output_fasta"
    fi
done < "${input_fasta}.fai"

echo "Extraction complete. Extracted regions are saved in $output_fasta. Modern sequences containing 'N' are saved in $output_modern_ns_fasta. Historical sequences containing 'N' are saved in $output_historical_ns_fasta."
