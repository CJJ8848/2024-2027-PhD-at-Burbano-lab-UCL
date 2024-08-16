#!/bin/bash

# Define the input FASTA file and the output file for the extracted region
output='/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/tapemeasure'
mkdir -p $output
step2=$output/step2
mkdir -p $step2

input_fasta="$output/step1/modernNs.fasta"

# Function to extract and process each sample
process_sample() {
    local sample=$1
    local paf_file=$2
    local raw_fasta=$3
    tape_real=-1
    tape_start=-1
    tape_end=0
    tape_chr=0
    strand=0
    while read -r line; do
        ref_name=$(echo "$line" | awk '{print $1}')
        ref_start=$(echo "$line" | awk '{print $3}')
        ref_end=$(echo "$line" | awk '{print $4}')
        strand=$(echo "$line" | awk '{print $5}')
        contig_name=$(echo "$line" | awk '{print $6}')
        contig_start=$(echo "$line" | awk '{print $8}')
        contig_end=$(echo "$line" | awk '{print $9}')
       if [[ $strand == "+" ]]; then
        if [[ $ref_name == "tailocin" ]]; then
            if (( (ref_end - 10249) <= 300 && (ref_end - 10249) >= -300 )); then
                #tape_start="$contig_end"
                tape_start=$((10249 - ref_end + contig_end))
                #tape_end=$((tape_real + 631)) wrong eg in p6.A10 the gap length is not the same in contig
                tape_chr="$contig_name"
                strand='+'
            fi
            if (( (ref_start - 10838) <= 300 && (ref_start - 10838) >= -300 )); then
                
                tape_end=$((10838 - ref_start + contig_start))
                tape_chr2="$contig_name"
            fi
        fi
       else
        if [[ $ref_name == "tailocin" ]]; then
            if (( (ref_end - 10249) <= 300 && (ref_end - 10249) >= -300 )); then
                #tape_end="$contig_start"
                tape_end=$((ref_end-10249 + contig_start))
               # tape_real=$((tape_end - 631))
                tape_chr="$contig_name"
                strand='-'
            fi
            if (( (ref_start - 10838) <= 300 && (ref_start - 10838) >= -300 )); then
                tape_start=$((ref_start-10838 + contig_end)) 
                #the start is from 10838 to 10249, check p13.C1, and then I reverse it to 10249 to 10838.
                #tape_real=$((ref_start - 10838 + contig_end))
                #tape_end=$((tape_real + 631))               
                tape_chr2="$contig_name"
            fi
        fi
       fi
    done < "$paf_file"



    if [[ $tape_chr == "$tape_chr2" ]]; then
        tape_coor="${tape_chr}:${tape_start}-${tape_end}"
        if [[ "$strand" == "+" ]]; then
          echo ">${sample}_tape" >> "$step2/modernNs30_tape_MSA.fasta"
          samtools faidx "$raw_fasta" "$tape_coor" |seqtk seq >> "$step2/modernNs30_tape_MSA.fasta"
        else
          echo '-' $sample
          #samtools faidx "$raw_fasta" "$tape_coor" > "$step2/${sample}_tape.fasta"
          echo ">${sample}_tape" >> "$step2/modernNs30_tape_MSA.fasta"
          samtools faidx "$raw_fasta" "$tape_coor" | seqtk seq -r - >> "$step2/modernNs30_tape_MSA.fasta"
        fi
       # Save individual sample tape to a separate file
    fi
}

# Clear the modernNs30_tape_MSA.fasta if it already exists
> "$step2/modernNs30_tape_MSA.fasta"

# Extract sample names from modernNs.fasta and process each sample
grep ">" "$input_fasta" | sed 's/>//' | while read -r sample_tape; do
    sample=$(echo "$sample_tape" | sed 's/_tape//')
    echo $sample
    # Define paths for the PAF file and raw FASTA (you need to adjust these paths as necessary)
    raw_fasta="/SAN/ugi/plant_genom/jiajucui/phylogeny/phylogeny_read2tree/read2treeinput/pankmerwithpan85modernraw/rawfasta30nonOTU5_55OTU5/${sample}.fasta.bgz"
    paf_file="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_modern85/mappings/${sample}_mapped.paf"
    

    # Process each sample
    process_sample "$sample" "$paf_file" "$raw_fasta"
done

echo "Processing complete. Extracted regions are saved in $step2."
