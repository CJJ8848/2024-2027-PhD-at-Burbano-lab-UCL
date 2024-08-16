#!/bin/bash

# Input FASTA file
output_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_modern85/"
rm -r $output_dir/haplotype_selected/merged_hmfa/

mkdir -p $output_dir/haplotype_selected/merged_hmfa/
input_fasta1="$output_dir/haplotype_selected/all_modern76fa_markexceptlongest.fasta"
input_fasta2=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/haplotype_selected/all_historicalfa_markexceptlongest.fasta
input_fasta="$output_dir/haplotype_selected/merged_hmfa/all_merged_hmfa_samples_tailocin_markexceptlongest.fasta"
cat /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_modern85/haplotype_selected/*markexceptlongest.fasta > $input_fasta1
cat /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/haplotype_selected/*markexceptlongest.fasta > $input_fasta2
cat $input_fasta1 $input_fasta2 > $input_fasta
cp $input_fasta /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/haplotype_selected/
summary_file="$output_dir/haplotype_selected/merged_hmfa/all_87fa_samples_tailocin_nonN_markexceptlongest.txt"
# Clear summary file if it exists
> "$summary_file"


# Parse the FASTA file and extract coordinates of non-N regions
awk '/^>/{if (seq) {print header"\t"seq}; header=$0; seq=""} !/^>/{seq=seq""$0} END{print header"\t"seq}' "$input_fasta" | while IFS=$'\t' read -r header sequence; do
  sample=$(echo "$header" | sed 's/^>//')
  position=0
  echo "$sequence" | grep -o . | while read -r base; do
    position=$((position + 1))
    if [[ "$base" != "N" ]]; then
      echo -e "${sample}\t${position}\t1" >> "$summary_file"
    fi
  done
done

echo "Summary file created at $summary_file"
