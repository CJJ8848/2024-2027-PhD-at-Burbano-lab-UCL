#!/bin/bash

# Define the base directories
bam_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46"
fastq_dir="/SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/2024_413_withcatmerge143"
link_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46_fastq"

# Create the directory for symbolic links
mkdir -p "$link_dir"

# Loop through each BAM file to get the sample names
for bam in "$bam_dir"/*.bam; do
  # Extract sample name from BAM file name
  samplename=$(basename "$bam" .mapped_to_Pseudomonas.dd.q20.bam)
  
  # Determine the correct raw FASTQ file name pattern
  if [[ "$samplename" == *.* ]]; then
    pattern="${samplename//./_}.fastq.gz"
  else
    pattern="${samplename}.fastq.gz"
  fi

  # Find the FASTQ file using the determined pattern
  fastq_file=$(ls "$fastq_dir" | grep  "$pattern")

  # Check if the FASTQ file exists and create the symbolic link
  if [ -n "$fastq_file" ]; then
    ln -s "$fastq_dir/$fastq_file" "$link_dir/"
  else
    echo "FASTQ file for sample $samplename not found."
  fi
done
fastq_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46_fastq"

# Change to the directory
cd "$fastq_dir" || exit

# Loop through each file in the directory
for file in *_*.fastq.gz; do
  # Create the new filename by replacing underscores with dots
  new_file=$(echo "$file" | sed 's/_/./1')
  
  # Rename the file
  mv "$file" "$new_file"
done

echo "Symbolic links created in $link_dir."
