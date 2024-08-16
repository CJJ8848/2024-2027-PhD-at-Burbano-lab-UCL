#!/bin/bash

# Define the input FASTA file and the output file for the extracted region
output='/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/tapemeasure'
mkdir -p $output
step3=$output/step3
step2=$output/step2
step1=$output/step1
mkdir -p $step3
reftape=$step3/alltapemeasureforfish.fasta
step1fasta="$step1/noNs_extracted_region.fasta"
step2fasta="$step2/modernNs30_tape_MSA.fasta"

cat $step1fasta $step2fasta | grep -v "^>[^:]*:[^-]*-.*$" > $reftape

