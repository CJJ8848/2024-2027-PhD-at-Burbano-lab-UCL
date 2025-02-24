#!/bin/bash

# Define input and output files
input_fasta="HTF_GC00000079_r1.fasta"
output_cleaned_fasta="HTF_GC00000079_r1noNs.fasta"
output_lengths="HTF_GC00000079_r1noNs.txt"
#rm all newlines in sequences and clean up the FASTA file
awk 'BEGIN {RS=">"; ORS=""} NR > 1 {
    header = $1;
    seq = $0;
    sub(header "\\n", "", seq);  # Remove the header line from the sequence block
    gsub("\\n", "", seq);  # Remove all newlines from the sequence
    gsub("-", "", seq);  # Remove dashes
    print ">" header "\n" seq "\n" >> "'"$output_cleaned_fasta"'";
}' "$input_fasta"

awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next} {seqlen += length($0)} END {if (seqlen) print seqlen}' $output_cleaned_fasta | sort | uniq -c > $output_lengths

# Notify user of completion
echo "cleaned FASTA written to $output_cleaned_fasta, and sequence lengths with counts written to $output_lengths." 
