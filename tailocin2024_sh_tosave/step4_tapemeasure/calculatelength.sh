#!/bin/bash

# Replace 'filename.fasta' with the path to your FASTA file
FASTA_FILE="../../tailocin_extract/tapemeasure/step3/alltapemeasureforfish.fasta"

awk '/^>/ {if (seq) print name, length(seq); name=$0; seq=""} !/^>/ {seq=seq""$0} END {print name, length(seq)}' "$FASTA_FILE" > ../../tailocin_extract/tapemeasure/step3/alltapemeasureforfishlength.txt
