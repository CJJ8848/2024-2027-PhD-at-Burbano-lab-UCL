#!/bin/bash
#$ -l tmem=100G
#$ -l h_vmem=100G
#$ -l h_rt=24:30:0
#$ -S /bin/bash
#$ -N vcftofastaandclustalotree
#$ -o /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/logs/
mkdir -p /SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/logs



# Activate the conda environment
source /home/jiajucui/miniconda3/bin/activate phylogeny_snp

# Define the base directories
bam_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46"
fastq_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_46_fastq"
output_dir="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract"
mkdir -p $output_dir
tmp_dir="${output_dir}/tmp"
readids_dir="${output_dir}/readids"

trimmed_fastq_dir="${output_dir}/trimmed_tailocin_fastq/"
assembly_dir="${output_dir}/assemblies"
mapping_dir="${output_dir}/mappings"
msa_dir="${output_dir}/msa"
new_tree_dir="${output_dir}/tree"
nocontigs_file="${msa_dir}/nocontigs.txt"
contigmapping_file="${mapping_dir}/contigmapping.txt"
haplotype_dir="${output_dir}/haplotype_selected"
nonN_file="${haplotype_dir}/nonN_TFAandHTF.txt"

# rm -r msa/ mappings/ assemblies/ contig_stats/ tree/
rm -r $assembly_dir $mapping_dir $msa_dir $new_tree_dir $haplotype_dir
#rm -r $mapping_dir $msa_dir $new_tree_dir $vcf_dir

mkdir -p $assembly_dir  $mapping_dir $msa_dir $new_tree_dir $haplotype_dir

tailocin_fasta="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/tailocin_region.fa"
reference_genome="/SAN/ugi/plant_genom/jiajucui/1_initial_data/reference_genome_Ps_with_tailocin_haplotypes/Pseudomonas.plate25.C2.pilon.contigs_renamed.with_Tail_Fiber_Haps.fasta"
tailocin_region="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_extract/regions.txt"

#22:15502-33559 but all the 40 genes
#TFA_p23.B8:1-513
#TFA_p26.D6:1-513
#TFA_p21.F9:1-498
#TFA_p25.C2:1-513
#TFA_p5.D5:1-513
#TFA_p25.A12:1-546
#TFA_p7.G11:1-546
#HTF_p23.B8:1-1803
#HTF_p26.D6:1-1803
#HTF_p21.F9:1-1245
#HTF_p25.C2:1-1803
#HTF_p5.D5:1-1803
#HTF_p25.A12:1-1383
#HTF_p7.G11:1-1830
#in the formated_region, all the TFA and HTF are reversed, make sure all of them are strand + like tailocin first
#I double checked with omega and found the *p25.C2 are the reference indeed but in reverse strand, so do the reverse and dont need to include the extracted ref chunk for TFA and HTF. only 7 and 7
# Read the tailocin region and process it
>"$tailocin_fasta"
while read -r region name; do
  if [ "$name" == "tailocin" ]; then
    samtools faidx "$reference_genome" "$region" | sed "1s/.*/>$name/" >> "$tailocin_fasta"
  else
    samtools faidx "$reference_genome" "$region" | seqtk seq -r | sed "1s/.*/>$name/" >> "$tailocin_fasta"
  fi
done < "$tailocin_region"
# Initialize or clear the nocontigs_file and contigmapping_file
> "$nocontigs_file"
> "$contigmapping_file" 
> "$nonN_file"

# Tools path
tools=/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/shfortailocin/tools




# Define regions and sequences
declare -A segment_lengths
segment_lengths["tailocin"]=18057
#segment_lengths["TFA_22_fragment_1_21_sp|P03740|TFA_LAMBD_Tail_fiber_assembly_protein"]=$((18230 - 17718 + 1))
segment_lengths["TFA_p23.B8"]=513
segment_lengths["TFA_p26.D6"]=513
segment_lengths["TFA_p21.F9"]=498
segment_lengths["TFA_p25.C2"]=513
segment_lengths["TFA_p5.D5"]=513
segment_lengths["TFA_p25.A12"]=546
segment_lengths["TFA_p7.G11"]=546
#segment_lengths["HTF_22_fragment_1_22_hypothetical_protein"]=$((20043 - 18241 + 1))
segment_lengths["HTF_p23.B8"]=1803
segment_lengths["HTF_p26.D6"]=1803
segment_lengths["HTF_p21.F9"]=1245
segment_lengths["HTF_p25.C2"]=1803
segment_lengths["HTF_p5.D5"]=1803
segment_lengths["HTF_p25.A12"]=1383
segment_lengths["HTF_p7.G11"]=1830
segment_names=(
    "tailocin"
    "TFA_p23.B8"
    "TFA_p26.D6"
    "TFA_p21.F9"
    "TFA_p25.C2"
    "TFA_p5.D5"
    "TFA_p25.A12"
    "TFA_p7.G11"
    "HTF_p23.B8"
    "HTF_p26.D6"
    "HTF_p21.F9"
    "HTF_p25.C2"
    "HTF_p5.D5"
    "HTF_p25.A12"
    "HTF_p7.G11"
)



#SPAdes
for bam in "$bam_dir"/*.bam; do
  samplename=$(basename "$bam" .mapped_to_Pseudomonas.dd.q20.bam)
  subset_fastq="${trimmed_fastq_dir}/${samplename}_subset.fastq.gz"
  # for 143 samples_cat_merged use --12 
  # --12 <file_name> File with interlaced forward and reverse paired-end reads.
  # --merged <file_name> File with merged paired reads. If the properties of the library permit, overlapping paired-end reads can be merged using special software.
  # Non-empty files with (remaining) unmerged left/right reads (separate or interlaced) must be provided for the same library for SPAdes to correctly detect the original read length.
  # but no additional files, --12 is ok
  # --isolate - isolate (standard) bacterial data;
  # --meta The --meta mode is designed for metagenomic data, which typically involves a heterogeneous mixture of reads. also good for dealing with a mixture of read qualities
  # --only-error-correction Performs read error correction only.
  # --only-assembler Runs assembly module only. if you have high quality reads, but here we have a mix of qualities, better run error correction
  # checked --meta gave more info in low quality reads like in120, it was nothing but here 30 contigs 
  # Determine the correct raw FASTQ file name pattern
  #but the PL0042 process is frozen even using -t 2 and require 50Gb, so remove --meta for it, and it works
  if [[ "$samplename" == *.* ]]; then
    subset_fastq1="${trimmed_fastq_dir}/${samplename}_subsetR1.fastq.gz"
    subset_fastq2="${trimmed_fastq_dir}/${samplename}_subsetR2.fastq.gz"
    
    spades.py --merge "$subset_fastq" -1 $subset_fastq1 -2 $subset_fastq2  --careful -o "$assembly_dir/$samplename" -k 21,33
  else
  # single end has not meta just use multicell/isolate as default
    spades.py -s "$subset_fastq"  -o "$assembly_dir/$samplename"  --careful -k 21,33
  fi
  # Step 2: Check the number of contigs in each assembly
  contig_file="$assembly_dir/$samplename/contigs.fasta"
#  contig_count=$(grep -c "^>" "$contig_file")
#  echo "$samplename: $contig_count contigs" >> "$contig_stats_dir/contig_counts.txt"
  # Check if the contig file exists
  if [[ ! -f "$contig_file" ]]; then
    echo "$samplename" >> "$nocontigs_file"
    echo "Contig file $contig_file does not exist. Sample name $samplename added to $nocontigs_file."
    continue
  fi

  reference_genome="${output_dir}/tailocin_region.fa"


  # Define your output files
  paf_file="${mapping_dir}/${samplename}_mapped.paf"
  fasta_out="${msa_dir}/${samplename}_tailocin_region.fasta"
  #  Run minimap2 to get the mapping
  #sergio's idea: map the reference short chunk of hyplotypes to the long contig, in principle works for both historical and modern samples.
  #$tools/minimap2/minimap2 -cx asm5 "$reference_genome" "$contig_file" > "$paf_file"
  mkdir -p ${msa_dir}/fake
  fasta_outfake="${msa_dir}/fake/${samplename}_fake.fasta"
  fasta_out1="${msa_dir}/${samplename}_tailocin_region_allconcatenated.fasta"
  $tools/minimap2/minimap2 -cx asm5 "$contig_file" "$reference_genome" > "$paf_file"


  if [[ $(less "$paf_file" | wc -l) -eq 0 ]]; then
    echo "$samplename nomappedcontig" >> "$nocontigs_file"
    echo "Contig file $contig_file does not exist. Sample name $samplename added to $nocontigs_file."
    continue
  fi
  # Calculate cumulative start positions
  declare -A cumulative_starts
  cumulative_starts["tailocin"]=0
  for i in "${!segment_names[@]}"; do
    if [[ $i -gt 0 ]]; then
      prev_segment="${segment_names[$((i-1))]}"
      cumulative_starts["${segment_names[$i]}"]=$((cumulative_starts["$prev_segment"] + segment_lengths["$prev_segment"]))
    fi
  done

  total_length=0
  for length in "${segment_lengths[@]}"; do
    total_length=$((total_length + length))
  done
  final_sequence=$(printf 'N%.0s' $(seq 1 $total_length))
  fakefinal_sequence=$(printf 'N%.0s' $(seq 1 $total_length))


  replace_sequence() {
    local start=$1
    local seq=$2
    final_sequence="${final_sequence:0:start}${seq}${final_sequence:$(($start + ${#seq}))}"
  }

  tailocin_count=0
  htf_count=0
  taf_count=0
  declare -a htf_list
  declare -a tfa_list
  declare -A nonN_counts

  while read -r line; do
    ref_name=$(echo "$line" | awk '{print $1}')
    ref_start=$(echo "$line" | awk '{print $3}')
    ref_end=$(echo "$line" | awk '{print $4}')
    strand=$(echo "$line" | awk '{print $5}')
    contig_name=$(echo "$line" | awk '{print $6}')
    contig_start=$(echo "$line" | awk '{print $8}')
    contig_end=$(echo "$line" | awk '{print $9}')

    if [[ "$strand" == "-" ]]; then
      contig_seq=$(samtools faidx "$contig_file" "$contig_name:$contig_start-$contig_end" | seqtk seq -r - | tail -n +2 | tr -d '\n' )
    else
      contig_seq=$(samtools faidx "$contig_file" "$contig_name:$contig_start-$contig_end" | tail -n +2 | tr -d '\n')
    fi

    segment_start=$((ref_start + cumulative_starts["$ref_name"]))
    replace_sequence $segment_start "$contig_seq"

    if [[ "$ref_name" == "tailocin" ]]; then
        tailocin_count=$((tailocin_count + 1))
    elif [[ "$ref_name" == HTF* ]]; then
        htf_count=$((htf_count + 1))
        htfname=${ref_name%%:*}
        htf_list+=("$htfname")
    elif [[ "$ref_name" == TFA* ]]; then
        tfa_count=$((tfa_count + 1))
        tfaname=${ref_name%%:*}
        tfa_list+=("$tfaname")
    fi
  done < "$paf_file"

  # Write the final sequence to fasta_out with individual segment headers
  echo ">${samplename}" > "$fasta_out1"
  echo "$final_sequence" >> "$fasta_out1"
  
  ## Write the final sequence to fasta_out with individual segment headers
  {
    for segment in "${segment_names[@]}"; do
      start=${cumulative_starts[$segment]}
      length=${segment_lengths[$segment]}
      segment_sequence=${final_sequence:$start:$length}
      echo ">$segment"
      echo "$segment_sequence"
    done
  } > "$fasta_out"

 # if [[ $tailocin_count -eq 0 && $htf_count -eq 0 && $tfa_count -eq 0 ]]; then
 #   rm $fasta_out $fasta_out1
 #   echo "$samplename no mapped contig"
 # fi

  # Append to the combined MSA file before formatting individual segments
  cat "$fasta_out1" >> "$msa_dir/all_historicalfa_samples_tailocin.fasta"

  # Output the mapping summary and lists to the contigmapping file
  {
    echo "$samplename: tailocin mapped contigs = $tailocin_count, TFA mapped contigs = $tfa_count, HTF mapped contigs = $htf_count"
    
    echo "HTF List:"
    for htf_segment in "${htf_list[@]}"; do
        echo "$htf_segment"
    done

    echo "TFA List:"
    for tfa_segment in "${tfa_list[@]}"; do
        echo "$tfa_segment"
    done
  } >> "$contigmapping_file"

  # Process segments to calculate non-N counts and proportions
  # Process segments to calculate non-N counts and proportions
  fakefornonN_replace_sequence() {
    local start=$1
    local seq=$2
    fakefinal_sequence="${fakefinal_sequence:0:start}${seq}${fakefinal_sequence:$(($start + ${#seq}))}"
   }

  while read -r line; do
    ref_name=$(echo "$line" | awk '{print $1}')
    ref_start=$(echo "$line" | awk '{print $3}')
    ref_end=$(echo "$line" | awk '{print $4}')
    strand=$(echo "$line" | awk '{print $5}')
    contig_name=$(echo "$line" | awk '{print $6}')
    contig_start=$(echo "$line" | awk '{print $8}')
    mappedbase=$(echo "$line" | awk '{print $10}')
    contig_end=$((contig_start + mappedbase))
    if [[ "$strand" == "-" ]]; then
      contig_seq=$(samtools faidx "$contig_file" "$contig_name:$contig_start-$contig_end" | seqtk seq -r - | tail -n +2 | tr -d '\n' )
    else
      contig_seq=$(samtools faidx "$contig_file" "$contig_name:$contig_start-$contig_end" | tail -n +2 | tr -d '\n')
    fi

    segment_start=$((ref_start + cumulative_starts["$ref_name"]))
    fakefornonN_replace_sequence $segment_start "$contig_seq"

  done < "$paf_file"
 ## Write the final sequence to fasta_out with individual segment headers
  {
    for segment in "${segment_names[@]}"; do
      start=${cumulative_starts[$segment]}
      length=${segment_lengths[$segment]}
      segment_sequence=${fakefinal_sequence:$start:$length}
      echo ">$segment"
      echo "$segment_sequence"
    done
  } > "$fasta_outfake"
  
  { 
    for segment in "${tfa_list[@]}"; do
      sequence=$(grep -A1 ">${segment}" "$fasta_outfake" | tail -n1)
      nonN_count=$(echo "$sequence" | tr -cd 'ATCGatcg' | wc -c)
      proportion=$(echo "scale=5; $nonN_count/${segment_lengths[$segment]}" | bc)
      nonN_counts["$segment"]=$proportion
      echo "$samplename >${segment}:$proportion"
    done
    for segment in "${htf_list[@]}"; do
      sequence=$(grep -A1 ">${segment}" "$fasta_outfake" | tail -n1)
      nonN_count=$(echo "$sequence" | tr -cd 'ATCGatcg' | wc -c)
      proportion=$(echo "scale=5; $nonN_count/${segment_lengths[$segment]}" | bc)
      nonN_counts["$segment"]=$proportion
      echo "$samplename >${segment}:$proportion"
    done
  } >> "$nonN_file"


  longest_tfa=$(printf "%s\n" "${!nonN_counts[@]}" | grep "^TFA" | while read segment; do echo "$segment ${nonN_counts[$segment]}"; done | sort -k2,2nr | head -n1 | awk '{print $1}')
  longest_htf=$(printf "%s\n" "${!nonN_counts[@]}" | grep "^HTF" | while read segment; do echo "$segment ${nonN_counts[$segment]}"; done | sort -k2,2nr | head -n1 | awk '{print $1}')
  # Filter out segments that are still all Ns
  if [[ -n "$longest_tfa" ]]; then
    longest_tfa_seq=$(grep -A1 ">${longest_tfa%%:*}" "$fasta_out" | tail -n1)
    longest_tfa_nonN=$(echo "$longest_tfa_seq" | tr -cd 'ATCGatcg' | wc -c)
    if [[ "$longest_tfa_nonN" -eq 0 ]]; then
      longest_tfa=""
    fi
  fi

  if [[ -n "$longest_htf" ]]; then
    longest_htf_seq=$(grep -A1 ">${longest_htf%%:*}" "$fasta_out" | tail -n1)
    longest_htf_nonN=$(echo "$longest_htf_seq" | tr -cd 'ATCGatcg' | wc -c)
    if [[ "$longest_htf_nonN" -eq 0 ]]; then
      longest_htf=""
    fi
  fi

  # Create the final FASTA file with the selected segments
  final_fasta="${haplotype_dir}/${samplename}.final.fasta"
  {
    echo ">tailocin"
    grep -A1 ">tailocin" "$fasta_out" | tail -n1

    if [[ -n "$longest_tfa" ]]; then
      echo ">${longest_tfa%%:*}"
      grep -A1 ">${longest_tfa%%:*}" "$fasta_out" | tail -n1
    fi

    if [[ -n "$longest_htf" ]]; then
      echo ">${longest_htf%%:*}"
      grep -A1 ">${longest_htf%%:*}" "$fasta_out" | tail -n1
    fi
  } > "$final_fasta"

 #only longest all other Ns
 # Create the final concatenated FASTA file
  final2_fasta="${haplotype_dir}/${samplename}.markexceptlongest.fasta"
  {
  echo ">${samplename}"
  final_sequence=$(printf 'N%.0s' $(seq 1 $total_length))
  # Retain the tailocin sequence
  tailocin_sequence=$(grep -A1 ">tailocin" "$fasta_out" | tail -n1)
  tailocin_start=${cumulative_starts["tailocin"]}
  replace_sequence $tailocin_start "$tailocin_sequence"
  # Retain the longest TFA sequence and mark others as Ns
  for segment in "${tfa_list[@]}"; do
    segment_start=${cumulative_starts["$segment"]}
    if [[ "$segment" == "${longest_tfa%%:*}" ]]; then
      tfa_sequence=$(grep -A1 ">${segment}" "$fasta_out" | tail -n1)
      replace_sequence $segment_start "$tfa_sequence"
    fi
  done

  # Retain the longest HTF sequence and mark others as Ns
  for segment in "${htf_list[@]}"; do
    segment_start=${cumulative_starts["$segment"]}
    if [[ "$segment" == "${longest_htf%%:*}" ]]; then
      htf_sequence=$(grep -A1 ">${segment}" "$fasta_out" | tail -n1)
      replace_sequence $segment_start "$htf_sequence"
    fi
  done

  echo "$final_sequence"
  } > "$final2_fasta"

 # if [[ $tailocin_count -eq 0 && $htf_count -eq 0 && $tfa_count -eq 0 ]]; then
 #   rm $final_fasta 
 #   echo "$samplename no mapped contig"
 # fi

  echo "The final sequence has been saved to $final_fasta"
done


# Create output files for HTF and TFA
htf_fasta="${haplotype_dir}/all_HTF_samples.fasta"
tfa_fasta="${haplotype_dir}/all_TFA_samples.fasta"
> "$htf_fasta"
> "$tfa_fasta"

# Iterate over each final fasta file and extract HTF and TFA sequences
for final_fa in "${haplotype_dir}"/*.final.fasta; do
  samplename=$(basename "$final_fa" .final.fasta)
  
  # Extract HTF sequence
  grep -A1 ">HTF" "$final_fa" | sed "s/^>/>${samplename}|/" >> "$htf_fasta"
  
  # Extract TFA sequence
  grep -A1 ">TFA" "$final_fa" | sed "s/^>/>${samplename}|/" >> "$tfa_fasta"
done

echo "HTF and TFA multi-sample FASTA files have been generated in $haplotype_dir."


#filter with 50% proportion covered
#!/bin/bash

# Define input and output directories

# Create filtered output files for HTF and TFA
filtered_htf_fasta="${haplotype_dir}/filtered_HTF_samples.fasta"
filtered_tfa_fasta="${haplotype_dir}/filtered_TFA_samples.fasta"
> "$filtered_htf_fasta"
> "$filtered_tfa_fasta"

# Threshold for non-N proportion
threshold=0.65

# Filter HTF sequences
while read -r line; do
  if [[ $line == ">"* ]]; then
    header=$line
    sequence=""
  else
    sequence=$line
    nonN_count=$(echo "$sequence" | tr -cd 'ATCGatcg' | wc -c)
    total_length=${#sequence}
    proportion=$(echo "scale=5; $nonN_count / $total_length" | bc)
    if (( $(echo "$proportion >= $threshold" | bc -l) )); then
      echo "$header" >> "$filtered_htf_fasta"
      echo "$sequence" >> "$filtered_htf_fasta"
    fi
  fi
done < "${haplotype_dir}/all_HTF_samples.fasta"

# Filter TFA sequences
while read -r line; do
  if [[ $line == ">"* ]]; then
    header=$line
    sequence=""
  else
    sequence=$line
    nonN_count=$(echo "$sequence" | tr -cd 'ATCGatcg' | wc -c)
    total_length=${#sequence}
    proportion=$(echo "scale=5; $nonN_count / $total_length" | bc)
    if (( $(echo "$proportion >= $threshold" | bc -l) )); then
      echo "$header" >> "$filtered_tfa_fasta"
      echo "$sequence" >> "$filtered_tfa_fasta"
    fi
  fi
done < "${haplotype_dir}/all_TFA_samples.fasta"

echo "Filtered HTF and TFA multi-sample FASTA files have been generated in $haplotype_dir."

#try use directly the msa by cat merge
# Step 5: Build a phylogenetic tree
#iqtree -s $htf_fasta  -nt AUTO -bb 1000 -alrt 1000 -pre "$new_tree_dir/h_HTF"


