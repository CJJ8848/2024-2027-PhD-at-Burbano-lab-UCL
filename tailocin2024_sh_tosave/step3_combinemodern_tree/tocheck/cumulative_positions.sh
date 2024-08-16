#!/bin/bash

# Define regions and sequences
declare -A segment_lengths
segment_lengths["tailocin"]=18057
segment_lengths["TFA_p23.B8"]=513
segment_lengths["TFA_p26.D6"]=513
segment_lengths["TFA_p21.F9"]=498
segment_lengths["TFA_p25.C2"]=513
segment_lengths["TFA_p5.D5"]=513
segment_lengths["TFA_p25.A12"]=546
segment_lengths["TFA_p7.G11"]=546
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


# Define output directory
output_dir=./
# Define output directory

# Define output path for cumulative positions file
cumulative_positions_file="${output_dir}/cumulative_positions.tsv"

# Calculate cumulative start and end positions
declare -A cumulative_starts
declare -A cumulative_ends

cumulative_starts["tailocin"]=1
cumulative_ends["tailocin"]=$((segment_lengths["tailocin"] ))

for i in "${!segment_names[@]}"; do
  if [[ $i -gt 0 ]]; then
    prev_segment="${segment_names[$((i-1))]}"
    cumulative_starts["${segment_names[$i]}"]=$((cumulative_starts["$prev_segment"] + segment_lengths["$prev_segment"]))
    cumulative_ends["${segment_names[$i]}"]=$((cumulative_starts["${segment_names[$i]}"] + segment_lengths["${segment_names[$i]}"] - 1))
  fi
done

# Output the cumulative start and end positions to a TSV file
echo -e "Segment\tStart\tEnd" > "$cumulative_positions_file"

for segment in "${segment_names[@]}"; do
  start=${cumulative_starts[$segment]}
  end=${cumulative_ends[$segment]}
  echo -e "${segment}\t${start}\t${end}" >> "$cumulative_positions_file"
done

echo "Cumulative start and end positions have been written to $cumulative_positions_file."
