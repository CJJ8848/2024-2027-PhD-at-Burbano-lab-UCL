#!/bin/bash

# Load environment
source /home/jiajucui/miniconda3/bin/activate phylogeny_snp

# Define directories
SAMPLE_LIST="/SAN/ugi/plant_genom/jiajucui/phylogeny/samplelist/tailocin46_85_m2.txt"
GENE_COORDS="/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/after_step2_coordinates.bed"

SUBSET_FASTQ_DIR="/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/results/method2/subset_fastq"
ASSEMBLY_DIR="/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/results/method2/assemblies"
MINIMAP_OUTPUT_DIR="/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/results/method2/minimap_results"
REFERENCE_GENOME="/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/results/output.fasta"
TOOLS="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/shfortailocin/tools"
mkdir -p /SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/results/method2/tmp
# Output file
#COVERAGE_MATRIX_FILE="/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/results/method2/gene_coverage_matrix.tsv"

# Create necessary directories
mkdir -p "$ASSEMBLY_DIR" "$MINIMAP_OUTPUT_DIR"

# Load SPAdes and Minimap2 paths
SPADES_CMD="spades.py"
MINIMAP_CMD="$TOOLS/minimap2/minimap2"

# Extract unique gene names
genes=($(awk '{print $4}' "$GENE_COORDS" | sort -u))

# Initialize the matrix file with header
echo -e "Gene\t$(cat "$SAMPLE_LIST" | tr '\n' '\t')" > "$COVERAGE_MATRIX_FILE"

# Process each sample in a for loop
for sample in $(cat "$SAMPLE_LIST"); do
    echo "Processing Sample: $sample"

    # Temporary file to store results for this sample
    TMP_FILE="/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/results/method2/tmp/${sample}_coverage.tsv"
    echo -e "Gene\t$sample" > "$TMP_FILE"

    # Loop through each gene in the BED file
    while read -r chrom start end gene strand; do
        subset_fastq="${SUBSET_FASTQ_DIR}/${sample}_${gene}_subset.fastq.gz"
        sample_assembly_dir="${ASSEMBLY_DIR}/${sample}_${gene}"
        contig_file="${sample_assembly_dir}/contigs.fasta"
        paf_file="${MINIMAP_OUTPUT_DIR}/${sample}_${gene}.paf"

        # Check if subset FASTQ file exists
        if [[ ! -f "$subset_fastq" ]]; then
            echo "WARNING: No subset FASTQ found for $sample - $gene"
            echo -e "$gene\t0" >> "$TMP_FILE"
            continue
        fi

        # Run SPAdes assembly (optional, uncomment if needed)
        echo "Running SPAdes for Sample: $sample, Gene: $gene"
        mkdir -p "$sample_assembly_dir"
        #$SPADES_CMD -s "$subset_fastq" -o "$sample_assembly_dir" --memory 8 -k 9,11,15

        # Run Minimap2
        if [[ -f "$contig_file" ]]; then
            echo "Running Minimap2 for Sample: $sample, Gene: $gene"
#            $MINIMAP_CMD -cx asm5 "$REFERENCE_GENOME" "$contig_file" > "$paf_file"

            # Compute covered proportion
            if [[ -s "$paf_file" ]]; then
		#$11 is the total length of alignment on the reference, $10 is the aligned block length but if there are mismatch (mutations) it would not be included in the $10, so we use $11 and the total length of the ref is $7
                covered_proportion=$(awk '{covered += $11; total += $7} END {if (total > 0) print covered/total * 100; else print 0}' "$paf_file")
            else
                covered_proportion=0
            fi
            echo -e "$gene\t$covered_proportion" >> "$TMP_FILE"
        else
            echo "WARNING: No contigs found for $sample - $gene"
            echo -e "$gene\t0" >> "$TMP_FILE"
        fi
    done < "$GENE_COORDS"


    echo "Finished processing $sample."
done

# Define directories and files
TMP_DIR="/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/results/method2/tmp"
OUTPUT_MATRIX="/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/results/method2/method2_afterassembly_gene_coverage_matrix.tsv"

# Find all TMP files
tmp_files=($(ls "$TMP_DIR"/*_coverage.tsv 2>/dev/null))

# Check if TMP files exist
if [[ ${#tmp_files[@]} -eq 0 ]]; then
    echo "❌ No TMP files found in $TMP_DIR"
    exit 1
fi

echo "✅ Found ${#tmp_files[@]} TMP files. Merging..."

# Extract gene names from the first TMP file
genes=($(awk 'NR>1 {print $1}' "${tmp_files[0]}" | sort -u))

# Extract sample names from TMP file names
samples=($(for file in "${tmp_files[@]}"; do basename "$file" | sed 's/_coverage.tsv//'; done))

# Initialize output file with header
echo -e "Gene\t${samples[*]}" > "$OUTPUT_MATRIX"

# Loop through each gene and collect data from all samples
for gene in "${genes[@]}"; do
    row="$gene"
    for sample in "${samples[@]}"; do
        tmp_file="${TMP_DIR}/${sample}_coverage.tsv"
        value=$(awk -v gene="$gene" '$1 == gene {print $2}' "$tmp_file")
        row+="\t${value:-0}"  # Default to 0 if missing
    done
    echo -e "$row" >> "$OUTPUT_MATRIX"
done

echo "✅ Gene coverage matrix saved to: $OUTPUT_MATRIX"
