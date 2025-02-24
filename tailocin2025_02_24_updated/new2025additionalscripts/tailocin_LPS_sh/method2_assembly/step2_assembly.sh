#!/bin/bash
#$ -l tmem=8G
#$ -l h_vmem=8G
#$ -l h_rt=3:30:0
#$ -S /bin/bash
#$ -N s2_assemble_matrix
#$ -t 1-131
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/

source /home/jiajucui/miniconda3/bin/activate phylogeny_snp

# Define directories
SAMPLE_LIST="/SAN/ugi/plant_genom/jiajucui/phylogeny/samplelist/tailocin46_85_m2.txt"
GENE_COORDS="/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/after_step2_coordinates.bed"

SUBSET_FASTQ_DIR="/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/results/method2/subset_fastq"
ASSEMBLY_DIR="/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/results/method2/assemblies"
MINIMAP_OUTPUT_DIR="/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/results/method2/minimap_results"
REFERENCE_GENOME="/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/results/output.fasta"
TOOLS="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/tailocin_2024_TF_Tapemeasure/shfortailocin/tools"

# Output matrix file
COVERAGE_MATRIX_FILE="/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/results/method2/gene_coverage_matrix.tsv"

# Create necessary directories
mkdir -p "$ASSEMBLY_DIR" "$MINIMAP_OUTPUT_DIR"

# Load SPAdes and Minimap2 paths
SPADES_CMD="spades.py"
MINIMAP_CMD="$TOOLS/minimap2/minimap2"

# Extract unique sample and gene names
samples=($(cat "$SAMPLE_LIST"))
genes=($(awk '{print $4}' "$GENE_COORDS" | sort -u))

# Initialize the matrix file with header (only on the first job)
if [[ "$SGE_TASK_ID" -eq 1 ]]; then
    echo -e "Gene\t${samples[*]}" > "$COVERAGE_MATRIX_FILE"
fi

# Get the current sample based on task ID
i=$SGE_TASK_ID
sample=$(sed -n "${i}p" "$SAMPLE_LIST")

# Temporary file to store results for this sample
TMP_FILE="/tmp/${sample}_coverage.tsv"
echo -e "Sample\tGene\tCovered_Proportion (%)" > "$TMP_FILE"

# Loop through each gene in the BED file
while read -r chrom start end gene strand; do
    # Define subset FASTQ filename
    subset_fastq="${SUBSET_FASTQ_DIR}/${sample}_${gene}_subset.fastq.gz"

    # Define output directories and files
    sample_assembly_dir="${ASSEMBLY_DIR}/${sample}_${gene}"
    contig_file="${sample_assembly_dir}/contigs.fasta"
    paf_file="${MINIMAP_OUTPUT_DIR}/${sample}_${gene}.paf"

    # Check if subset FASTQ file exists
    if [[ ! -f "$subset_fastq" ]]; then
        echo "WARNING: No subset FASTQ found for $sample - $gene"
        continue
    fi

    # Run SPAdes assembly (uncomment to enable)
    echo "Running SPAdes for Sample: $sample, Gene: $gene"
    mkdir -p "$sample_assembly_dir"
    $SPADES_CMD -s "$subset_fastq" -o "$sample_assembly_dir" --memory 8 -k 9,11,15

done < "$GENE_COORDS"


echo "Processing complete for $sample. Matrix updated in $COVERAGE_MATRIX_FILE."
