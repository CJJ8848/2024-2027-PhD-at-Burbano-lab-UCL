#!/bin/bash
#$ -l tmem=4G
#$ -l h_vmem=4G
#$ -l h_rt=3:30:0
#$ -S /bin/bash
#$ -t 1-46
#$ -N vcftofastaandclustalotree
#$ -o /SAN/ugi/plant_genom/jiajucui/logs/
#$ -e /SAN/ugi/plant_genom/jiajucui/logs/
echo "Task id is $SGE_TASK_ID"
# Sample list file
SAMPLE_LIST="/SAN/ugi/plant_genom/jiajucui/phylogeny/samplelist/tailocin46_method2.txt"

# Define directories
BAM_DIR="/SAN/ugi/plant_genom/jiajucui/4_mapping_to_pseudomonas/2024_413_but280withnsqs10_233_HB37_ps"
FASTQ_DIR=/SAN/ugi/plant_genom/jiajucui/2_trimmed_merged/2024_413_withcatmerge143/
OUTPUT_DIR="/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/results/method2"
SUBSET_FASTQ_DIR="${OUTPUT_DIR}/subset_fastq"
GENE_COORDS="/SAN/ugi/plant_genom/jiajucui/tailocin_info/mydata_LPSgenesPA/after_step2_coordinates.bed"

# Create output directories
mkdir -p "$OUTPUT_DIR" "$SUBSET_FASTQ_DIR"

# Loop through 46 samples
#for i in $(seq 1 46); do

i=$SGE_TASK_ID    
# Extract sample name
    sample=$(sed -n "${i}p" "$SAMPLE_LIST")

    # Find BAM file
    BAM_FILE=$(find "$BAM_DIR" -maxdepth 1 -name "${sample}*.bam" | head -n 1)

    # Check if BAM file exists
    if [[ -z "$BAM_FILE" || ! -f "$BAM_FILE" ]]; then
        echo "ERROR: No BAM file found for sample ${sample}" >&2
        continue  # Skip this sample and move to the next
    fi

    echo "Processing BAM file: $BAM_FILE"

    # Ensure BAM is indexed
    samtools index "$BAM_FILE"

    # Find the raw FASTQ file
    RAW_FASTQ="${FASTQ_DIR}/${sample}.fastq.gz"

    # Check if the FASTQ file exists
    if [[ ! -f "$RAW_FASTQ" ]]; then
        echo "ERROR: No FASTQ file found for sample ${sample}" >&2
        continue
    fi

    # Read gene coordinates from BED file and extract read IDs
    while read chrom start end gene strand; do
        [[ "$chrom" == "#"* ]] && continue  # Skip header lines

        # Define read ID and subset FASTQ output files
        read_id_file="${OUTPUT_DIR}/${sample}_${gene}_readid.txt"
        subset_fastq="${SUBSET_FASTQ_DIR}/${sample}_${gene}_subset.fastq.gz"

        # Extract read IDs for the gene region
        samtools view "$BAM_FILE" "$chrom:$start-$end" | awk '{print $1}' | sort | uniq > "$read_id_file"

        echo "Extracted read IDs for gene: $gene in sample: $sample"

        # Extract subset FASTQ reads using seqtk
        /SAN/ugi/plant_genom/jiajucui/tools/seqtk/seqtk subseq "$RAW_FASTQ" "$read_id_file" | gzip > "$subset_fastq"

        echo "Subset FASTQ created: $subset_fastq"
    done < "$GENE_COORDS"


echo "Processing complete. Subset FASTQ files are stored in $SUBSET_FASTQ_DIR"
