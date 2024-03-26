#!/bin/bash

# Directory containing BAM files
BAM_DIR="/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/"

# Output directory for coverage files and average coverage
OUTPUT_DIR="/mnt/Partition_2/Interns/Aly_Hissam/Coverage_hisat_files/"
mkdir -p "$OUTPUT_DIR"

# Loop through each BAM file in the directory
for BAM_FILE in "$BAM_DIR"/*.bam; do
    # Extract the base name of the file (without path and extension)
    BASE_NAME=$(basename "$BAM_FILE" .bam)

    # Define output file name for coverage
    OUTPUT_FILE="$OUTPUT_DIR/${BASE_NAME}_coverage.txt"
    # Define output file name for average coverage
    AVG_COV_FILE="$OUTPUT_DIR/${BASE_NAME}_average_coverage.txt"

    # Run samtools depth and output the coverage to a file
    samtools depth "$BAM_FILE" > "$OUTPUT_FILE"
    echo "Coverage for $BAM_FILE written to $OUTPUT_FILE"

    # Calculate average coverage using awk and output to a separate file
    awk '{sum+=$3} END {print "Average coverage: ",sum/NR}' "$OUTPUT_FILE" > "$AVG_COV_FILE"
    echo "Average coverage for $BAM_FILE written to $AVG_COV_FILE"
done

