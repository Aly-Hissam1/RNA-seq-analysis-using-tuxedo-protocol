#!/bin/bash

# Array of BAM file paths
BAM_FILES=(
"/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302632.TOPHAT.ba/accepted_hits.sorted.bam"
"/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302634.TOPHAT.ba/accepted_hits.sorted.bam"
"/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302635.TOPHAT.ba/accepted_hits.sorted.bam"
"/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302636.TOPHAT.ba/accepted_hits.sorted.bam"
"/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302640.TOPHAT.ba/accepted_hits.sorted.bam"
"/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302643.TOPHAT.ba/accepted_hits.sorted.bam"
"/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302631.TOPHAT.ba/accepted_hits.sorted.bam"
"/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302633.TOPHAT.ba/accepted_hits.sorted.bam"
"/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302639.TOPHAT.ba/accepted_hits.sorted.bam"
"/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302641.TOPHAT.ba/accepted_hits.sorted.bam"
"/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302642.TOPHAT.ba/accepted_hits.sorted.bam"
)

# Output directory for average coverage files
OUTPUT_DIR="/mnt/Partition_2/Interns/Aly_Hissam/Coverage_results/"
mkdir -p "$OUTPUT_DIR"

# Loop through each BAM file
# Loop through each BAM file
for BAM_FILE in "${BAM_FILES[@]}"; do
    # Get the directory name of the BAM file
    DIR_NAME=$(dirname "$BAM_FILE")
    # Extract the name of the directory, not the full path
    BASE_NAME=$(basename "$DIR_NAME")

    # Define output file name for average coverage, using the directory name
    AVG_COV_FILE="$OUTPUT_DIR/${BASE_NAME}.txt"

    # Calculate average coverage using samtools depth and awk, and write to the output file
    samtools depth "$BAM_FILE" | awk '{sum+=$3} END {print "Average coverage for '"$BASE_NAME"':",sum/NR}' > "$AVG_COV_FILE"
    echo "Average coverage calculated for $BAM_FILE, output in $AVG_COV_FILE"
done

echo "All average coverage calculations are complete."

