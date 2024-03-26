#!/bin/bash

To_Fastp()
{
	if [ -f "$cwd/Fastp_files/$ID""_R1.fil.fastq.gz" ] && [ -f "$cwd/Fastp_files/$ID""_R2.fil.fastq.gz" ] && [ -f "$cwd/Fastp_files/$ID.json" ] && [ -f "$cwd/Fastp_files/$ID.html" ];
	then
		echo "FastP file exists."
	else
		echo "Fastp the following files: "$forward, $reverse
		fastp -i $forward -I $reverse -o $cwd/Fastp_files/$ID"_R1.fil.fastq.gz" -O $cwd/Fastp_files/$ID"_R2.fil.fastq.gz"  --json $cwd/Fastp_files/$ID.json  --html $cwd/Fastp_files/$ID.html
		echo "Reads are filtered successfully"
	fi
}

TOPHAT_ALIGN()
{
        if [ -f "$Ref/BTgenome.2.bt2" ];
        then
                echo "Reference Indexed."
        else
                echo "Indexing the reference."
                bowtie2-build $ref BTgenome
        fi
        if [ -f "$cwd/Aligned_files/$ID.TOPHAT.ba/accepted_hits.bam" ];
        then
                echo "Bam file exists."
        else
		echo "Aligning the following files: "$ID"_R1.fil.fastq.gz", $ID"_R2.fil.fastq.gz"
		reads_1=$cwd/Fastp_files/$ID"_R1.fil.fastq.gz"
		reads_2=$cwd/Fastp_files/$ID"_R2.fil.fastq.gz"
		output=$cwd/Aligned_files/$ID.TOPHAT.ba
		transcript_annotation=/mnt/Partition_2/Interns/Aly_Hissam/genes.gtf
                tophat -p 4 -G ${transcript_annotation} -o ${output} --no-novel-juncs BTgenome ${reads_1} ${reads_2}
        fi
}

SortIndex()
{
        if [ -f "$cwd/Aligned_files/$ID.TOPHAT.ba/accepted_hits.sorted.bam" ];
        then
                echo "The processed file already indexed."
        else
                echo "Sorting the processed files"
		samtools sort $cwd/Aligned_files/$ID.TOPHAT.bam/accepted_hits.bam -@ 8 -o $cwd/Aligned_files/$ID.TOPHAT.ba/accepted_hits.sorted.bam
		echo "Indexing the processed files"
		samtools index $cwd/Aligned_files/$ID.TOPHAT.ba/accepted_hits.sorted.bam
        fi
}


Cufflinks()
{
	mkdir -p ${cwd}/cufflinks_output/${ID}_Cuff
        if [ -f ${cwd}/cufflinks_output/${ID}_Cuff/transcripts.gtf ];
        then
                echo "The cufflinks file already exist."
        else
                echo "Cufflinks is working on the processed files"
		cufflinks -p 8 -g $cwd/genes.gtf -o ${cwd}/cufflinks_output/${ID}_Cuff/ $cwd/Aligned_files/$ID.TOPHAT.ba/accepted_hits.sorted.bam
        fi
	# Path to the assembly file
	ASSEMBLY_FILE="${cwd}/assemblies.txt"

	# Path to the new file you want to add
	NEW_FILE_PATH="${cwd}/cufflinks_output/${ID}_Cuff/transcripts.gtf"

	# Create the assembly file if it does not exist
	if [ ! -f "$ASSEMBLY_FILE" ]; then
		touch "$ASSEMBLY_FILE"
		echo "Created $ASSEMBLY_FILE"
	fi

	# Check if the new file path is already in the assembly file
	if ! grep -Fxq "$NEW_FILE_PATH" "$ASSEMBLY_FILE"; then
		# If the path is not found, append it to the assembly file
		echo "$NEW_FILE_PATH" >> "$ASSEMBLY_FILE"
		echo "Added $NEW_FILE_PATH to $ASSEMBLY_FILE"
	else
		echo "$NEW_FILE_PATH is already in $ASSEMBLY_FILE"
	fi
}
HISAT2_ALIGN()
{
	if [ -f "$Ref/genome.6.ht2" ];
	then
		echo "Reference Indexed."
	else
		echo "Indexing the reference."
		hisat2-build $ref $Ref/genome
	fi
        if [ -f "$cwd/Aligned_files/$ID.bam" ];
        then
                echo "Bam file exists."
        else
                echo "Aligning the following files: "$ID"_R1.fil.fastq.gz", $ID"_R2.fil.fastq.gz"
                hisat2 -p 8 -x $Ref/genome -1 $cwd/Fastp_files/$ID"_R1.fil.fastq.gz" -2 $cwd/Fastp_files/$ID"_R2.fil.fastq.gz" | samtools sort -@ 8 -o $cwd/Aligned_files/$ID.bam
	fi
}

To_Index()
{
	if [ -f "$cwd/Aligned_files/"$ID".bam.bai" ];
	then
		echo "The processed file already indexed."
	else
		echo "Indexing the processed files"
		samtools index  "$cwd/Aligned_files/"$ID".bam"
	fi
}
HTSeq_QUANT(){
	if [ -f "$cwd/Quantification_files/$ID.txt" ];
	then
		echo "Quantification file exists."
	else
		echo "Quantifying the following files: "$ID".bam"
		htseq-count -f bam -s no -r pos -m union "$cwd/Aligned_files/"$ID".bam" /mnt/Partition_2/Interns/Aly_Hissam/genes.gtf > "$cwd/Quantification_files/$ID.txt"
	fi
}
cufmerge()
{
        if [ -f "$cwd/merged_asm/merged.gtf" ];
        then
                echo "The merged file already exist."
        else
                echo "Cuffmerge is running now..."
                cuffmerge -g ${cwd}/genes.gtf -s ${ref} -p 8 ${ASSEMBLY_FILE}
        fi
}
To_Combine()
{
	if [ -f "$cwd/Quantification_Combination/HTSeq_combined_counts.csv" ];
        then
                echo "The count matrix already exist."
        else
                echo "Combining quantification files..."
                python3 to_merge.py
        fi
}
Protocol1()
{
	To_Fastp
	TOPHAT_ALIGN
	SortIndex
	Cufflinks
}
Protocol2(){
	To_Fastp
	HISAT2_ALIGN
	To_Index
	HTSeq_QUANT
}
Call_Protocol()
{
	if [ "$call" == 1 ];
	then
		echo "Running First protocol TopHat on "$ID"..."
		Protocol1
	else
		echo "Running Second protocol HISAT2 on "$ID"..."
		Protocol2
	fi
}
filename=$1
while read line || [ -n "$line" ]; do
	IFS=$'\t' read -r forward reverse ref ID Script<<< "$line"
	Ref="$(dirname "$ref")"
	cwd=$(pwd)
	call=$Script
	mkdir -p $cwd/Fastp_files
	mkdir -p $cwd/Aligned_files
	mkdir -p $cwd/Quantification_files
	mkdir -p $cwd/Quantification_Combination
        mkdir -p $cwd/cufflinks_output
	Call_Protocol
	echo "$ID"
done < $filename
if [ "$call" == 1 ];
then
	cufmerge
	echo "Running First protocol TopHat on "$ID"..."
else
	To_Combine
fi
echo "We are Done..."

# Cuffdiff command:
# cuffdiff -o cuffdif_output -b /mnt/Partition_2/Interns/Aly_Hissam/genome.fa -p 8 -L HCM,AS -u /mnt/Partition_2/Interns/Aly_Hissam/merged_asm/merged.gtf /mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302632.TOPHAT.ba/accepted_hits.sorted.bam,/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302634.TOPHAT.ba/accepted_hits.sorted.bam,/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302635.TOPHAT.ba/accepted_hits.sorted.bam,/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302636.TOPHAT.ba/accepted_hits.sorted.bam,/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302640.TOPHAT.ba/accepted_hits.sorted.bam,/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302643.TOPHAT.ba/accepted_hits.sorted.bam /mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302631.TOPHAT.ba/accepted_hits.sorted.bam,/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302633.TOPHAT.ba/accepted_hits.sorted.bam,/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302639.TOPHAT.ba/accepted_hits.sorted.bam,/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302641.TOPHAT.ba/accepted_hits.sorted.bam,/mnt/Partition_2/Interns/Aly_Hissam/Aligned_files/SRR24302642.TOPHAT.ba/accepted_hits.sorted.bam
