#!/usr/bin/env python3

import pandas as pd
import os

# Directory where your HTSeq count files are stored
count_files_dir = '/mnt/Partition_2/Interns/Aly_Hissam/Quantification_files/'
count_files = os.listdir(count_files_dir)
combined_counts = pd.DataFrame()
for file in count_files:
    file_path = os.path.join(count_files_dir, file)
    sample_name = file.split('.')[0]
    counts = pd.read_csv(file_path, sep='\t', header=None, names=['GeneID', sample_name], index_col=0)
    # Remove any rows that aren't gene counts (e.g., "__no_feature", "__ambiguous", etc.)
    counts = counts[~counts.index.str.startswith("__")]
    # Add the counts to the combined DataFrame
    if combined_counts.empty:
        combined_counts = counts
    else:
        combined_counts = combined_counts.join(counts, how='outer')
# Save the combined counts to a CSV file
combined_counts.to_csv('Quantification_Combination/HTSeq_combined_counts.csv')
