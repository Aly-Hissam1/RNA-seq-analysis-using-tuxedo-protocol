library(biomaRt)

# Function to convert Ensembl IDs (with version) to gene symbols
convert_ensembl_to_symbol <- function(ensembl_ids) {
  # Connect to the Ensembl BioMart database
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  # Remove version numbers from Ensembl IDs (if present)
  ensembl_ids_no_version <- gsub("\\..*$", "", ensembl_ids)
  
  # Get gene symbols from Ensembl IDs
  genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id',
                 values = ensembl_ids_no_version,
                 mart = ensembl)
  
  # Return the mapping as a data frame
  return(genes)
}

# Load the TSV file (path needs to be updated to actual file location)
data <- read.table("/Users/ALY HISSAM/OneDrive/Desktop/Magdi Intern/R normalization/GSE206978_HCM_vs_stenosis_raw_counts.tsv", header = TRUE, sep = "\t")

# Assuming your dataframe is named 'df' and the column with Ensembl IDs is 'Ensembl_version'
# This line will loop through each value in 'Ensembl_version' and remove the version number after the dot
data$Ensembl_version <- gsub("\\..*$", "", data$Ensembl_version)

# Convert the first column from Ensembl IDs to gene symbols
ensembl_ids <- data[,1]
gene_symbols_df <- convert_ensembl_to_symbol(ensembl_ids)

# Merge the gene symbols back into the original dataframe
data_with_symbols <- merge(data, gene_symbols_df, by.x = 1, by.y = 'ensembl_gene_id', all.x = TRUE)

# Replace the first column with gene symbols
data_with_symbols[,1] <- data_with_symbols$external_gene_name

# Optionally, remove the now redundant 'external_gene_name' column
data_with_symbols$external_gene_name <- NULL

# Rename the first column to 'Gene_Symbol'
colnames(data_with_symbols)[1] <- 'Gene_Symbol'

duplicates <- data[duplicated(data$Ensembl_version), "Ensembl_version"]
data <- data[!duplicated(data$Ensembl_version), ]
duplicates <- data[duplicated(data$Ensembl_version), "Ensembl_version"]
print(duplicates)
write.table(data_with_symbols,file= "GSE206978_combined_with_symbols.tsv", sep="\t", row.names = TRUE,col.names = NA,quote = FALSE)

