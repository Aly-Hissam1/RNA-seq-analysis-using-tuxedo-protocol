# Load the edgeR package
library(edgeR)
counts <- read.table("/Users/ALY HISSAM/OneDrive/Desktop/Magdi Intern/R normalization/GSE206978_combined_with_symbols.tsv", header = TRUE, sep = "\t")

duplicates <- counts[duplicated(counts$Gene_Symbol), "Gene_Symbol"]
print(duplicates)
counts <- counts[!duplicated(counts$Gene_Symbol), ]
duplicates <- counts[duplicated(counts$Gene_Symbol), "Gene_Symbol"]

counts <- counts[, -1]
gene_symbols <- counts[, 1]
# Convert Null values to 0
counts[is.na(counts)] <- 0
# Assuming 'samples' is a vector containing the condition for each sample
groups <- factor(c("HCM", "HCM","HCM","HCM","HCM","HCM","HCM","HCM","AS", "AS","AS","AS","AS"))
# Convert the count matrix to a DGEList object
y <- DGEList(counts=counts, group=groups)
## normalize using TMM
y <- calcNormFactors(y)
## obtain normalized data on the log2-scale
logCPMs <- cpm(y, log = TRUE, prior.count = 1)
logCPMs <- cbind(Gene_Symbol = gene_symbols, logCPMs)

#Write normalised result to a file
write.table(logCPMs,file= "GSE206978_TMM_Normalized.tsv", sep="\t", row.names = TRUE,col.names = NA,quote = FALSE)

