# Load the edgeR package
library(edgeR)
counts <- read.csv("/Users/ALY HISSAM/OneDrive/Desktop/Magdi Intern/HTSeq_combined_counts.csv", header = TRUE, row.names = 1)
counts <- counts[, -ncol(counts)]
duplicates <- counts[duplicated(counts$GeneID), "GeneID"]
print(duplicates)
# Convert Null values to 0
counts[is.na(counts)] <- 0
# Assuming 'samples' is a vector containing the condition for each sample
groups <- factor(c("HCM", "AS", "HCM", "AS", "HCM", "HCM", "AS", "AS", "HCM", "HCM", "AS"))
# Convert the count matrix to a DGEList object
y <- DGEList(counts=counts, group=groups)
## normalize using TMM
y <- calcNormFactors(y)
## obtain normalized data on the log2-scale
logCPMs <- cpm(y, log = TRUE, prior.count = 1)
#Write normalised result to a file
write.table(logCPMs,file= "HTSeq_combined_counts_TMM_Normalized.tsv", sep="\t", row.names = TRUE,col.names = NA,quote = FALSE)
