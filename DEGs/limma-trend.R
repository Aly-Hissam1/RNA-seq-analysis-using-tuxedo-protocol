if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "pathview"))
BiocManager::install("ReactomePA")
BiocManager::install('DESeq2')
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)  # This is for human gene annotation; use the appropriate package for your organism
library(edgeR)
library(limma)
library(DESeq2)

# Data preprocessing and normalization
count <- read.table("/Users/ALY HISSAM/OneDrive/Desktop/Magdi Intern/R normalization/HTSeq_combined_counts.csv", header=TRUE, sep=',')
count <- count[, -ncol(count)]
# Remove duplication if exist...
count <- count [rev(order(apply(count[,-1], 1, sd))),]
count <- count [!duplicated(count[,1]),]
# Set row names to gene names (or IDs)
rownames(count) <- count[,1]
# Filter gene names
count <- count[,-1]
# Enforce all counts to be integers
count <- round(count, 0)
# Convert the count matrix to a DGEList object
y <- DGEList(counts=count)
logCPM <- cpm(y, log= TRUE)
plotDensities(logCPM, legend = FALSE,  main = "Before Filtering")
#Use function filterByExpr to go over count matrix and decide what are low counts
keep <- filterByExpr(y)
table(keep)
#Create a new object y that will have only the genes that are informative and do not have very low counts
y <- y[keep, , keep.lib.sizes=FALSE]
logCPM <- cpm(y, log= TRUE)
plotDensities(logCPM, legend = FALSE,  main = "After Filtering")
# Normalize using TMM
y <- calcNormFactors(y, method = "TMM")
cpm.normalized <- cpm(y, log = FALSE, prior.count = 0, normalized.lib.sizes = TRUE)
write.csv(cpm.normalized, file = "/Users/ALY HISSAM/OneDrive/Desktop/Magdi Intern/R normalization/HTSeq_combined_counts_normalized.csv", row.names = TRUE)
normalized_matrix <- as.matrix(cpm.normalized)
# Provide the label of the samples...
annotation<-read.csv("/Users/ALY HISSAM/OneDrive/Desktop/Magdi Intern/DEGs/Annotation.csv", header=TRUE)
feature <- annotation$Condition
metadata <- data.frame(Condition = feature)
group <- as.character(metadata$Condition)
colData <- data.frame(sample = colnames(count), groups = group)
design <- model.matrix(~group, data = y) 
head(design)
##############################################################################
# limma-trend
fit <- lmFit(normalized_matrix, design, trend=TRUE)
fit <- eBayes(fit)
# Extract the results and apply the thresholds
results <- topTable(fit, coef=2, sort.by="P", number=Inf, adjust.method="BH")
results$logFC <- as.numeric(results$logFC)
topGenestrend <- results[abs(results$logFC) > log2(1.0) & results$adj.P.Val < 0.25,]
gene_names <- rownames(topGenestrend)
# Print gene names
print(gene_names)
#############################################################################3
# Pathway enrichment analysis
# Update ge variable with the intersected genes between different tools...
ge <- c('SSPN', 'ARHGAP24', 'PLEKHA7', 'BBS12', 'HDAC9', 'TNIK')
entrez_ids <- bitr(ge, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
ego <- enrichGO(gene         = entrez_ids$ENTREZID,
                OrgDb        = org.Hs.eg.db,
                keyType      = "ENTREZID",
                ont          = "BP",  # For Biological Processes; use "CC" for Cellular Components or "MF" for Molecular Functions
                pAdjustMethod = "BH",  # Method for adjusting p-values
                qvalueCutoff  = 0.2,  # Threshold for adjusted p-values
                readable      = TRUE)  # To get readable gene names
ek <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                 organism     = "hsa",  
                 pAdjustMethod = "BH",
                 qvalueCutoff  = 0.05)
reactome_results <- enrichPathway(gene = entrez_ids$ENTREZID,
                                  organism = "human", # use "human" for human data
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff = 0.05)
dotplot(ego)
dotplot(ek)
dotplot(reactome_results)
