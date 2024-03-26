if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "pathview"))
BiocManager::install("ReactomePA")
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)  # This is for human gene annotation; use the appropriate package for your organism
BiocManager::install('DESeq2')
library(edgeR)
library(DESeq2)


# Data preprocessing and normalization
count <- read.table("/Users/ALY HISSAM/OneDrive/Desktop/Magdi Intern/R normalization/HTSeq_combined_counts.csv", header=TRUE, sep=',')
count <- count[, -ncol(count)]
# Don't include 1st column since it is the gene names (or ids)
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
# Filter low-expressed genes
# Keep the genes that have Count-Per-Million more than k = 0.5 in n = 1 libraries
# It is pretty similar to filterByExpr(y, min.count = 0.5) from edgeR [1] but its choice of n is different.
logCPM <- cpm(y, log= TRUE)
plotDensities(logCPM, legend = FALSE,  main = "Before Filtering")
#Use function filterByExpr to go over count matrix and decide what are low counts
keep <- filterByExpr(y)
table(keep)
#Create a new object y that will have only the genes that are informative and do not have very low counts
y <- y[keep, , keep.lib.sizes=FALSE]
logCPM <- cpm(y, log= TRUE)
plotDensities(logCPM, legend = FALSE,  main = "After Filtering")
## normalize using TMM
y <- calcNormFactors(y, method = "TMM")
cpm.normalized <- cpm(y, log = FALSE, prior.count = 0, normalized.lib.sizes = TRUE)
write.csv(cpm.normalized, file = "/Users/ALY HISSAM/OneDrive/Desktop/Magdi Intern/R normalization/HTSeq_combined_counts_normalized.csv", row.names = TRUE)

normalized_matrix <- as.matrix(cpm.normalized)

annotation<-read.csv("/Users/ALY HISSAM/OneDrive/Desktop/Magdi Intern/DEGs/Annotation.csv", header=TRUE)
feature <- annotation$Condition
metadata <- data.frame(Condition = feature)
group <- as.character(metadata$Condition)
colData <- data.frame(sample = colnames(count), groups = group)
design <- model.matrix(~group, data = y) 
head(design)

##############################################################################
# limma-trend
library(limma)
# Perform the voom transformation and fit the linear model
fit <- lmFit(normalized_matrix, design, trend=TRUE)
# Apply empirical Bayes moderation
fit <- eBayes(fit)
# Extract the results and apply the thresholds
results <- topTable(fit, coef=2, sort.by="P", number=Inf, adjust.method="BH")
results$logFC <- as.numeric(results$logFC)
topGenestrend <- results[abs(results$logFC) > log2(1.0) & results$adj.P.Val < 0.25,]

#########################################################################
# Voom
# Perform the voom transformation and fit the linear model
v <- voom(normalized_matrix, design, plot=TRUE)
fit <- lmFit(v, design)
# Apply empirical Bayes moderation
fit <- eBayes(fit)
# Extract the results and apply the thresholds
results <- topTable(fit, coef=2, sort.by="P", number=Inf, adjust.method="BH")
results$logFC <- as.numeric(results$logFC)
topGenes <- results[abs(results$logFC) > log2(1.0) & results$adj.P.Val < 0.25,]


gene_names <- rownames(topGenestrend)

# Print gene names
print(gene_names)
#############################################################################3
# Pathway enrichment analysis

ge <- c('SSPN', 'ARHGAP24', 'PLEKHA7', 'BBS12', 'HDAC9', 'TNIK')
ge <- c('PLIN2', 'SSPN', 'LGR4', 'ARHGAP24', 'EZR', 'AMFR', 'TNIK', 'PLEKHA7', 'AGPAT1')


entrez_ids <- bitr(ge, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

entrez_ids
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





ego
ek
reactome_results
dotplot(ego)
dotplot(ek)
dotplot(reactome_results)

# For Reactome and other databases, similar functions are available in clusterProfiler or other packages.
