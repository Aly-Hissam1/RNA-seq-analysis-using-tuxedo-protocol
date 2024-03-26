install.packages("DESeq2")
library(dplyr)
library(limma)
library(edgeR)
library(paletteer)


count.matrix <- read.csv("/Users/ALY HISSAM/OneDrive/Desktop/Magdi Intern/HTSeq_combined_counts.csv")

count.matrix <- count.matrix [rev(order(apply(count.matrix[,-1], 1, sd))),]
count.matrix <- count.matrix [!duplicated(count.matrix[,1]),]

# Set row names to gene names (or IDs)
rownames(count.matrix) <- count.matrix[,1]

# Filter gene names
count.matrix <- count.matrix[,-1]

# Enforce all counts to be integers
count.matrix <- round(count.matrix, 0)

# Filter low-expressed genes
# Keep the genes that have Count-Per-Million more than k = 0.5 in n = 1 libraries
# It is pretty similar to filterByExpr(y, min.count = 0.5) from edgeR [1] but its choice of n is different.
count.matrix <- count.matrix [rowSums(cpm(count.matrix) >= 0.5) >= 1,]


# Example meta.data object
meta.data <- data.frame(
  SampleID = c("SRR24302636",	"SRR24302631",	"SRR24302634",	"SRR24302641",	"SRR24302632",	"SRR24302643",	"SRR24302642",	"SRR24302639",	"SRR24302635",	"SRR24302640",	"SRR24302633"),
  Status = c("HCM", "AS", "HCM", "AS", "HCM", "HCM", "AS", "AS", "HCM", "HCM", "AS")
)


# Create your desired groups
# Remember these parameters (Status, Age, Gender) depends on your own analyses!
group <- meta.data$Status

# Assign each sample to its group assuming the column names of 
# count matrix are the sample names.
colData <- data.frame(sample = colnames(count.matrix), groups = group)
colnames(colData)  = c("sample", "groups")





# Create a DGEList object with just the count data
dge <- DGEList(counts = count.matrix)

# Add group information
dge$samples$group <- factor(colData$groups)




# Perform TMM normalization
dge <- calcNormFactors(dge, method = "TMM")

# Get the normalized counts
count.matrix.normalized <- cpm(dge, log = TRUE)




color.palette <- paletteer_d("ggsci::nrc_npg")



barplot(colSums(count.matrix)/1e6, 
        las  = 3, 
        col  = sample(color.palette, size = 1), 
        main ="Total read counts (millions)")  

barplot(colSums(count.matrix.normalized)/1e6, 
        las  = 3, 
        col  = sample(color.palette, size = 1), 
        main ="Total read counts (millions)")  

#####################################################

plot(density(count.matrix.normalized[,1]), 
     lwd  = 2,
     col  = sample(color.palette, size = 1),
     xlab = "Expression values", ylab="Density", 
     main = "Distribution of transformed data")
for (i in seq(2,ncol(count.matrix.normalized))){
  lines(density(count.matrix.normalized[,i]), 
        lwd = 2, 
        col = sample(color.palette, size = 1))
}

#####################################################


FDR.cutoff <- 0.1
LFC.cutoff <- 0


# edgeR
# Ensure group is a factor with correct levels
group <- factor(group)
levels(group) # Check the levels




# Create DGEList object
y <- DGEList(counts = count.matrix)

# Calculate normalization factors for library sizes with TMM
y <- calcNormFactors(y, method = "TMM")

# Add intercept term for multiple comparisons
design <- model.matrix(~ 0 + group) 
rownames(design) <- colnames(count.matrix)
colnames(design) <- levels(group)

# Estimate dispersion for genes with Bayesian Shrinkage
y <- estimateDisp(y,design)

# Fit the model
fit.glm <- glmQLFit(y,design)

# Since you have two conditions: "HCM" and "AS", let's create a contrast to compare these two
contrasts <- makeContrasts(
  HCMvsAS = HCM - AS,
  levels = design
)

# Fit the model with the contrast
fit.contrast <- glmQLFTest(fit.glm, contrast=contrasts[,"HCMvsAS"])

# Extract the top differentially expressed genes (DEGs)
top.degs <- topTags(fit.contrast)

# View the top DEGs
top.degs
###############################################

# Limma trend
library(limma)

fit.trend <- lmFit(count.matrix.normalized, design)
fit.trend2 <- eBayes(contrasts.fit(fit.trend, contrasts), trend = T)

# Create DE gene list for limma-trend
DE.genes.trend <- list()

for (i in seq_len(ncol(contrasts))){
  contrast.name <- colnames(contrasts)[i]
  top.table <- topTable(fit.trend2, 
                        coef = colnames(contrasts)[i],
                        p.value = FDR.cutoff,
                        lfc = LFC.cutoff,
                        number = Inf)
  DE.genes.trend[[contrast.name]] <- top.table
}


# Assuming the contrast of interest is named "HCMvsAS"
top.degs.limma.trend <- topTable(fit.trend2, coef="HCMvsAS", sort.by="P", number=10)

# Print the top differentially expressed genes
print(top.degs.limma.trend)

#################################################
#limma voom
v <- voom(count.matrix, design, plot=F)
fit.voom <- lmFit(v, design)
fit.voom2 <- eBayes(contrasts.fit(fit.voom, contrasts))
# summary(decideTests(fit.voom2, method="separate", lfc = 0, p.value = 0.1))

# Create DE gene list for limma-voom
DE.genes.voom <- list()

for (i in seq_len(ncol(contrasts))){
  contrast.name <- colnames(contrasts)[i]
  top.table <- topTable(fit.voom2, 
                        coef = colnames(contrasts)[i],
                        p.value = FDR.cutoff, 
                        lfc = LFC.cutoff, 
                        number = Inf)
  DE.genes.voom[[contrast.name]] <- top.table
}


# Assuming the contrast of interest is named "HCMvsAS"
top.degs.limma.voom <- topTable(fit.voom2, coef="HCMvsAS", sort.by="P", number=10)

# Print the top differentially expressed genes
print(top.degs.limma.voom)
############################################
# Deseq2
library( "DESeq2" )

dds <- DESeqDataSetFromMatrix(countData = count.matrix, colData = colData, design = ~ groups)


dge <- DGEList(counts=counts)

keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]


dge <- calcNormFactors(dge)

logCPM <- cpm(dge, log=TRUE, prior.count=3)

fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
topTable(fit, coef=ncol(design))

# Load the updated count data from the uploaded file
file_path_updated <- "/Users/ALY HISSAM/OneDrive/Desktop/Magdi Intern/HTSeq_combined_counts.csv" # Update the file path accordingly
counts_updated <- read.csv(file_path_updated, row.names = 1)

# Extract labels from the last row and then remove it from the count data
labels_updated <- as.character(tail(counts_updated, n=1))
counts_updated <- counts_updated[-nrow(counts_updated), ]

# Convert count data to numeric as they might be read as characters
counts_updated <- apply(counts_updated, 2, as.numeric)

# Display the first few rows of the updated count data and the labels to ensure they are correctly extracted
head(counts_updated)
labels_updated



library(limma)
library(edgeR)

# Assuming 'counts_updated' is your count matrix and 'labels_updated' contains your sample labels

# 1. Create DGEList object
dge <- DGEList(counts=counts_updated)

# 2. Filter out lowly expressed genes

# Assuming 'labels_updated' contains your sample labels
# Create a design matrix
group <- factor(labels_updated)
design <- model.matrix(~0+group)  # '0+' omits the intercept, creating a binary design matrix

# Create a DGEList object
dge <- DGEList(counts=counts_updated)

# Filter lowly expressed genes with the design matrix
keep <- filterByExpr(dge, design=design)
dge <- dge[keep,]


# 3. Calculate normalization factors
dge <- calcNormFactors(dge)

# 4. Create the design matrix
design <- model.matrix(~0 + labels_updated)
colnames(design) <- levels(labels_updated)

# 5. Fit the linear model and apply eBayes
fit <- lmFit(dge, design)
fit <- eBayes(fit)

# 6. Extract top differentially expressed genes
topGenes <- topTable(fit, coef=2)  # Assuming you want to compare the second condition to the first

# View top differentially expressed genes
head(topGenes)
