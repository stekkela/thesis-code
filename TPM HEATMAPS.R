# Read CSV (first column = gene names, rest = samples)
fpkm_df <- read.csv("~/RNA-seqAdipo/pair1genes.csv", check.names = FALSE)
head(fpkm_df)
library(dyprl)
# Suppose your first column is gene names
expr <- fpkm_df[, -1]   # remove gene_name column

# Make sure it's numeric
expr <- as.data.frame(lapply(expr, as.numeric))

# Now calculate TPM
tpm <- apply(expr, 2, function(x) (x / sum(x)) * 1e6)
rownames(tpm) <- fpkm_df$gene_name
library(pheatmap)
tpm


# 1️⃣ Remove genes with very low expression in all samples (optional)
tpm_filtered <- tpm[rowMeans(tpm) > 1, ]   # keep genes with mean TPM > 1

# 2️⃣ Compute mean TPM per sample
sample_means <- colMeans(tpm_filtered)

# 3️⃣ Compute overall mean
overall_mean <- mean(sample_means)

# 4️⃣ Keep samples that are not below the overall mean
tpm_final <- tpm_filtered[, sample_means >= overall_mean]

adipo_FACE_GENES <- c("AKRIN2", "CTSB", "SMILR", "GNG2", "MYOM1", "EIF3FP3", 
                 "SMIM10", "ELFN1", "H3F3AP4", "MGARP", "ABCC4", "KRT7", "PDE7B")

ppargactivators <- c("GOLGA8A","PHACTR3", "MEOX2","POU3F3", "MTIG","COL4A3","CXCL14","NRXN3","ADRB1", "FABP3")


WNT <- c( "DAB2","DISC1","NKD1","RSPO1","RSPO3","TLE1","TRABD2B","WNT11","WNT9A","CALCOCO1","DAKK1","DAAM2","DACT1", "DACT3","DACT3","TCF7","ZBED3")
ecm<- c("ADAMTS15","ADAMTSL5","SMOC2","TNFRF11B","COL4A6","COL8A2","MMP1","MMP15","MMP19","MMP27","MMP3")

  
TF <- c("ARDI5B","ELF4","ATF5","ARNT2","EPC2","IRF7","MORF4L2","PPARG","PARP9","STAT1","SPIN1","TP53","TRIM8","ETV7","SMURF2","BCL6","AEBP1","ARID5B","WWP2","CBX7","HIPK2","NFKB1","NR1D2","PER1","STAT1","ZBTB2","ZBTB14","ZEB2")
  
STAT <- c("LIF","CLCF1","CSF2","CRLF1","IGF1","IL6")
 
  
  subset_fpkm <- tpm_final %>%
  filter(rownames(tpm_final) %in% ppargactivators)
subset_fpkm <- fpkm_df %>%
  filter(gene_name %in% ppargactivators) %>%
  column_to_rownames("gene_name")

log_fpkm <- log2(subset_fpkm + 1)

pheatmap(subset_fpkm,
         cluster_rows = TRUE,      # cluster genes
         cluster_cols = FALSE,     # keep samples in original order
         scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "STAT")

log_tpm <- log2(subset_tpm + 1)



# Plot heatmap
pheatmap(log_tpm,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         main = "Expression of EBF2, ZNF423, PPARG, and CEBPA (TPM)")




###Do the rest


# 5️⃣ Save the filtered TPMs
write.csv(tpm_final, "filtered_TPMs.csv")
log_tpm <- log2(tpm + 1)

var_genes <- head(order(apply(log_tpm, 1, var), decreasing=TRUE), 50)
mat_top <- log_tpm[var_genes, ]
mat_top

# Heatmap
pheatmap(mat_top,
         scale = "row",             # standardize expression per gene
         fontsize_row = 8,
         fontsize_col = 10,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         angle_col = 45)
install.packages("pheatmap")  
library(pheatmap)

# Suppose mat_subset is your genes x samples matrix
pheatmap(mat_top,
         scale = "row",
         fontsize_row = 8,
         fontsize_col = 10,
         cluster_rows = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         angle_col = "45") # quotes! must be a string


pheatmap(mat_top,
         scale = "row",
         fontsize_row = 8,
         fontsize_col = 10,
         cluster_rows = TRUE,       # keep clustering for rows
         cluster_cols = FALSE,      # disable clustering for columns
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         angle_col = "45")













subset_fpkm <- tpm[rownames(tpm) %in% STAT, ]




# Example: 'tpm' is a data frame of TPMs with gene names as rownames
# and columns: AN1ADAD_rep1, AN1ADAD_rep2, AK1DAD_rep1, AK1DAD_rep2, etc.

# First, subset samples for each group
group1 <- grep("AN1ADAD", colnames(tpm), value = TRUE)
group2 <- grep("AK1DAD", colnames(tpm), value = TRUE)
grep("AN1ADAD", colnames(tpm), value = TRUE)
tpm_group1 <- tpm[, group1]
tpm_group2 <- tpm[, group2]



# Calculate variance between the two groups (e.g. absolute log2 fold change)
logFC <- log2(tpm_group1) - log2(tpm_group2 + 1)

# Or if you want actual variability measure:
# compute variance of expression values across all replicates in both conditions
var_total <- apply(tpm[, c(tpm_group1, group2)], 1, var)

# Combine into a data frame
var_genes <- data.frame(
  gene = rownames(tpm),
  logFC = logFC,
  variance = var_total
)

# Rank by variability or absolute fold change
var_genes <- var_genes[order(-abs(var_genes$logFC)), ]

# Top 100 most variable genes
head(var_genes, 100)



library(pheatmap)
varvar_genes_filtered
tpm_filtered <- tpm[!grepl("^(MT|RNU|RIF|RF|MIR|SNOR|AC|AP|AL|RP|RN|RA|RX)|RNA", rownames(tpm), ignore.case = TRUE), ]
tpm_filtered <- tpm_filtered[rowSums(tpm_filtered[, c("AN1ADAD", "AK1ADAD")] >= 1) > 0, ]
logFC <- log2(tpm_filtered[, "AN1ADAD"] + 1) - log2(tpm_filtered[, "AK1ADAD"] + 1)

var_genes_filtered <- data.frame(
  gene = rownames(tpm_filtered),
  logFC = logFC
)
var_genes_filtered <- var_genes_filtered[order(-abs(var_genes_filtered$logFC)), ]
top_genes <- head(var_genes_filtered$gene, 50)
top_genes <- head(var_genes_filtered$gene, 20)
top_genes
# Step 6: Prepare heatmap data
heatmap_data <- log2(tpm_filtered[top_genes, c("AN1ADAD", "AK1ADAD")] + 1)

# Step 7: Draw heatmap

pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  main = "Top 20 Most Different Genes (RNA/RNU/RIF filtered)",fontsize_row = 5,   # smaller row labels (genes)
  fontsize_col = 10
)
# 1️⃣ Compute log2 fold change (if not already done)
logFC <- log2(tpm_filtered[, "AN1ADAD"] + 1) - log2(tpm_filtered[, "AK1ADAD"] + 1)

# 2️⃣ Keep only genes higher in AK1DAD (logFC < 0)
genes_higher_in_AK1DAD <- names(logFC)[logFC < 0]

# 3️⃣ Take top 20 most differentially higher genes
top_genes_AK1DAD <- head(genes_higher_in_AK1DAD[order(logFC[genes_higher_in_AK1DAD])], 50)

# 4️⃣ Subset TPM for heatmap
heatmap_data_AK1DAD <- log2(tpm_filtered[top_genes_AK1DAD, c("AN1ADAD", "AK1ADAD")] + 1)

# 5️⃣ Draw heatmap
library(pheatmap)
pheatmap(
  heatmap_data_AK1DAD,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "Top 20 Genes Higher in AK1DAD")



# 1️⃣ Filter lowly expressed genes (TPM ≥ 1 in at least one sample)
tpm_filtered <- tpm[rowSums(tpm[, c("AN1ADUN", "AK1ADUN")] >= 1) > 0, ]
tpm_filtered <- tpm[!grepl("^(MT|RNU|RIF|RF|MIR|SNOR|AC|AP|AL|RP|RN|RA|RX)|RNA", rownames(tpm), ignore.case = TRUE), ]
tpm_filtered2 <- tpm_filtered[rowSums(tpm_filtered[, c("AN1ADUN", "AK1ADUN")] >= 1) > 0, ]

# 2️⃣ Compute log2 fold change for the filtered genes
logFC <- log2(tpm_filtered2[, "AN1ADUN"] + 1) - log2(tpm_filtered[, "AK1ADUN"] + 1)

# 3️⃣ Create data frame
var_genes <- data.frame(
  gene = rownames(tpm_filtered2),
  logFC = logFC
)

# 4️⃣ Sort by absolute logFC
var_genes_sorted <- var_genes[order(-abs(var_genes$logFC)), ]
colnames(tpm_filtered2)

# 5️⃣ Take top 100 variable genes
top100_genes <- head(var_genes_sorted, 50)
Top_genes_vector <- top100_genes$gene

# Subset TPM for heatmap
heatmap_data <- log2(tpm_filtered2[Top_genes_vector, c("AN1ADUN", "AK1ADUN")] + 1)

# Draw heatmap
library(pheatmap)
pheatmap(
  heatmap_data,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "Top 100 Variable Genes"
)

heatmap_data_AK1DAD <- log2(tpm_filtered[top100_genes, c("AN1ADUN", "AK1ADUN")] + 1)

# 5️⃣ Draw heatmap
library(pheatmap)
pheatmap(
  heatmap_data_AK1DAD,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "NORMAL CELLS VARIABLE GENES"
)
