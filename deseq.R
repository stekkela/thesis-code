library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(viper)
library(DESeq2)

BiocManager::install("viper")

data <- read.csv("C:/Users/hmz255/OneDrive - Queen Mary, University of London/Documents/RNA-seqAdipo/gene_count.csv", header=T)
data


data <- data[!duplicated(data$gene_name),]
data <- data[!is.na(data$gene_name),]

fibsdata <- read.csv("C:/Users/hmz255/OneDrive - Queen Mary, University of London/Documents/RNA-seqAdipo/DESeq2_results.csv", header=T)
data


rownames(fibsdata)<- data$gene_name

fibsdata$gene_name <- NULL
head(fibsdata)
header <- read.csv("C:/Users/hmz255/OneDrive - Queen Mary, University of London/Documents/RNA-seqAdipo/metadata.csv", header=TRUE)

header

dds <- DESeqDataSetFromMatrix(
  countData = data,
  colData = header,
  design = ~ 1      # design is irrelevant for normalization
)

vsd <- vst(dds)

#########DO ALL THE DIFFERENT DDS 
# Subset face adipocytes
dds_adipocyte<- dds[,dds$media == "adipocyte"]

# Drop unused levels
dds_adipocyte$body_site <- droplevels(dds_adipocyte$body_site)
dds_adipocyte$media     <- droplevels(dds_adipocyte$media)
dds_adipocyte$cell_type <- droplevels(dds_adipocyte$cell_type)
# Make sure 'normal' is the reference
dds_adipocyte$cell_type <- relevel(dds_adipocyte$cell_type, "normal")

design(dds_adipocyte) <- ~ cell_type
dds_adipocyte <- DESeq(dds_adipocyte)
res_adipocyte <- results(dds_adipocyte, contrast = c("cell_type", "keloid", "normal"))
dds_adipocyte
levels(dds_adipocyte$cell_type)

ds
resultsNames(dds_adipocyte)
keep <- rowSums(counts(dds_adipocyte) >= 10) >= 2
dds <- dds_adipocyte[keep, ]

# 3. Run DESeq
dds_adipocyte <- DESeq(dds_adipocyte)

# 4. Get results
res <- results(dds_adipocyte)

dds <- dds[keep, ]

# 3. Apply rlog
rld <- rlog(dds, blind = FALSE)

##################################### PCA


vsd <- vst(fibsdata, blind=TRUE) 
mat <- assay(dds)
mat <- as.matrix(fibsdata)
mat
# assume net has columns: source, target, weight
res <- decouple(
  mat = mat,
  net = regulons,
  .source = "source",
  .target = "target")

res <- decouple(
  mat = mat,
  net = regulons,
  .source = "source",
  .target = "target",
  minsize = 1
)
f# Load DoRothEA human regulons
regulons <- decoupleR::get_dorothea(organism = "human", levels = c("A", "B", "C"))

regulons

write.csv(res,"~/REAADIPO.csv")

res1 <- run_viper(
  mat = expr_mat,
  net = regulons,
  .source = "source",
  .target = "target",
  .mor = "mor"
)

res2 <- decouple(
  mat = mat,          # your VST or rlog matrix
  net = regulons,     # your TF→target network
  method = "mlm"      # <- specify the method explicitly
)

sig_tfs <- res %>%
  filter(p_value < 0.05, abs(score) > 3)

res

pheatmap::pheatmap(mat = sig_tfs,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 11,
                   cellheight = 11,
                   treeheight_row = 11,
                   treeheight_col = 11)


vars <- rowVars(mat)

# Keep top 1000 variable genes
top_genes <- names(sort(vars, decreasing = TRUE))[1:1000]

pheatmap(mat[top_genes, ])

res_gsea <- decoupleR::run_fgsea(mat = mat, 
                                 network = network, 
                                 .source = 'source', 
                                 .target = 'target', 
                                 nproc = 1, 
                                 minsize = 0)

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)
pheatmap::pheatmap(mat = mat,
                   color = colors.use,
                   border_color = "white",
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 0,
                   treeheight_col = 0)



pheatmap(sig_tfs,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE)




sig_numeric <- sig_tfs[, sapply(sig_tfs, is.numeric)]
sig_numeric <- sig_numeric[complete.cases(sig_numeric), ]



sig_only <- sig_tfs %>%
  filter(p_value < 0.01, abs(score) > 5)
nrow(sig_only)
mat_sig <- sig_only %>%
  group_by(source, condition) %>%
  summarize(score = mean(score, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()
mat_sig
mat_sig[is.na(mat_sig)] <- 0

row_vars <- apply(mat_sig, 1, var)
mat_sig <- mat_sig[row_vars > 0, ]

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

pheatmap::pheatmap(
  mat = mat_sig,
  color = colors.use,
  border_color = "white",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 15,
  cellheight = 15,
  treeheight_row = 0,
  treeheight_col = 0
)




top_tfs <- sig_tfs %>%
  group_by(source) %>%
  summarize(max_abs_score = max(abs(score))) %>%
  arrange(desc(max_abs_score)) %>%
  slice(1:30)  








mat

mat <- as.matrix(sig_numeric)
pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = TRUE)



library(dplyr)




sig_tfs <- res %>%
  filter(p_value < 0.05, abs(score) > 2)

library(tidyr)
library(tibble)

heatmap_mat <- sig_tfs %>%
  filter(source %in% sig_tfs$source) %>%
  pivot_wider(names_from = condition, values_from = statistic) %>%
  column_to_rownames("source") %>%
  as.matrix()


heatmap_mat <- sig_tfs %>%
  group_by(source, condition) %>%   # aggregate if duplicates
  summarize(statistic = mean(statistic, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = condition, values_from = statistic) %>%
  column_to_rownames("source") %>%
  as.matrix()



# Convert expression matrix to long format
mat_long <- as.data.frame(mat) %>%
  rownames_to_column("target") %>%
  pivot_longer(-target, names_to = "sample", values_to = "expr")


# Filter significant TFs (optional, if you have thresholds)
# Here just using activity > 2 or < -2 as example
sig_tf_activity <- tf_activity %>%
  filter(abs(activity) > 2)

heatmap_mat <- sig_tf_activity %>%
  pivot_wider(names_from = sample, values_from = activity) %>%
  column_to_rownames("source") %>%
  as.matrix()
heatmap_mat

library(pheatmap)
pheatmap(heatmap_mat, scale = "row", cluster_rows = TRUE, cluster_cols = TRUE)

heatmap_mat_clean <- heatmap_mat[complete.cases(heatmap_mat), ]
heatmap_mat[is.na(heatmap_mat)] <- 0
heatmap_mat_clean <- heatmap_mat
pheatmap(
  heatmap_mat_clean,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Significant TF activity"
)


heatmap_mat_numeric <- apply(heatmap_mat, 2, as.numeric)
rownames(heatmap_mat_numeric) <- rownames(heatmap_mat)


colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)
heatmap_mat_numeric[is.na(heatmap_mat_numeric)] <- 0

heatmap_mat_clean

pheatmap(
  heatmap_mat_numeric,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Significant TF activity"
)

pheatmap::pheatmap(mat = heatmap_mat,
                   color = colors.use,
                   border_color = "white",
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   cellwidth = 11,
                   cellheight = 11,
                   treeheight_row = 3,
                   treeheight_col = 3)


# Heatmap
pheatmap::pheatmap(mat = mat,
                   color = colors.use,
                   border_color = "white",
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   cellwidth = 15,
                   cellheight = 15,
                   treeheight_row = 0,
                   treeheight_col = 0)


str(res$statistic)

res$statistic <- as.numeric(as.character(res$statistic))




activity_mat <- res %>%
  select(source, condition, score) %>%
  pivot_wider(names_from = condition, values_from = score) %>%
  column_to_rownames("source") %>%
  as.matrix()
var_genes <- apply(expr_mat, 1, function(x) sd(x, na.rm = TRUE))
deg <- names(sort(var_genes, decreasing = TRUE))[1:50]  # top 50 variable

activity_deg <- activity_mat[deg, ]

pheatmap(
  activity_deg,
  scale = "row",
  cluster_cols = FALSE, cluster_rows = TRUE)

pheatmap(
  heatmap_mat_clean,
  scale = "row",
  clustering_method = "ward.D2",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 6
)


# Transform to wide matrix
sample_acts_mat <- expr_mat %>%
  tidyr::pivot_wider(id_cols = 'condition', 
                     names_from = 'source',
                     values_from = 'score') %>%
  tibble::column_to_rownames('condition') %>%
  as.matrix()

# Scale per feature
sample_acts_mat <- scale(sample_acts_mat)

# Color scale
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05,2, length.out = floor(100 / 2)))

# Plot
pheatmap::pheatmap(mat = sample_acts_mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 2,
                   cellheight = 2,
                   treeheight_row = 2,
                   treeheight_col = 2)


# Calculate variance per feature
feature_var <- apply(sample_acts_mat, 2, var)

# Keep only features with high variance
threshold <- 0.5  # adjust based on your dataset
diff_features <- names(feature_var[feature_var > threshold])

# Subset matrix
sample_acts_mat_diff <- sample_acts_mat[, diff_features]

# Calculate variance per feature
feature_var <- apply(sample_acts_mat, 2, var)

pheatmap::pheatmap(mat = sample_acts_mat_diff,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 2,
                   cellheight = 2,
                   treeheight_row = 2,
                   treeheight_col = 2)

regulons






library(dplyr)

library(dplyr)
library(tidyr)

# Keep only 'condition' and numeric columns
numeric_cols <- c("condition", names(sample_acts)[sapply(sample_acts, is.numeric)])

sample_long <- sample_acts %>%
  select(all_of(numeric_cols)) %>%
  pivot_longer(cols = -condition, names_to = "source", values_to = "score")


# Initialize a vector to store p-values
pvals <- sapply(unique(sample_long$source), function(feature) {
  scores <- sample_long %>% filter(source == feature)
  aov_res <- aov(score ~ condition, data = scores)
  summary(aov_res)[[1]][["Pr(>F)"]][1]
})

# Select features with significant difference
diff_features <- names(pvals[pvals < 0.05])
sample_acts_mat_diff <- sample_acts_mat[, diff_features, drop = FALSE]
# Remove columns with any NA/NaN/Inf
sample_acts_mat_diff <- sample_acts_mat_diff[, 
                                             apply(sample_acts_mat_diff, 2, function(x) all(is.finite(x)))]


sample_acts_mat_diff <- sample_acts_mat_diff[
  apply(sample_acts_mat_diff, 1, function(x) all(is.finite(x))),
]
pheatmap::pheatmap(mat = sample_acts_mat_diff,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 2,
                   cellheight = 2,
                   treeheight_row = 2,
                   treeheight_col = 2)








sample_acts <- decoupleR::run_mlm(mat = data, 
                                  net = regulons, 
                                  .source = 'source', 
                                  .target = 'target',
                                  .mor = 'mor', 
                                  minsize = 5)


pheatmap(
  heatmap_mat_clean,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Significant TF activity"
)





# Transform to wide matrix
sample_acts_mat <- sample_acts %>%
  tidyr::pivot_wider(id_cols = 'condition', 
                     names_from = 'source',
                     values_from = 'score') %>%
  tibble::column_to_rownames('condition') %>%
  as.matrix()

any(is.nan(heatmap_filtered))
any(is.infinite(heatmap_filtered))

heatmap_filtered[is.nan(heatmap_filtered)] <- 0
heatmap_filtered[is.infinite(heatmap_filtered)] <- 0


pheatmap(
  heatmap_filtered,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Significant TF activity"
)

# Scale per feature
sample_acts_mat <- scale(sample_acts_mat)

# Color scale
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
colors.use <- grDevices::colorRampPalette(colors = colors)(100)

my_breaks <- c(seq(-2, 0, length.out = ceiling(100 / 2) + 1),
               seq(0.05,2, length.out = floor(100 / 2)))

# Plot
pheatmap::pheatmap(mat = mat,
                   color = colors.use,
                   border_color = "white",
                   breaks = my_breaks,
                   cellwidth = 11,
                   cellheight = 11,
                   treeheight_row = 11,
                   treeheight_col = 11)



write.csv(sample_acts_mat,"TFpathways.csv")
library(ggplot2)
library(RColorBrewer)

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p <- ggplot(top_paths, aes(x = stats::reorder(source, score), y = score, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black") +
  scale_fill_manual(values = c("fibroblast" = colors[2], "fibroblast" = colors[1])) +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways") +
  ylab("Score")

p



uytlibrary(dplyr)

top_paths <- contrast_acts %>%
  arrange(desc(abs(score))) %>%
  slice(1:50)  # top 20 pathways


ggplot(top_paths, aes(x = reorder(source, score), y = score, fill = score)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_gradient2(low = colors[1], mid = "whitesmoke", high = colors[2], midpoint = 0) +
  theme_minimal() +
  theme(axis.title = element_text(face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  xlab("Pathways") +
  ylab("Score")


TF <- read.csv("~/TFpathways.csv", row.names = 1, check.names = FALSE)


TF
mat <- as.matrix(TF)
pheatmap(
  mat,
  scale = "row",        # προαιρετικό
  cluster_cols = FALSE,  # κάνε το FALSE αν δεν θέλεις clustering
  cluster_rows = TRUE
)

mat <- as.matrix(TF)


dim(mat)
rownames(mat)[1:5]
colnames(mat)[1:5]
mat <- t(mat)



tf_sd <- apply(mat, 1, sd)

selected_tfs <- names(sort(tf_sd, decreasing = TRUE))[1:50]

pheatmap(mat[selected_tfs, ], scale = "row",cluster_cols = FALSE)


tf_means <- rowMeans(mat, na.rm = TRUE)

upper <- quantile(tf_means, 0.90)
lower <- quantile(tf_means, 0.10)

selected_tfs <- names(tf_means)[tf_means > upper | tf_means < lower]

selected_tfs



groupA_cols <- 1:8
groupB_cols <- 9:16

pvals <- apply(mat, 1, function(x) {
  wilcox.test(x[groupA_cols], x[groupB_cols])$p.value
})


meanA <- rowMeans(mat[, groupA_cols], na.rm = TRUE)
meanB <- rowMeans(mat[, groupB_cols], na.rm = TRUE)
diff <- meanA - meanB

contrast_selected <- contrast_acts[contrast_acts$source %in% names(pvals)[pvals < 0.05], ]
# Keep only log2FoldChange values
plot_df <- contrast_acts[contrast_acts$measurement == "log2FoldChange", ]

contrast_selected

write.csv(contrast_selected,"cs.csv")
pheatmap(mat[selected, c(groupA_cols, groupB_cols)], scale = "row", cluster_cols = FALSE)

selected <- contrast_selected$source
selected
# Plot
colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p <- ggplot2::ggplot(data = selected, 
                     mapping = ggplot2::aes(x = stats::reorder(source, score), 
                                            y = score)) + 
  ggplot2::geom_bar(mapping = ggplot2::aes(fill = score),
                    color = "black",
                    stat = "identity") +
  ggplot2::scale_fill_gradient2(low = colors[1], 
                                mid = "whitesmoke", 
                                high = colors[2], 
                                midpoint = 0) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.title = element_text(face = "bold", size = 12),
                 axis.text.x = ggplot2::element_text(angle = 45, 
                                                     hjust = 1, 
                                                     size = 10, 
                                                     face = "bold"),
                 axis.text.y = ggplot2::element_text(size = 10, 
                                                     face = "bold"),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank()) +
  ggplot2::xlab("Pathways")

p




# Create a proper data frame
plot_df <- contrast_selected[, c("source", "score")]

# Optional: reorder by score
plot_df$source <- factor(plot_df$source, levels = plot_df$source[order(plot_df$score)])


plot_df$source <- make.unique(as.character(plot_df$source))
plot_df$source <- factor(plot_df$source, levels = plot_df$source[order(plot_df$score)])


library(dplyr)

plot_df <- plot_df %>%
  group_by(source) %>%
  summarise(score = mean(score, na.rm = TRUE)) %>%
  ungroup()

plot_df$source <- factor(plot_df$source, levels = plot_df$source[order(plot_df$score)])



library(ggplot2)
colors <- c("blue", "red")

ggplot(plot_df, aes(x = source, y = score)) +
  geom_bar(aes(fill = score), color = "black", stat = "identity") +
  scale_fill_gradient2(low = colors[1], mid = "whitesmoke", high = colors[2], midpoint = 0) +
  theme_minimal() +
  theme(
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  xlab("Pathways")



colnames(contrast_acts)
plot_df <- contrast_acts[contrast_acts$condition == "log2FoldChange", ]






colnames(df)


pathway <- 'RUNX1'

df <- net %>%
  dplyr::filter(source == pathway) %>%
  dplyr::arrange(target) %>%
  dplyr::mutate(ID = target, 
                color = "3") %>%
  tibble::column_to_rownames('target')

inter <- sort(dplyr::intersect(rownames(deg), rownames(df)))

df <- df[inter, ]

df['t_value'] <- deg[inter, ]

df <- df %>%
  dplyr::mutate(color = dplyr::if_else(weight > 0 & t_value > 0, '1', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight > 0 & t_value < 0, '2', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight < 0 & t_value > 0, '2', color)) %>%
  dplyr::mutate(color = dplyr::if_else(weight < 0 & t_value < 0, '1', color))

colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])

p <- ggplot2::ggplot(data =TF, 
                     mapping = ggplot2::aes(x = weight, 
                                            y = t_value, 
                                            color = color)) + 
  ggplot2::geom_point(size = 2.5, 
                      color = "black") + 
  ggplot2::geom_point(size = 1.5) +
  ggplot2::scale_colour_manual(values = c(colors[2], colors[1], "grey")) +
  ggrepel::geom_label_repel(mapping = ggplot2::aes(label = ID)) + 
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position = "none") +
  ggplot2::geom_vline(xintercept = 0, linetype = 'dotted') +
  ggplot2::geom_hline(yintercept = 0, linetype = 'dotted') +
  ggplot2::ggtitle(pathway)

p
#> Warning: ggrepel: 447 unlabeled data points (too many overlaps). Consider
#> increasing max.overlaps
#> 
#> 
#> 

