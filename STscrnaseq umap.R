library(RColorBrewer)

n_idents   <- length(unique(combined.all_2$orig.ident))
n_clusters <- length(levels(combined.all_2$seurat_clusters))

cols_idents <- colorRampPalette(brewer.pal(12, "Set3"))(n_idents)
cols_clusters <- colorRampPalette(brewer.pal(12, "Set3"))(n_clusters)

p1 <- DimPlot(
  combined.all_2,
  reduction = "umap",
  group.by = "orig.ident",
  cols = cols_idents
)

p2 <- DimPlot(
  combined.all_2,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  cols = cols_clusters
)

p1 + p2

nature_muted <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999",
  "#8DA0CB", "#FC8D62", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3", "#66C2A5", "#FC8D62"
)

p2 <- DimPlot(
  combined.all_2,
  reduction = "umap",
  label = FALSE,
  repel = TRUE,
  cols = nature_muted[1:n_clusters],raster=FALSE
)

p2




FeaturePlot(
  combined.all_2,
  features = c("COL1A1", "COL1A2", "DCN", "LUM", "COL3A1"),
  min.cutoff = "q10",
  max.cutoff = "q90",
  cols = c("#F5E6EB", "#7A1E48")
)



FeaturePlot(
  combined.all_2,
  features = c("COL1A1", "COL1A2", "VIM", "LUM", "COL3A1"),
  min.cutoff = "q10",
  max.cutoff = "q90",
  cols = c("#FAEEF2", "#9A3B63")
)


FeaturePlot(
  combined.all_2,
  features = c("COL1A1", "COL1A2", "DCN", "LUM", "COL3A1","VIM"),
  min.cutoff = "q10",
  max.cutoff = "q90",
  cols = c("# ", "#7A1E48")
)

library(Seurat)
library(ggplot2)

library(Seurat)
library(ggplot2)

# Define 11 custom colors
custom_colors <- c(
  "#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e",
  "#e6ab02","#a6761d","#666666","#e41a1c","#377eb8","#984ea3"
)

# UMAP with cluster labels and custom colors
p2 <- DimPlot(
  combined.all_2,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  cols = custom_colors,
  raster = FALSE
)

p2




library(Seurat)
library(ggplot2)

custom_colors <- c(
  "#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e",
  "#e6ab02","#a6761d","#666666","#e41a1c","#377eb8",
  "#984ea3","#ff7f00","#ffff33","#a65628","#f781bf",
  "#999999","#8dd3c7","#bebada","#fb8072","#80b1d3",
  "#b3de69","#fccde5"
)


# UMAP with orig.ident, labels, and custom colors
p1 <- DimPlot(
  combined.all_2,
  reduction = "umap",
  group.by = "orig.ident",
  label = TRUE,
  repel = TRUE,
  cols = custom_colors,
  raster = FALSE
)

p1




cluster_colors <- c(
  "#1b9e77","#d95f02","#7570b3","#e7298a",
  "#66a61e","#e6ab02","#a6761d","#666666"
)

# UMAP colored by clusters with your custom colors
DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "seurat_clusters",
  cols = cluster_colors,
  label = TRUE,       # optional: show cluster labels
  repel = TRUE,       # avoid overlap
  raster = FALSE
)


# Create a temporary folder outside OneDrive
dir.create("C:/R_temp", showWarnings  = FALSE)

# Write the CSV there
write.csv(direder.fb.markers,
          "C:/R_temp/all.kfb.markers_2.csv",
          row.names = FALSE)



# Get counts
counts <- GetAssayData(direder_fb_only, slot = "counts")  # genes x cells

# Get samples that actually exist
samples_present <- intersect(unique(direder_fb_only$orig.ident), samples)



pseudo_bulk <- matrix(0, nrow = nrow(counts), ncol = length(samples_present))
rownames(pseudo_bulk) <- rownames(counts)
colnames(pseudo_bulk) <- samples_present

for (s in samples_present) {
  # Get column indices of cells for this sample
  cell_idx <- which(direder_fb_only$orig.ident == s)
  
  if (length(cell_idx) > 0) {
    # Use column indices to subset counts safely
    pseudo_bulk[, s] <- rowSums(counts[, cell_idx, drop = FALSE])
  }
}









