# Subset face adipocytes
dds_fibroblast <- dds[dds$media == "fibroblast"]

# Drop unused levels
dds_fibroblast$body_site <- droplevels(dds_fibroblast$body_site)
dds_fibroblast$media     <- droplevels(dds_fibroblast$media)
dds_fibroblast$cell_type <- droplevels(dds_fibroblast$cell_type)
# Make sure 'normal' is the reference
dds_fibroblast$cell_type <- relevel(dds_fibroblast$cell_type, "normal")
design(dds_fibroblast) <- ~ cell_type
dds_fibroblast <- DESeq(dds_fibroblast)
res_fibroblast <- results(dds_fibroblast, contrast = c("cell_type", "keloid", "normal"))
dds_fibroblast
levels(dds_fibroblast$cell_type)

library(EnhancedVolcano)
EnhancedVolcano(res_fibroblast,
                lab = rownames(res_fibroblast),   # gene names
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                title = 'KDF vs NDF fibroblast media face',
                subtitle = 'DESeq2 Results',
                legendPosition = 'right',xlim = c(-5, 5),   # set same x-axis limits
                ylim = c(0, 10))

# Convert to a regular data frame
res_adipocyte_abdomen_df <- as.data.frame(res_fibroblast)

write.csv(res_fibroblast,
          file = "C:/Users/hmz255/OneDrive - Queen Mary, University of London/Documents/RNA-seqAdipo/res_fibroblast_NDFVSKDF.csv",
          row.names = TRUE)



