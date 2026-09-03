library(ggplot2)
library(dplyr)
library(ggrepel)

# Calculate 95th percentile of KDF_NDF_change (binding score)
score_cutoff <- quantile(df_plot$KDF_NDF_change, 0.95, na.rm = TRUE)

# Add significance column
df_plot <- df_plot %>%
  mutate(
    sig = case_when(
      KDF_NDF_change > score_cutoff ~ "Up",
      KDF_NDF_change < -score_cutoff ~ "Down",
      TRUE ~ "NS"
    )
  )

# Map "Up" → "Keloid", "Down" → "Normal" for legend labels
legend_labels <- c("Up" = "Keloid", "Down" = "Normal", "NS" = "NS")

# Define colors
cols <- c("Up" = "red", "Down" = "blue", "NS" = "grey")

# Create plot with labels for significant points
vol <- ggplot(df_plot, aes(x = KDF_NDF_change, y = -log10(padj))) +
  geom_point(aes(color = sig), alpha = 0.75, size = 2.2) +
  scale_color_manual(
    values = cols,
    labels = legend_labels
  ) +
  geom_vline(xintercept = c(-score_cutoff, score_cutoff), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(
    data = df_plot %>% filter(sig != "NS"),
    aes(label = Gene), # Replace 'Gene' with the name column in your df
    size = 3,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = 20
  ) +
  labs(
    title = "Volcano-like Plot: Binding Score vs p-value",
    subtitle = paste0("Cutoffs: |KDF_NDF_change| > 95th percentile"),
    x = "Binding Score (KDF_NDF_change)",
    y = "-log10(adj. p-value)",
    color = "Expression"
  ) +
  theme_minimal()

# Save plot
ggsave("binding_score_vs_pvalue_labeled.png", plot = vol, width = 8, height = 6, dpi = 300)
