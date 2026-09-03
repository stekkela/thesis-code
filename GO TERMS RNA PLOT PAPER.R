library(readr)
library(dplyr)
library(ggplot2)

go_df <- read_csv("RNA_GO_file.csv",show_col_types = FALSE)

go_df
go_df <- go_df %>%
  mutate(
    Category = case_when(
      Term == "GOTERM_BP_DIRECT" ~ "Biological Process",
      Term == "GOTERM_CC_DIRECT" ~ "Cellular Component",
      Term == "GOTERM_MF_DIRECT" ~ "Molecular Function"
    )
  )

go_df <- go_df %>%
  group_by(Category) %>%
  arrange(Score) %>%
  mutate(`Term name` = factor(`Term name`, levels = `Term name`)) %>%
  ungroup()

ggplot(go_df, aes(x = Score, y = `Term name`, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(
    "Biological Process" = "#5AAE61",
    "Cellular Component" = "#E08214",
    "Molecular Function" = "#4393C3"
  )) +
  labs(
    x = expression(-log[10]~p~value),
    y = NULL
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.text.y = element_text(size = 9)
  )+ labs(title = "GO Enrichment of DEGs in Keloids") +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 14
    )
  )

