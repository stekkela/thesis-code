

library(dplyr)
library(ggplot2)

# Read your file
DF <- read.csv("~/cellchar_2025-Oct-09-17-09-18_Single_Target_Data.csv", check.names = FALSE)
# 2️⃣ Automatically extract the row letter (A–H) from the "Well Label" column
DF <- DF %>%
  mutate(Row = str_extract(`WELL LABEL`, "^[A-H]"))

# 3️⃣ Assign cell types based on row letter
DF <- DF %>%
  mutate(Row = str_extract(`WELL LABEL`, "^[A-H]")) %>% 
  mutate(CellType = case_when(
    Row == "A" ~ "NDFC",
    Row == "B" ~ "KDFC",
    Row == "C" ~ "NDFT",
    Row == "D" ~ "KDFT",
    Row == "E" ~ "NDFC",
    Row == "F" ~ "KDFC",
    Row == "G" ~ "NDFT",
    Row == "H" ~ "KDFT",
    TRUE ~ "Unknown"
  ))
names(DF) <- trimws(names(DF), which = "both")
names(DF) <- make.names(names(DF))
well_means <- DF %>%
  group_by(`WELL LABEL`, CellType) %>%
  summarise(mean_cell_area = mean(`Nuclei Area wv1`, na.rm = TRUE), .groups = "drop")
            
well_means <- DF %>%
  group_by(WELL.LABEL, CellType) %>%
  summarise(mean_cell_area = mean(Nuclei.Area.wv1, na.rm = TRUE), .groups = "drop")
celltype_stats <- well_means %>%
  group_by(CellType) %>%
  summarise(
    mean_area = mean(mean_cell_area, na.rm = TRUE),
    sd_area = sd(mean_cell_area, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
celltype_stats <- well_means %>%
  group_by(CellType) %>%
  summarise(
    mean_area = mean(mean_cell_area, na.rm = TRUE),
    sd_area = sd(mean_cell_area, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )


DF <- read.csv("~/cellchar_2025-Oct-09-17-09-18_Single_Target_Data.csv", check.names = FALSE)

# ==========================
# 3️⃣ Clean column names
# ==========================
# Remove leading/trailing spaces and make names syntactically valid
names(DF) <- make.names(trimws(names(DF)))

# ==========================
# 4️⃣ Create CellType column based on Row letter
# ==========================
DF <- DF %>%
  mutate(Row = str_extract(WELL.LABEL, "^[A-H]")) %>%
  mutate(CellType = case_when(
    Row %in% c("A","E") ~ "NDFC",
    Row %in% c("B","F") ~ "KDFC",
    Row %in% c("C","G") ~ "NDFT",
    Row %in% c("D","H") ~ "KDFT",
    TRUE ~ "Unknown"
  ))

# ==========================
# 5️⃣ Compute mean Nuclei.Area.wv1 per well
# ==========================
well_means <- DF %>%
  group_by(WELL.LABEL, CellType) %>%
  summarise(mean_cell_area = mean(Nuclei.Area.wv1, na.rm = TRUE), .groups = "drop")


celltype_stats <- well_means %>%
  group_by(CellType) %>%
  summarise(
    mean_area = mean(mean_cell_area, na.rm = TRUE),
    sd_area = sd(mean_cell_area, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
dput(colnames(DF))

























