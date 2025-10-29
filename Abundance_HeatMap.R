# Load required packages
library(tidyverse)
library(reshape2)
library(pheatmap)   # for heatmaps
library(vegan)      # for ordination like PCA/NMDS
library(ggplot2)    # for custom plots

combined2 <- read_csv(file = "tables/combined_functional.csv")

subset_df <- combined2 %>% filter(Treatment == "Basal Diet" | Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Essential oils")

# Group by Sample and KEGG pathway
library(tidyverse)

# Group by Sample and KEGG pathway
kegg_plot_data <- subset_df %>%
  group_by(Treatment, subgroup) %>%
  summarize(total_abundance = sum(abundance)) %>%
  ungroup()

# Plot top 10 pathways
top_paths <- kegg_plot_data %>%
  group_by(subgroup) %>%
  summarize(mean_abundance = mean(total_abundance)) %>%
  top_n(15, mean_abundance) %>%
  pull(subgroup)

ggplot(kegg_plot_data %>% filter(subgroup %in% top_paths),
       aes(x = Treatment, y = total_abundance, fill = subgroup)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Sample", y = "Pathway Abundance", fill = "KEGG Pathway") +
  scale_color_manual() +
  theme_classic() 


library(reshape2)
library(pheatmap)

# Summarize by Sample and KO
ko_matrix <- subset_df %>%
  group_by(Treatment, ko_id) %>%
  summarize(abundance = sum(abundance)) %>%
  pivot_wider(names_from = ko_id, values_from = abundance, values_fill = 0) %>%
  column_to_rownames("Treatment")

# Draw heatmap
pheatmap(ko_matrix, scale = "row", clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean")


#heatmap
library(tidyverse)

# Keep only your four treatments; ensure numeric abundance and clean keys
subset_df <- combined2 %>%
  filter(Treatment %in% c("Basal Diet","Probiotic","BMD","Essential oils")) %>%
  mutate(
    abundance = as.numeric(abundance),
    type = na_if(type, ""),          # turn empty strings into NA
    Age_Days = as.character(Age_Days)
  ) %>%
  drop_na(Treatment, Age_Days, type)  # drop rows missing keys

# Summarize, complete the full grid, and fill missing totals with 0
kegg_plot_data <- subset_df %>%
  group_by(Phase, Treatment, Age_Days, type) %>%
  summarise(total_abundance = sum(abundance, na.rm = TRUE), .groups = "drop") %>%
  # Ensure only the ages you want, in order
  mutate(Age_Days = factor(Age_Days, levels = c("1","10","21"))) %>%
  # Complete the grid so every Treatment x Age_Days x type exists
  tidyr::complete(Treatment, Age_Days, type, fill = list(total_abundance = 0L))

# Choose top 25 types (exclude any accidental NA)
top_types <- kegg_plot_data %>%
  filter(!is.na(type)) %>%
  group_by(type) %>%
  summarise(mean_abundance = mean(total_abundance), .groups = "drop") %>%
  slice_max(mean_abundance, n = 25) %>%
  pull(type)

# Plot (no NAs left; but set a safety na.value anyway)
ggplot(
  kegg_plot_data %>% filter(type %in% top_types),
  aes(x = Age_Days, y = type, fill = total_abundance)
) +
  geom_tile() +
  facet_grid(~ Treatment, drop = FALSE) +
  scale_fill_gradient(name = "Total Abundance",
                      low = "#FFFFFF", high = "#D95E31FF",
                      na.value = "#FFFFFF") +
  labs(x = "Age (Days)", y = NULL, title = "Bacterial Functional Type Abundance") +
  theme_bw() +
  theme(
    strip.placement = "outside",
    strip.background = element_rect(fill = "#EEEEEE", color = "#FFFFFF"),
    plot.title = element_text(hjust = 0.5)
  )

## z-score heatmap
library(tidyverse)
library(pheatmap)

df0 <- combined2 %>%
  filter(Treatment %in% c("Basal Diet","Probiotic","BMD","Essential oils")) %>%
  mutate(
    abundance = as.numeric(abundance),
    type = na_if(type, "")
  ) %>%
  drop_na(Treatment, type, abundance)

# ---- choose a sample ID column (robust) ----
id_col <- intersect(c("Sample_ID","sample_id","file"), names(df0))[1]
stopifnot(length(id_col) == 1)

# ---- 1) % within each sample ----
per_sample <- df0 %>%
  group_by(.data[[id_col]]) %>%
  mutate(pct = 100 * abundance / sum(abundance, na.rm = TRUE)) %>%
  ungroup()

# ---- 2) average % by Treatment × type ----
avg_by_trt <- per_sample %>%
  group_by(Treatment, type) %>%
  summarise(mean_pct = mean(pct, na.rm = TRUE), .groups = "drop")

# ---- 3) make matrix type × Treatment ----
mat <- avg_by_trt %>%
  pivot_wider(names_from = Treatment, values_from = mean_pct, values_fill = 0) %>%
  filter(!is.na(type)) %>%
  mutate(type = make.unique(type)) %>%
  column_to_rownames("type") %>%
  as.matrix()

# ---- 4) keep informative rows (top 50 by variance across treatments) ----
if (nrow(mat) > 35) {
  v <- apply(mat, 1, var, na.rm = TRUE)
  keep <- names(sort(v, decreasing = TRUE))[1:50]
  mat <- mat[keep, , drop = FALSE]
}

# ---- 5) row-wise Z-score and cap to [-2, 2] to avoid blown-out colors ----
z <- t(scale(t(mat)))                       # Z-score per pathway
z[is.nan(z)] <- 0                           # rows with 0 variance -> 0
z <- pmax(pmin(z, 2), -2)                   # cap extremes

# ---- 6) column annotation (optional) ----
ann_col <- data.frame(Treatment = colnames(z))
rownames(ann_col) <- colnames(z)

# ---- 7) draw heatmap ----
pheatmap(
  z,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  color = colorRampPalette(c("#343854FF","#8C384DFF","#CF2438FF"))(200),
  border_color = NA,
  main = "KEGG pathways — mean % per Treatment",
  annotation_col = ann_col,
  fontsize_row = 7,
  fontsize_col = 10
)


