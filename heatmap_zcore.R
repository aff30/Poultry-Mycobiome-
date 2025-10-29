##############################################################
# Functional KEGG analysis — Fungal Metabolic & Core Processes
# Author: Ana Fonseca
# Date: 2025-10-27
##############################################################

# ---- Load packages ----
library(tidyverse)
library(pheatmap)
library(stringr)

# ---- Load data ----
combined2 <- read_csv("tables/combined_functional.csv")

# ---- Filter relevant treatments and KEGG Level 1 categories ----
subset_df <- combined2 %>%
  filter(Treatment %in% c("Basal Diet", "Probiotic", "BMD", "Essential oils")) %>%
  filter(broadgroup %in% c(
    "Metabolism",
    "Genetic Information Processing"
  )) %>%
  mutate(
    abundance = as.numeric(abundance),
    type = na_if(type, ""),
    Age_Days = as.character(Age_Days)
  ) %>%
  drop_na(Treatment, type, abundance)

# ---- Clean pathway names (remove species suffixes) ----
subset_df <- subset_df %>%
  mutate(
    type = str_remove(type, " ?- ?(Escherichia coli|yeast|fly|multiple species|HIV-1)$")
  )


##############################################################
# Z-score heatmap — KEGG functional profiles by Treatment
##############################################################

# Identify sample ID column
id_col <- intersect(c("Sample_ID", "sample_id", "file"), names(subset_df))[1]
stopifnot(length(id_col) == 1)

# ---- 1) Convert to % within each sample ----
per_sample <- subset_df %>%
  group_by(.data[[id_col]]) %>%
  mutate(pct = 100 * abundance / sum(abundance, na.rm = TRUE)) %>%
  ungroup()

# ---- 2) Average % by Treatment × type ----
avg_by_trt <- per_sample %>%
  group_by(Treatment, type) %>%
  summarise(mean_pct = mean(pct, na.rm = TRUE), .groups = "drop")

# ---- 3) Build matrix (type × Treatment) ----
mat <- avg_by_trt %>%
  pivot_wider(names_from = Treatment, values_from = mean_pct, values_fill = 0) %>%
  filter(!is.na(type)) %>%
  mutate(type = make.unique(type)) %>%
  column_to_rownames("type") %>%
  as.matrix()

# ---- 4) Keep top 50 most variable pathways ----
if (nrow(mat) > 30) {
  v <- apply(mat, 1, var, na.rm = TRUE)
  keep <- names(sort(v, decreasing = TRUE))[1:50]
  mat <- mat[keep, , drop = FALSE]
}

# ---- 5) Z-score per row, cap at [-2, 2] ----
z <- t(scale(t(mat)))
z[is.nan(z)] <- 0
z <- pmax(pmin(z, 2), -2)

# ---- 6) Rename “BMD” → “Antibiotic” ----
colnames(z) <- gsub("^BMD$", "Antibiotic", colnames(z))

# ---- 7) Column annotation + custom colors ----
ann_col <- data.frame(Treatment = colnames(z))
rownames(ann_col) <- colnames(z)

treatment_colors <- c(
  "Basal Diet"     = "#996666",  
  "Probiotic"      = "#24868EFF",  
  "Antibiotic"     = "#C7E020FF", 
  "Essential oils" = "#35B779FF"   
)
ann_colors <- list(Treatment = treatment_colors)

# --- normalize column names first ---
colnames(z) <- trimws(colnames(z))
colnames(z) <- sub("^BMD$", "Antibiotic", colnames(z), ignore.case = FALSE)

# ---- 4) Keep top most variable pathways (robust) ----
if (nrow(mat) > 30) {
  v <- apply(mat, 1, var, na.rm = TRUE)
  n_keep <- min(50, nrow(mat))
  keep <- names(sort(v, decreasing = TRUE))[seq_len(n_keep)]
  mat <- mat[keep, , drop = FALSE]
}

# ---- 5) Z-score per row, cap at [-2, 2] ----
z <- t(scale(t(mat)))
z[is.nan(z)] <- 0
z <- pmax(pmin(z, 2), -2)

# ---- 6) Rename “BMD” → “Antibiotic” and normalize names ----
colnames(z) <- trimws(colnames(z))
colnames(z) <- sub("^BMD$", "Antibiotic", colnames(z), ignore.case = FALSE)

# ---- 7) Column annotation + custom colors ----
desired_order <- c("Basal Diet", "Probiotic", "Antibiotic", "Essential oils")
present <- desired_order[desired_order %in% colnames(z)]

# enforce order
z <- z[, present, drop = FALSE]

ann_col <- data.frame(Treatment = present, row.names = present)

treatment_colors <- c(
  "Basal Diet"     = "#996666",
  "Probiotic"      = "#24868EFF",
  "Antibiotic"     = "#C7E020FF",
  "Essential oils" = "#35B779FF"
)
# keep only colors for present columns
ann_colors <- list(Treatment = treatment_colors[present])

# ---- 8) Plot (turn OFF column clustering!) ----
pheatmap(
  z,
  cluster_cols = FALSE,                 # <--- keep your specified order
  cluster_rows = TRUE,
  clustering_distance_cols = "euclidean",
  clustering_distance_rows = "euclidean",
  color = colorRampPalette(c("#A3BAC2FF", "#449DB3FF", "#8C6E5DFF"))(200),
  border_color = NA,
  main = "",
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  fontsize_row = 7,
  fontsize_col = 12,
  fontsize = 9,
  show_rownames = TRUE,
  show_colnames = TRUE
)

