##############################################################
# Z-score heatmap — KEGG functional profiles by Treatment
##############################################################

library(tidyverse)
library(pheatmap)
library(stringr)
library(grid)

# ---- Load data ----
combined2 <- readr::read_csv("tables/combined_functional.csv", show_col_types = FALSE)

# ---- Filter treatments + desired KEGG Level-1 categories ----
subset_df <- combined2 %>%
  filter(Treatment %in% c("Basal Diet","Probiotic","BMD","Essential oils"),
         broadgroup %in% c("Metabolism",
                           "Genetic Information Processing")) %>%
  mutate(
    abundance = as.numeric(abundance),
    type      = na_if(type, ""),
    Age_Days  = as.character(Age_Days)
  ) %>%
  drop_na(Treatment, type, abundance)

# ---- Clean function names ----
subset_df <- subset_df %>%
  mutate(
    # remove species suffixes from KEGG names
    type = str_remove(type, " ?- ?(Escherichia coli|yeast|fly|multiple species)$"),
    # normalize weird whitespace/newlines/nbsp and trim
    type = str_replace_all(type, "[\\r\\n]", " "),
    type = str_replace_all(type, "[\\u00A0\\u2007\\u202F]", " "),
    type = str_squish(type)
  ) %>%
  # remove rows that are literally just "N", "n", or "NA" strings
  filter(!str_detect(type, "^\\s*[A-Za-z]\\s*$"),
         type != "NA")

# ---- Identify sample ID column ----
id_col <- intersect(c("Sample_ID","sample_id","file"), names(subset_df))[1]
stopifnot(length(id_col) == 1)

# ---- % within each sample ----
per_sample <- subset_df %>%
  group_by(.data[[id_col]]) %>%
  mutate(pct = 100 * abundance / sum(abundance, na.rm = TRUE)) %>%
  ungroup()

# ---- Average % by Treatment × type ----
avg_by_trt <- per_sample %>%
  group_by(Treatment, type) %>%
  summarise(mean_pct = mean(pct, na.rm = TRUE), .groups = "drop")

# ---- Build matrix (type × Treatment) ----
mat <- avg_by_trt %>%
  pivot_wider(names_from = Treatment, values_from = mean_pct, values_fill = 0) %>%
  filter(!is.na(type), type != "NA", type != "") %>%
  mutate(type = make.unique(type)) %>%
  column_to_rownames("type") %>%
  as.matrix()

# ---- Rename BMD → Antibiotic in columns (if present) ----
colnames(mat) <- trimws(colnames(mat))
colnames(mat) <- sub("^BMD$", "Antibiotic", colnames(mat), ignore.case = FALSE)

# ---- Remove any NA cells / NA rownames just in case ----
mat[is.na(mat)] <- 0
mat <- mat[!(rownames(mat) %in% c("NA","N","")), , drop = FALSE]

# ---- Keep top most-variable rows (robust) ----
if (nrow(mat) > 30) {
  v <- apply(mat, 1, var, na.rm = TRUE)
  n_keep <- min(50, nrow(mat))
  keep <- names(sort(v, decreasing = TRUE))[seq_len(n_keep)]
  mat <- mat[keep, , drop = FALSE]
}

# ---- Z-score by row and cap to [-2, 2] ----
z <- t(scale(t(mat)))
z[is.na(z)] <- 0
z <- pmax(pmin(z, 2), -2)

# ---- Enforce desired treatment column order ----
desired_order <- c("Basal Diet", "Probiotic", "Antibiotic", "Essential oils")
present <- desired_order[desired_order %in% colnames(z)]
z <- z[, present, drop = FALSE]

# ---- Column annotation + custom colors (only for present cols) ----
ann_col <- data.frame(Treatment = present, row.names = present)

treatment_colors <- c(
  "Basal Diet"     = "#996666",
  "Probiotic"      = "#24868EFF",
  "Antibiotic"     = "#C7E020FF",
  "Essential oils" = "#35B779FF"
)
ann_colors <- list(Treatment = treatment_colors[present])

# ---- Plot (no column clustering, fixed order) ----
p <- pheatmap(
  z,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_distance_cols = "euclidean",
  clustering_distance_rows = "euclidean",
  color = colorRampPalette(c("#469D76FF", "#3C4B99FF", "#924099FF"))(200),
  border_color = NA,
  main = "",
  annotation_col = ann_col,
  annotation_colors = ann_colors,
  fontsize_row = 7,
  fontsize_col = 12,
  fontsize = 9,
  show_rownames = TRUE,
  show_colnames = TRUE,
  legend = TRUE
)

