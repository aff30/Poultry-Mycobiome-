library(tidyverse)

func <- read_csv("tables/combined_functional.csv", show_col_types = FALSE)
names(func)
func %>%
  filter(!is.na(ko_id)) %>%
  summarise(total_ASVs_annotated = n_distinct(ko_id))

# ======================================================
# ALDEx2 analysis for KEGG functional profiles
# ======================================================
library(tidyverse)
library(readxl)
library(ALDEx2)

# ---------- 1) Load data ----------
func <- read_csv(file = "tables/combined_functional.csv")
meta <- read_excel("met_Chicken/Functional_metadata.xlsx")

# ---------- 2) Prepare metadata ----------
# Adjust this rename if your metadata column is named differently
meta <- meta %>%
  rename(Sample_ID = Sample_ID) %>%     # change if column name differs
  mutate(Treatment = as.factor(Treatment))

# ---------- 3) Summarize KEGG functional table ----------
# Aggregate abundance per KO per sample
func_sum <- func %>%
  group_by(ko_id, Sample_ID) %>%
  summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")

# Pivot to wide format (features = rows, samples = columns)
func_wide <- func_sum %>%
  pivot_wider(names_from = Sample_ID, values_from = abundance, values_fill = 0)

# Convert to data frame and clean
count_matrix <- as.data.frame(func_wide)
rownames(count_matrix) <- count_matrix$ko_id
count_matrix <- dplyr::select(count_matrix, -ko_id)

# ---------- 4) Match order of samples ----------
meta <- meta %>% filter(Sample_ID %in% colnames(count_matrix))
count_matrix <- count_matrix[, meta$Sample_ID]

# Create condition vector for ALDEx2
conds <- meta$Treatment
names(conds) <- meta$Sample_ID

# ---------- 5) Ensure correct data types ----------
# Convert all count columns to numeric
count_matrix[] <- lapply(count_matrix, function(x) as.numeric(as.character(x)))

# Convert condition vector to character
conds <- as.character(conds)

# Confirm structure
cat("\nMatrix check (first 5 columns):\n")
print(sapply(count_matrix[1:5], class))
cat("\nCondition check:\n")
print(table(conds))
cat("\nAll samples match? ", all(colnames(count_matrix) == names(conds)), "\n")

# ======================================================
# 5.5) Convert to integer-like pseudo-counts for ALDEx2
# ======================================================

# Scale abundances to pseudo-counts (ALDEx2 requires integers)
count_matrix_int <- round(count_matrix * 1e6)
count_matrix_int[count_matrix_int == 0] <- 1  # avoid zeros

# Verify counts are integers
cat("\nValue summary after scaling:\n")
print(summary(as.vector(as.matrix(count_matrix_int))[1:100]))

# ---------- 6) Run ALDEx2 ----------
x <- aldex.clr(count_matrix_int,
               conds = conds,
               mc.samples = 128,     # increase for accuracy (128–256)
               denom = "all",
               verbose = TRUE)

# Kruskal–Wallis test for multi-group differences
kw <- aldex.kw(x)

# ---------- 7) Extract significant KEGG features ----------
sig_kw <- kw %>%
  tibble::rownames_to_column("ko_id") %>%
  filter(kw.eBH < 0.05) %>%          # use the BH-adjusted p-value
  arrange(kw.eBH)

# ---------- 8) Annotate with KEGG metadata ----------
func_ann <- func %>%
  distinct(ko_id, path_id, map_ID, type, broadgroup, subgroup)

sig_kw_annot <- sig_kw %>%
  left_join(func_ann, by = "ko_id")

# ---------- 9) Save results ----------
write.csv(sig_kw_annot, "tables/ALDEx2_KEGG_results.csv", row.names = FALSE)
cat("\n✅ ALDEx2 complete — results saved as ALDEx2_KEGG_results.csv\n")


# ======================================================
# 10) Visualize Top 20 Differential KEGG Pathways
# ======================================================

# Reload annotated results if needed
sig_kw_annot <- read_csv("tables/ALDEx2_KEGG_results.csv", show_col_types = FALSE)

# Rank and select top 20
top20 <- sig_kw_annot %>%
  mutate(effect_size = abs(kw.eBH)) %>%
  arrange(kw.ep) %>%
  slice_head(n = 20) %>%
  mutate(
    label = coalesce(type, path_id, ko_id),
    label = str_trunc(label, 60, "right"),
    label = fct_reorder(label, -log10(kw.ep))
  )

# ---- Plot ----
p <- ggplot(top20, aes(x = -log10(kw.ep), y = label, fill = broadgroup)) +
  geom_col(width = 0.7) +
  scale_fill_brewer(palette = "Set2", name = "KEGG Broadgroup") +
  labs(
    x = expression(-log[10]("p-value")),
    y = "KEGG Pathway"
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )
p

## sugnificant broadgroups ALDEX2
library(tidyverse)
library(forcats)
library(readr)

combined2 <- read_csv("tables/combined_functional.csv")
sig_kw_annot <- read_csv("tables/ALDEx2_KEGG_results.csv", show_col_types = FALSE)

merged_meta <- merge(combined2, sig_kw_annot, by = "ko_id")  # replace "Sample_ID" with your actual common column

# Remove negative and positive controls (assuming column is "Treatment")
meta_clean <- merged_meta %>%
  mutate(
    Treatment = case_when(
      Treatment == "BMD" ~ "Antibiotic",
      TRUE ~ Treatment
    )
  ) %>%
  filter(!Treatment %in% c("N_C", "P_C", "Baseline", "Environment Control"))

merged_meta <- meta_clean

# --- Rank and select top 20 pathways globally (optional: per treatment) ---
top20_all <- merged_meta %>%
  group_by(Treatment) %>%
  arrange(kw.ep, .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup() %>%
  mutate(
    label = coalesce(type.y, path_id.y, ko_id),
    label = str_trunc(label, 60, "right"),
    label = fct_reorder(label, -log10(kw.ep))
  )

# --- Plot per Treatment ---
p <- ggplot(top20_all, aes(x = -log10(kw.ep), y = label, fill = broadgroup.y)) +
  geom_col(width = 0.7) +
  scale_fill_brewer(palette = "Set2", name = "KEGG Broadgroup") +
  labs(
    x = expression(-log[10]("p-value")),
    y = "KEGG Pathway"
  ) +
  facet_wrap(~Treatment, scales = "free_y") +
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom"
  )

library(viridis)

p <- ggplot(top20_all, aes(x = -log10(kw.ep), y = label, fill = broadgroup.y)) +
  geom_col(width = 0.7) +
  scale_fill_viridis_d(option = "D") +  # Options: "A", "B", "C", "D", "E"
  labs(
    x = expression(-log[10]("p-value")),
    y = "KEGG Pathway"
  ) +
  facet_wrap(~Treatment, scales = "free_y") +
  theme_classic() +
  theme(
    strip.text = element_text(face = "bold", size = 12),
    axis.text.y = element_text(size = 8),
    legend.position = "bottom"
  )


p
