# Load required libraries
library(phyloseq)
library(tidyverse)
library(pheatmap)
library(grid)


ps_fungi <- readRDS("ps/ps_taxfungi.RDS")
ps_fungi

ps_fungi <- subset_samples(ps_fungi, Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Basal Diet" | Treatment == "Essential oils")

ps_bacteria <- readRDS("ps/ps_taxbac.RDS")
ps_bacteria

ps_bacteria <- subset_samples(ps_bacteria, Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Basal Diet" | Treatment == "Essential oils")

# 1. Transform to relative abundance
ps_fungi_rel <- transform_sample_counts(ps_fungi, function(x) x / sum(x))
ps_bacteria_rel <- transform_sample_counts(ps_bacteria, function(x) x / sum(x))

# 2. Aggregate to Genus level
tax_glom_safe <- function(ps, rank) {
  tryCatch({ tax_glom(ps, taxrank = rank) }, error = function(e) ps)
}
ps_fungi_genus <- tax_glom_safe(ps_fungi_rel, "Genus")
ps_bacteria_genus <- tax_glom_safe(ps_bacteria_rel, "Genus")

# 3. Extract taxonomic labels
fungi_labels <- tax_table(ps_fungi_genus)[, "Genus"]
bacteria_labels <- tax_table(ps_bacteria_genus)[, "Genus"]

# 4. Build abundance matrices and transpose
fungi_df <- as.data.frame(t(otu_table(ps_fungi_genus)))
bacteria_df <- as.data.frame(t(otu_table(ps_bacteria_genus)))

# 5. Select top 10 abundant taxa
fungi_top <- names(sort(colMeans(fungi_df), decreasing = TRUE))[1:15]
bacteria_top <- names(sort(colMeans(bacteria_df), decreasing = TRUE))[1:15]
fungi_df <- fungi_df[, fungi_top]
bacteria_df <- bacteria_df[, bacteria_top]

# 6. Match shared samples
shared_samples <- intersect(rownames(fungi_df), rownames(bacteria_df))
fungi_df <- fungi_df[shared_samples, ]
bacteria_df <- bacteria_df[shared_samples, ]

# 7. Load metadata
meta <- sample_data(ps_fungi) %>% as.data.frame()
meta <- meta[shared_samples, ]

# 8. Loop through Treatment x Age_Days groups
meta$Group <- paste(meta$Treatment, meta$Age_Days, sep = "_Day")
groups <- unique(meta$Group)

for (grp in groups) {
  grp_samples <- rownames(meta[meta$Group == grp, ])
  if (length(grp_samples) >= 5) {
    fungi_sub <- fungi_df[grp_samples, ]
    bacteria_sub <- bacteria_df[grp_samples, ]
    
    # Remove taxa with zero standard deviation
    fungi_sub <- fungi_sub[, apply(fungi_sub, 2, sd) != 0, drop = FALSE]
    bacteria_sub <- bacteria_sub[, apply(bacteria_sub, 2, sd) != 0, drop = FALSE]
    
    if (ncol(fungi_sub) > 0 && ncol(bacteria_sub) > 0) {
      combined_df <- cbind(fungi_sub, bacteria_sub)
      cor_matrix <- cor(combined_df, method = "spearman")
      fungi_bac_cor <- cor_matrix[colnames(fungi_sub), colnames(bacteria_sub)]
      
      # Clean and italicize labels
      clean_and_italic <- function(x) {
        cleaned <- gsub("^g__", "", x)           # Remove 'g__' prefix
        formatted <- paste0("italic('", cleaned, "')")  # Make parse-friendly
        return(formatted)
      }
      
      rownames(fungi_bac_cor) <- clean_and_italic(as.character(fungi_labels[colnames(fungi_sub)]))
      colnames(fungi_bac_cor) <- clean_and_italic(as.character(bacteria_labels[colnames(bacteria_sub)]))
      
# Title formatting
treatment_part <- sub("_Day.*", "", grp)
      
# Use dplyr's version explicitly
treatment_part <- dplyr::recode(treatment_part,
                  "BMD" = "Antibiotic",
                  "Probiotic" = "Probiotic",
                  "Basal Diet" = "Basal Diet",
                  "Essential oils" = "Essential oils")
      
      day_part <- sub(".*_Day", "d", grp)
      main_title <- paste(treatment_part, "-", day_part)
      
      #Plot 
      pheatmap(
        fungi_bac_cor,
        main = "",
        cluster_rows = TRUE,
        cluster_cols = TRUE,
        labels_row = parse(text = rownames(fungi_bac_cor)),   # ← italic
        labels_col = parse(text = colnames(fungi_bac_cor)),
        fontsize = 14,
        fontsize_row = 16,  # increase row label size (fungi genus)
        fontsize_col = 16   # increase column label size (bacteria genus)
      )
      
      # Add left-aligned title manually
      grid.text(
        label = main_title,
        x = unit(0.98, "npc"),  # left-side (normalized parent coordinates)
        y = unit(0.98, "npc"),  # top-side
        just = c("right", "top"),
        gp = gpar(fontsize = 12, fontface = "bold")
      )
    }
  }
}

# Initialize results list
all_results <- list()

for (grp in groups) {
  grp_samples <- rownames(meta[meta$Group == grp, ])
  if (length(grp_samples) >= 5) {
    fungi_sub <- fungi_df[grp_samples, ]
    bacteria_sub <- bacteria_df[grp_samples, ]
    
    # Remove zero-variance taxa
    fungi_sub <- fungi_sub[, apply(fungi_sub, 2, sd) != 0, drop = FALSE]
    bacteria_sub <- bacteria_sub[, apply(bacteria_sub, 2, sd) != 0, drop = FALSE]
    
    if (ncol(fungi_sub) > 0 && ncol(bacteria_sub) > 0) {
      cor_table <- data.frame()
      
      for (f_genus in colnames(fungi_sub)) {
        for (b_genus in colnames(bacteria_sub)) {
          # Spearman correlation
          cor_test <- cor.test(fungi_sub[[f_genus]], bacteria_sub[[b_genus]], method = "spearman")
          
          # Get genus names (cleaned, but not italicized for table)
          f_name <- as.character(fungi_labels[f_genus])
          b_name <- as.character(bacteria_labels[b_genus])
          
          cor_table <- rbind(cor_table, data.frame(
            Fungal_Genus = ifelse(is.na(f_name), f_genus, f_name),
            Bacterial_Genus = ifelse(is.na(b_name), b_genus, b_name),
            Group = grp,
            Spearman_Rho = cor_test$estimate,
            P_value = cor_test$p.value,
            Mean_Fungal_Abundance = mean(fungi_sub[[f_genus]]),
            Mean_Bacterial_Abundance = mean(bacteria_sub[[b_genus]])
          ))
        }
      }
      
      all_results[[grp]] <- cor_table
    }
  }
}

# Combine all results
final_table <- do.call(rbind, all_results)

# View top results
head(final_table)

# Optional: Filter strong correlations
strong_results <- final_table %>% filter(abs(Spearman_Rho) > 0.5 & P_value < 0.05)

# Save to CSV
write.csv(final_table, "tables/fungi_bacteria_spearman_results.csv", row.names = FALSE)


##Get the significant heatmap
library(dplyr)
library(reshape2)
library(pheatmap)

# ---- 1. Filter for significant results
sig_table <- final_table %>%
  filter(P_value < 0.05)

# ---- 2. Clean genus names (remove 'g__' prefix)
clean_taxa <- function(x) gsub("^g__", "", x)
sig_table <- sig_table %>%
  mutate(
    Fungal_Genus = clean_taxa(Fungal_Genus),
    Bacterial_Genus = clean_taxa(Bacterial_Genus),
    Pair = paste0(Fungal_Genus, " ~ ", Bacterial_Genus)
  )

# ---- 3. Reshape to matrix: rows = pairs, cols = Group (Treatment_Day), values = Spearman_Rho
heatmap_matrix <- dcast(sig_table, Pair ~ Group, value.var = "Spearman_Rho")
rownames(heatmap_matrix) <- heatmap_matrix$Pair
heatmap_matrix$Pair <- NULL

# ---- 4. Replace NA with 0 (optional: or keep as NA to leave blank in heatmap)
heatmap_matrix[is.na(heatmap_matrix)] <- 0

# ---- 5. Order columns by Day (1 → 10 → 21), ignoring treatment
# Extract Group names (e.g., "Probiotic_d1", "BMD_d10", ...)
column_names <- colnames(heatmap_matrix)

# Create data.frame with extracted Day and Treatment
column_info <- data.frame(
  name = column_names,
  Day = as.numeric(gsub(".*_d", "", column_names)),
  Treatment = gsub("_d.*", "", column_names),
  stringsAsFactors = FALSE
)

# Arrange by Day first, then Treatment (if multiple per day)
column_ordered <- column_info %>%
  arrange(Day, Treatment) %>%
  pull(name)

# Reorder columns in matrix
heatmap_matrix <- heatmap_matrix[, column_ordered]

# ---- 6. Create gaps between Day groups for visual separation in the heatmap
# Find positions where Day changes
day_vector <- column_info$Day[match(column_ordered, column_info$name)]
gaps_col <- which(diff(day_vector) != 0)

# ---- 7. Plot the heatmap
pheatmap(
  mat = heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize = 10,
  fontsize_row = 5,    # << smaller font for pair names
  fontsize_col = 10,
  gaps_col = gaps_col,
  main = "Significant Fungi–Bacteria Correlations (Grouped by Day)",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101)
)

