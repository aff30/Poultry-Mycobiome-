#graph article relative abundance 

library(ggplot2)
library(readxl)
library(dplyr)
library(tidyverse)
library(ggtext)

df <- readxl::read_excel("tables/R_graph_article_maingenera_day.xlsx")


# Convert to long format
df_long <- df %>%
  pivot_longer(cols = -c(Treatment, Days),
               names_to = "Genus",
               values_to = "Abundance")

# Convert Days to numeric (if it's not already)
df_long$Days <- as.numeric(df_long$Days)

# Assign manual colors (HEX codes from your image, approx top to bottom)
genus_colors <- c(
  "Candida" = "#A6CEE3",
  "Fusarium" = "#1F78B4",
  "Aspergillus" = "#B2DF8A",
  "Malassezia" = "#FB9A99",
  "Pyricularia" = "#E31A1C",
  "Saccharomyces" = "#99fff7",
  "Penicillium" = "#432e51"
)

# Plot bar chart per day
ggplot(df_long, aes(x = Treatment, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7) +
  facet_wrap(~Days, ncol = 3) +  # one plot per day
  scale_fill_manual(values = genus_colors) +  # optional: your custom color palette
  labs(title = "Fungal Genus Composition per Day",
       x = "Treatment",
       y = "Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold"),
        legend.position = "right")

# Plot: Line plots faceted by Day
ggplot(df_long, aes(x = Treatment, y = Abundance, group = Genus, color = Genus)) +
  geom_line(aes(group = Genus), size = 0.5) +
  geom_point(size = 2) +
  facet_wrap(~Days, ncol = 3) +
  scale_color_manual(values = genus_colors) +
  labs(title = "",
       x = "Treatment",
       y = "Relative Abundance (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right",
        strip.text = element_text(face = "bold"))

