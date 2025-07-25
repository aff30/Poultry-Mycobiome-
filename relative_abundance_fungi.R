#relative abundance

#load packages 
#generally useful ones 
library(BiocManager)
library(stringr)
library(tidyverse)
library(tidylog)

#ones we specifically need 
library(phyloseq)
library(ggplot2)
library(vegan) # version 2.6-4
library(decontam)
library(microViz)

#Load data
ps <- readRDS("ps/ps_taxfungi.RDS")
ps

ps <- subset_samples(ps, Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Basal Diet" | Treatment == "Essential oils")

ps1 <- subset_samples(ps, Treatment == "Probiotic")
ps2 <- subset_samples(ps, Treatment == "BMD")
ps3 <- subset_samples(ps, Treatment == "Basal Diet")
ps4 <- subset_samples(ps, Treatment == "Essential oils")
ps5 <- subset_samples(ps, Treatment == "Baseline")

hueRank <- "Genus"
hueRankPlural <- "Genera"
shadeRank <- "Species"

# Sort phyloseq at lower, and then higher ranks
pseq2 <- ps5 %>%
  tax_fix(unknowns = c("s__kudriavzevii", "s__oryzae", "s__phaffii")) %>% 
  tax_sort(by = sum, at = shadeRank) %>%
  tax_sort(by = sum, at = hueRank) %>%
  tax_agg(rank = shadeRank)

# Specify number of hues and shades desired
nHues <- 3 # "Other" phyla will be shades of grey
nShades <- 4 # "Other" families will be the lightest shade of each hue

hierarchicalPalInfo <- data.frame(
  hue = as.vector(tt_get(pseq2)[, hueRank]),
  shade = as.vector(tt_get(pseq2)[, shadeRank]),
  counts = taxa_sums(otu_get(pseq2))
)

hierarchicalPalInfo <- hierarchicalPalInfo %>%
  dplyr::mutate(
    hue = forcats::fct_other(
      f = hue, keep = unique(hue)[seq_len(nHues)],
      other_level = paste("Other", hueRankPlural)
    ),
    nChrHue = nchar(as.character(hue)), padHue = max(nChrHue) - nChrHue
  ) %>%
  dplyr::group_by(hue) %>%
  dplyr::mutate(
    shade = forcats::fct_other(
      f = shade, keep = unique(shade)[seq_len(nShades - 1)],
      other_level = "Other"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    nChrShade = nchar(as.character(shade)), padShade = max(nChrShade) - nChrShade,
    Taxa = paste0(hue, ": ", strrep(" ", padHue), shade, strrep(" ", padShade))
  )

hierarchicalPalMatrix <- matrix(
  data = sapply(
    X = seq(from = 30, to = 75, length.out = nShades),
    FUN = function(l) scales::hue_pal(l = l, h.start = 120)(n = nHues)
  ),
  byrow = TRUE, ncol = nHues
)
hierarchicalPalMatrix <- cbind(hierarchicalPalMatrix, grey.colors(n = nShades))

hierarchicalPal <- hierarchicalPalMatrix %>%
  as.vector() %>%
  setNames(unique(hierarchicalPalInfo$Taxa))
tax_palette_plot(hierarchicalPal) +
  theme(axis.text.y.left = element_text(species = "mono"))


## the normal 
sample_data(pseq2)$Age_Days <- factor(sample_data(pseq2)$Age_Days, levels = c("1", "10", "21"))

pseq2 %>%
  comp_barplot(tax_level = "Genus", n_taxa = 12, bar_width = 0.8) +
  theme(
    legend.text = element_text(family = "mono"),
    axis.text.y = element_blank(), 
    axis.ticks.y = element_blank(), 
    strip.background = element_rect(fill = "white")
  ) +
  coord_flip() +
  labs(x = NULL, y = NULL) +
  facet_wrap(~Treatment + Age_Days, scales = "free_y")

ggsave("plots/fungi_relative_abundance_BASIC_COLORED_TRT_plot1.pdf")

pseq2 %>%
  comp_barplot(tax_level = "Genus", n_taxa = 10, bar_width = 0.8) +
  theme(legend.text = element_text(family = "mono"),
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        strip.background = element_rect(fill = "white")) +
  coord_flip() + labs(x = NULL, y = NULL) + facet_wrap(~Age_Days, scales = "free_y")
ggsave("plots/fungi_relative_abundance_AGE_BASIC_COLORED_plot.pdf")
