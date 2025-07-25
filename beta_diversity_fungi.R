## Beta diversity 

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
library(pairwiseAdonis)

#Load data
ps <- readRDS("ps/ps_taxfungi.RDS")
ps

ps <- subset_samples(ps, Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Basal Diet" | Treatment == "Essential oils")

#principal coordinate analysis

#colors
trts <- c("Basal Diet" = "#996666", 
          "Probiotic" = "#24868EFF",
          "Essential oils" = "#35B779FF",
          "BMD" = "#C7E020FF")

AGE <- c("1" ="#996699", 
         "10" = "#006666", 
         "21" = "#333399")

## Beta Diversity
ps <- tax_fix(ps)
ps <- phyloseq_validate(ps, remove_undetected = TRUE)

ps %>% 
  tax_fix(unknowns = c("s__kudriavzevii", "s__oryzae", "s__phaffii")) %>% 
  tax_transform(trans = "clr", rank = "Species")

ps %>% 
  tax_transform(trans = "identity", rank = "unique") %>% 
  dist_calc("aitchison") 

sample_data(ps)

# Plot a PCA
# Principal Components Analysis is an unconstrained method that does not use a distance matrix.
# Each point is a sample, and samples that appear closer together are typically 
# more similar to each other than samples which are further apart.


ps %>% tax_fix() %>% 
  tax_fix(unknowns = c("s__kudriavzevii", "s__oryzae", "s__phaffii")) %>% 
  tax_transform("clr", rank = "Species") %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>% 
  ord_plot(color = "Age_Days", shape = "Treatment", size = 8) +
  stat_ellipse(aes(group = Age_Days, color=Age_Days, linetype = Age_Days)) +
  scale_linetype_manual(values = c("dashed", "solid", "dashed"))+
  theme_classic() +
  scale_color_manual(values = c("#996699", "#006666", "#333399")) +
  scale_shape_manual(values = c(16, 15, 23, 17)) +
  theme_classic()

ggsave("plots/fungi_betadiversity_plot1.pdf")


# roughly estimate which samples will contain more of that taxon 
ps %>% tax_fix() %>% 
  tax_fix(unknowns = c("s__kudriavzevii", "s__oryzae", "s__phaffii")) %>% 
  tax_transform("clr", rank = "Species") %>% 
  # when no distance matrix or constraints are supplied, PCA is the default/auto ordination method
  ord_calc() %>% 
  ord_plot(color = "Age_Days", plot_taxa = 1:7, size = 6, 
           taxa.text.size = 8, parse = TRUE) +
  stat_ellipse(aes(group = Age_Days, color=Age_Days, linetype = Age_Days)) +
  scale_linetype_manual(values = c("dashed", "solid", "dashed"))+
  scale_color_manual(values = c("#996699", "#006666", "#333399")) +
  scale_shape_manual(values = c(16, 15, 23, 17)) +
  theme_classic() +
  theme(axis.text= element_text(size = 22),
        legend.text = element_text(size = 22),
        legend.title = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.title.x = element_text(size = 22)
  )

ggsave(ggsave("plots/fungi_plot_PCA_age1-10-21.pdf", width = 25, height = 20, units = "cm"))

## PCoA
# Aitchison distance

ps %>% tax_fix() %>% 
  tax_transform("identity", rank = "Class") %>% # don't transform!
  dist_calc("aitchison") %>% 
  ord_calc("PCoA") %>% 
  ord_plot(color = "Treatment", shape = "Age_Days", size = 7) +
  stat_ellipse(aes(group = Treatment, color=Treatment, linetype = Treatment)) +
  scale_linetype_manual(values = c("dashed", "solid", "dashed", "solid"))+
  scale_color_manual(values = trts, labels =  c("Basal Diet" = "Basal Diet", "Essential oils"="Essential Oils", "Probiotic" = "Probiotic", "BMD" = "Antibiotic")) +
  theme_classic() 



ggsave("plots/fungi_treatment_MDS_plot24.pdf")

#PERMANOVA overall trt and age 
# clr transform phyloseq objects at Genus level
beta1 <- ps %>% 
  tax_fix(unknowns = c("s__kudriavzevii", "s__oryzae", "s__phaffii")) %>% 
  tax_transform(trans = "clr", rank = "Species") %>% 
  ps_get()

# generate distance matrix
psdist <- phyloseq::distance(beta1, method = "euclidean")

#ADONIS test
# age, treatment

vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Treatment, permutations = 10000) 
# p = 0.9799
vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Age_Days, permutations = 10000)
# p = 9.999e-05 ***

# Pairwise comparison using pairwiseAdonis 
data <- phyloseq::sample_data(beta1)
data
pairwiseAdonis::pairwise.adonis(psdist, data$Age_Days)

#pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
#1  1 vs 10  1  522.9703  8.876471 0.12349620   0.001      0.003   *
 # 2  1 vs 21  1  608.4942 10.743615 0.14769152   0.001      0.003   *
 # 3 10 vs 21  1  242.5319  3.595656 0.05399235   0.001      0.003   *

# testing if age is masking treatment
# ADONIS with both Age and Treatment
adonis_result <- vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Treatment + phyloseq::sample_data(beta1)$Age_Days, permutations = 10000)
print(adonis_result)

#Df SumOfSqs      R2      F    Pr(>F)    
#Model     5   2085.9 0.19153 4.3117 9.999e-05 ***
 # Residual 91   8804.9 0.80847                     
#Total    96  10890.9 1.00000  

vegan::adonis2(psdist ~ phyloseq::sample_data(beta1)$Treatment * phyloseq::sample_data(beta1)$Age_Days, permutations = 10000)
#Df SumOfSqs      R2      F    Pr(>F)    
#Model    11   2676.0 0.24571 2.5172 9.999e-05 ***
 # Residual 85   8214.9 0.75429                     
#Total    96  10890.9 1.00000

#The significant p-value suggests that the interaction between Treatment and Age_Days is meaningful in explaining variation in the mycobiome data.

# Assuming 'psdist' is your distance matrix and the groups are based on Age_Days and Treatment
betadisper_result <- betadisper(psdist, group = interaction(phyloseq::sample_data(beta1)$Age_Days, phyloseq::sample_data(beta1)$Treatment))

# To view the centroids and distances:
betadisper_result$centroids  # this gives the centroid locations for each group
betadisper_result$distances  # this gives the distances to the centroid for each sample

# Assuming 'betadisper_result$centroids' is a matrix or data frame with centroid coordinates
centroids <- betadisper_result$centroids
dist_centroids <- dist(centroids)
dist_centroids

#Calculating centroids - June 2025

library(phyloseq)
library(dplyr)
library(ggplot2)

# Use your existing Aitchison distance (Euclidean on CLR)
ord <- ape::pcoa(psdist)

ord_df <- as.data.frame(ord$vectors[, 1:2])
ord_df$SampleID <- rownames(ord_df)

meta_df <- as(sample_data(beta1), "data.frame")  # fix here
meta_df$SampleID <- rownames(meta_df)

plot_df <- left_join(ord_df, meta_df, by = "SampleID")
colnames(plot_df)[1:2] <- c("PCoA1", "PCoA2")

#By age
centroids_age <- plot_df %>%
  group_by(Age_Days) %>%
  summarise(Centroid1 = mean(PCoA1), Centroid2 = mean(PCoA2))

# By TRT
centroids_trt <- plot_df %>%
  group_by(Treatment) %>%
  summarise(Centroid1 = mean(PCoA1), Centroid2 = mean(PCoA2))

#Plot with samples, centroids, and distance lines

#Age centroids and distance lines
plot_df_age <- left_join(plot_df, centroids_age, by = "Age_Days")

ggplot(plot_df_age, aes(x = PCoA1, y = PCoA2, color = factor(Age_Days))) +
  geom_segment(aes(xend = Centroid1, yend = Centroid2), alpha = 1) +  # lines to centroid
  geom_point(size = 7, alpha = 1) +  # sample points
  geom_point(aes(x = Centroid1, y = Centroid2), shape = 1, size = 7, color = "black") +  # centroids
  geom_text(data = centroids_age, aes(x = Centroid1, y = Centroid2, label = Age_Days),
            vjust = -1, size = 7, fontface = "bold", color = "black") +
  theme_minimal() +
  scale_color_manual(values = c("#996699", "#006666", "#333399")) +
  labs(color = "Age (Days)") +
  theme(
    axis.title.x = element_text(size = 25),  # X-axis title font size
    axis.title.y = element_text(size = 24),  # Y-axis title font size
    axis.text.x = element_text(size = 19),   # X-axis tick labels
    axis.text.y = element_text(size = 19),
    legend.title = element_text(size = 25),
    legend.key.size = unit(2, "lines"),
    legend.text = element_text(size = 24))  

# Treatment centroids and distance lines
plot_df_trt <- left_join(plot_df, centroids_trt, by = "Treatment")


ggplot(plot_df_trt, aes(x = PCoA1, y = PCoA2, color = Treatment)) +
  geom_segment(aes(xend = Centroid1, yend = Centroid2), alpha = 1) +
  geom_point(size = 6, alpha = 1) +
  geom_point(aes(x = Centroid1, y = Centroid2), shape = 1, size = 9, color = "black") +
  geom_text(data = centroids_trt, aes(x = Centroid1, y = Centroid2, label = Treatment),
            vjust = -1.8, size = 5, fontface = "bold", color = "black") +
  scale_color_manual(values = trts, labels =  c("Basal Diet" = "Basal Diet", "Essential oils"="Essential Oils", "Probiotic" = "Probiotic", "BMD" = "Antibiotic")) +
  labs(color = "Treatment") +
  theme(
    axis.title.x = element_text(size = 25),  # X-axis title font size
    axis.title.y = element_text(size = 24),  # Y-axis title font size
    axis.text.x = element_text(size = 19),   # X-axis tick labels
    axis.text.y = element_text(size = 19),
    legend.title = element_text(size = 25),
    legend.key.size = unit(2, "lines"),
    legend.text = element_text(size = 22)) + 
  theme_minimal()
#If centroids for two groups are far apart, it suggests a strong composition difference between those groups.
#If centroids are close together, groups are similar in composition or function

#subset by age and plot centroids by tRt
#d1
ps <- subset_samples(ps, Age_Days == "1")
sample_data(ps)

# clr transform phyloseq objects at Genus level
beta1 <- ps %>% 
  tax_fix(unknowns = c("s__kudriavzevii", "s__oryzae", "s__phaffii")) %>% 
  tax_transform(trans = "clr", rank = "Species") %>% 
  ps_get()

# generate distance matrix
psdist <- phyloseq::distance(beta1, method = "euclidean")

# Assuming 'psdist' is your distance matrix and the groups are based on Age_Days and Treatment
betadisper_result <- betadisper(psdist, group = interaction(phyloseq::sample_data(beta1)$Age_Days, phyloseq::sample_data(beta1)$Treatment))

# To view the centroids and distances:
betadisper_result$centroids  # this gives the centroid locations for each group
betadisper_result$distances  # this gives the distances to the centroid for each sample

# Assuming 'betadisper_result$centroids' is a matrix or data frame with centroid coordinates
centroids <- betadisper_result$centroids
dist_centroids <- dist(centroids)
dist_centroids

#Calculating centroids - June 2025

library(phyloseq)
library(dplyr)
library(ggplot2)

# Use your existing Aitchison distance (Euclidean on CLR)
ord <- ape::pcoa(psdist)

ord_df <- as.data.frame(ord$vectors[, 1:2])
ord_df$SampleID <- rownames(ord_df)

meta_df <- as(sample_data(beta1), "data.frame")  # fix here
meta_df$SampleID <- rownames(meta_df)

plot_df <- left_join(ord_df, meta_df, by = "SampleID")
colnames(plot_df)[1:2] <- c("PCoA1", "PCoA2")

#By age
centroids_age <- plot_df %>%
  group_by(Age_Days) %>%
  summarise(Centroid1 = mean(PCoA1), Centroid2 = mean(PCoA2))

# By TRT
centroids_trt <- plot_df %>%
  group_by(Treatment) %>%
  summarise(Centroid1 = mean(PCoA1), Centroid2 = mean(PCoA2))


# Treatment centroids and distance lines
plot_df_trt <- left_join(plot_df, centroids_trt, by = "Treatment")


ggplot(plot_df_trt, aes(x = PCoA1, y = PCoA2, color = Treatment)) +
  geom_segment(aes(xend = Centroid1, yend = Centroid2), alpha = 1) +
  geom_point(size = 6, alpha = 1) +
  geom_point(aes(x = Centroid1, y = Centroid2), shape = 1, size = 9, color = "black") +
  geom_text(data = centroids_trt, aes(x = Centroid1, y = Centroid2, label = Treatment),
            vjust = -1.8, size = 5, fontface = "bold", color = "black") +
  scale_color_manual(values = trts, labels =  c("Basal Diet" = "Basal Diet", "Essential oils"="Essential Oils", "Probiotic" = "Probiotic", "BMD" = "Antibiotic")) +
  labs(color = "Treatment") +
  theme(
    axis.title.x = element_text(size = 25),  # X-axis title font size
    axis.title.y = element_text(size = 24),  # Y-axis title font size
    axis.text.x = element_text(size = 22),   # X-axis tick labels
    axis.text.y = element_text(size = 22),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 22),
    axis.text= element_text(size = 22)
    ) + 
  theme_classic2()
#If centroids for two groups are far apart, it suggests a strong compositional difference between those groups.
#If centroids are close together, groups are similar in composition or function


#subset by age and plot centroids by tRt
ps <- readRDS("ps/ps_taxfungi.RDS")
ps

ps <- subset_samples(ps, Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Basal Diet" | Treatment == "Essential oils")

#d10
ps <- subset_samples(ps, Age_Days == "10")
sample_data(ps)

# clr transform phyloseq objects at Genus level
beta1 <- ps %>% 
  tax_fix(unknowns = c("s__kudriavzevii", "s__oryzae", "s__phaffii")) %>% 
  tax_transform(trans = "clr", rank = "Species") %>% 
  ps_get()

# generate distance matrix
psdist <- phyloseq::distance(beta1, method = "euclidean")

# Assuming 'psdist' is your distance matrix and the groups are based on Age_Days and Treatment
betadisper_result <- betadisper(psdist, group = interaction(phyloseq::sample_data(beta1)$Age_Days, phyloseq::sample_data(beta1)$Treatment))

# To view the centroids and distances:
betadisper_result$centroids  # this gives the centroid locations for each group
betadisper_result$distances  # this gives the distances to the centroid for each sample

# Assuming 'betadisper_result$centroids' is a matrix or data frame with centroid coordinates
centroids <- betadisper_result$centroids
dist_centroids <- dist(centroids)
dist_centroids

#Calculating centroids - June 2025

library(phyloseq)
library(dplyr)
library(ggplot2)

# Use your existing Aitchison distance (Euclidean on CLR)
ord <- ape::pcoa(psdist)

ord_df <- as.data.frame(ord$vectors[, 1:2])
ord_df$SampleID <- rownames(ord_df)

meta_df <- as(sample_data(beta1), "data.frame")  # fix here
meta_df$SampleID <- rownames(meta_df)

plot_df <- left_join(ord_df, meta_df, by = "SampleID")
colnames(plot_df)[1:2] <- c("PCoA1", "PCoA2")

#By age
centroids_age <- plot_df %>%
  group_by(Age_Days) %>%
  summarise(Centroid1 = mean(PCoA1), Centroid2 = mean(PCoA2))

# By TRT
centroids_trt <- plot_df %>%
  group_by(Treatment) %>%
  summarise(Centroid1 = mean(PCoA1), Centroid2 = mean(PCoA2))


# Treatment centroids and distance lines
plot_df_trt <- left_join(plot_df, centroids_trt, by = "Treatment")


ggplot(plot_df_trt, aes(x = PCoA1, y = PCoA2, color = Treatment)) +
  geom_segment(aes(xend = Centroid1, yend = Centroid2), alpha = 1) +
  geom_point(size = 6, alpha = 1) +
  geom_point(aes(x = Centroid1, y = Centroid2), shape = 1, size = 9, color = "black") +
  geom_text(data = centroids_trt, aes(x = Centroid1, y = Centroid2, label = Treatment),
            vjust = -1.8, size = 5, fontface = "bold", color = "black") +
  scale_color_manual(values = trts, labels =  c("Basal Diet" = "Basal Diet", "Essential oils"="Essential Oils", "Probiotic" = "Probiotic", "BMD" = "Antibiotic")) +
  labs(color = "Treatment") +
  theme(
    axis.title.x = element_text(size = 25),  # X-axis title font size
    axis.title.y = element_text(size = 24),  # Y-axis title font size
    axis.text.x = element_text(size = 22),   # X-axis tick labels
    axis.text.y = element_text(size = 22),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 22),
    axis.text= element_text(size = 22)
  ) + 
  theme_classic2()
#If centroids for two groups are far apart, it suggests a strong compositional difference between those groups.
#If centroids are close together, groups are similar in composition or function


#subset by age and plot centroids by tRt
ps <- readRDS("ps/ps_taxfungi.RDS")
ps

ps <- subset_samples(ps, Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Basal Diet" | Treatment == "Essential oils")

#d21
ps <- subset_samples(ps, Age_Days == "21")
sample_data(ps)

# clr transform phyloseq objects at Genus level
beta1 <- ps %>% 
  tax_fix(unknowns = c("s__kudriavzevii", "s__oryzae", "s__phaffii")) %>% 
  tax_transform(trans = "clr", rank = "Species") %>% 
  ps_get()

# generate distance matrix
psdist <- phyloseq::distance(beta1, method = "euclidean")

# Assuming 'psdist' is your distance matrix and the groups are based on Age_Days and Treatment
betadisper_result <- betadisper(psdist, group = interaction(phyloseq::sample_data(beta1)$Age_Days, phyloseq::sample_data(beta1)$Treatment))

# To view the centroids and distances:
betadisper_result$centroids  # this gives the centroid locations for each group
betadisper_result$distances  # this gives the distances to the centroid for each sample

# Assuming 'betadisper_result$centroids' is a matrix or data frame with centroid coordinates
centroids <- betadisper_result$centroids
dist_centroids <- dist(centroids)
dist_centroids

#Calculating centroids - June 2025

library(phyloseq)
library(dplyr)
library(ggplot2)

# Use your existing Aitchison distance (Euclidean on CLR)
ord <- ape::pcoa(psdist)

ord_df <- as.data.frame(ord$vectors[, 1:2])
ord_df$SampleID <- rownames(ord_df)

meta_df <- as(sample_data(beta1), "data.frame")  # fix here
meta_df$SampleID <- rownames(meta_df)

plot_df <- left_join(ord_df, meta_df, by = "SampleID")
colnames(plot_df)[1:2] <- c("PCoA1", "PCoA2")

#By age
centroids_age <- plot_df %>%
  group_by(Age_Days) %>%
  summarise(Centroid1 = mean(PCoA1), Centroid2 = mean(PCoA2))

# By TRT
centroids_trt <- plot_df %>%
  group_by(Treatment) %>%
  summarise(Centroid1 = mean(PCoA1), Centroid2 = mean(PCoA2))


# Treatment centroids and distance lines
plot_df_trt <- left_join(plot_df, centroids_trt, by = "Treatment")


ggplot(plot_df_trt, aes(x = PCoA1, y = PCoA2, color = Treatment)) +
  geom_segment(aes(xend = Centroid1, yend = Centroid2), alpha = 1) +
  geom_point(size = 6, alpha = 1) +
  geom_point(aes(x = Centroid1, y = Centroid2), shape = 1, size = 9, color = "black") +
  geom_text(data = centroids_trt, aes(x = Centroid1, y = Centroid2, label = Treatment),
            vjust = -1.8, size = 5, fontface = "bold", color = "black") +
  scale_color_manual(values = trts, labels =  c("Basal Diet" = "Basal Diet", "Essential oils"="Essential Oils", "Probiotic" = "Probiotic", "BMD" = "Antibiotic")) +
  labs(color = "Treatment") +
  theme(
    axis.title.x = element_text(size = 25),  # X-axis title font size
    axis.title.y = element_text(size = 24),  # Y-axis title font size
    axis.text.x = element_text(size = 22),   # X-axis tick labels
    axis.text.y = element_text(size = 22),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 22),
    axis.text= element_text(size = 22)
  ) + 
  theme_classic2()
#If centroids for two groups are far apart, it suggests a strong compositional difference between those groups.
#If centroids are close together, groups are similar in composition or function


