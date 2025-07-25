#Differential relative abundance

set.seed(200789)
require(tidyverse)
require(phyloseq)
library(ALDEx2)
library(microViz)
library(writexl)
library(ggplot2)
library(ggrepel)
library(ggpubr)

#Load data
ps <- readRDS("ps/ps_taxfungi.RDS")
ps

ps <- subset_samples(ps, Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Basal Diet" | Treatment == "Essential oils")

# Inspect tables
head(otu_table(ps), n=3)
head(tax_table(ps), n=3)
head(sample_data(ps), n=3)

#transform the ps 
#transposed_otu_table <- t(otu_table)
#otu_table(ps) <- transposed_otu_table

dim(otu_table(ps))  # Check dimensions of the OTU table
head(otu_table(ps))  # View the first few rows of the OTU table

sample_data(ps)

#Two treatments only - it is a t-test  - Basal X Activo
##select two samples
ps <- subset_samples(ps, Treatment == "Basal Diet" | Treatment == "Probiotic")
sample_data(ps)

# aggregate at genus level for counts
ps1 <- ps %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

head(sample_data(ps1))

#transform to df
ps2 <- as.data.frame(otu_table(ps1))

# run AlDEx2 function
aldex2_da <- ALDEx2::aldex(ps2, phyloseq::sample_data(ps)$Treatment, taxa_rank = "all", gamma = 0.5, method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 1000, denom = "none")

sig_aldex2 <- aldex2_da %>%
  filter(wi.eBH < 0.05)

# setup tax table to be able to merge
taxa_info <- data.frame(tax_table(ps))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
#No differences 

#Check the other treatments - Basal x Magni
# load data
set.seed(200789)
ps <- readRDS("ps/ps_taxfungi.RDS")
ps

ps <- subset_samples(ps, Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Basal Diet" | Treatment == "Essential oils")

# Inspect tables
head(otu_table(ps), n=3)
head(tax_table(ps), n=3)
head(sample_data(ps), n=3)

ps <- subset_samples(ps, Treatment == "Basal Diet" | Treatment == "BMD")

sample_data(ps)
tail(sample_data(ps))

# aggregate at genus level for counts
ps1 <- ps %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

head(sample_data(ps1))

#transform to df
ps2 <- as.data.frame(otu_table(ps1))

# run AlDEx2 function
aldex2_da <- ALDEx2::aldex(ps2, phyloseq::sample_data(ps)$Treatment, taxa_rank = "all", gamma = 0.5, method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 1000, denom = "none")

sig_aldex2 <- aldex2_da %>%
  filter(wi.eBH < 0.05)

# setup tax table to be able to merge
taxa_info <- data.frame(tax_table(ps))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
#No differences 

#Check the other treatments - BMD X Prob
# load data
set.seed(200789)
ps <- readRDS("ps/ps_taxfungi.RDS")
ps

ps <- subset_samples(ps, Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Basal Diet" | Treatment == "Essential oils")

# Inspect tables
head(otu_table(ps), n=3)
head(tax_table(ps), n=3)
head(sample_data(ps), n=3)

ps <- subset_samples(ps, Treatment == "BMD" | Treatment == "Probiotic")

sample_data(ps)
tail(sample_data(ps))

# aggregate at genus level for counts
ps1 <- ps %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

head(sample_data(ps1))

#transform to df
ps2 <- as.data.frame(otu_table(ps1))

# run AlDEx2 function
aldex2_da <- ALDEx2::aldex(ps2, phyloseq::sample_data(ps)$Treatment, taxa_rank = "all", gamma = 0.5, method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 1000, denom = "none")

sig_aldex2 <- aldex2_da %>%
  filter(wi.eBH < 0.05)

# setup tax table to be able to merge
taxa_info <- data.frame(tax_table(ps))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
#No differences 

#Check the other treatments -  BMD X EO
# load data
set.seed(200789)
ps <- readRDS("ps/ps_taxfungi.RDS")
ps

ps <- subset_samples(ps, Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Basal Diet" | Treatment == "Essential oils")

# Inspect tables
head(otu_table(ps), n=3)
head(tax_table(ps), n=3)
head(sample_data(ps), n=3)

ps <- subset_samples(ps, Treatment == "BMD" | Treatment == "Essential oils")

sample_data(ps)
tail(sample_data(ps))

# aggregate at genus level for counts
ps1 <- ps %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

head(sample_data(ps1))

#transform to df
ps2 <- as.data.frame(otu_table(ps1))

# run AlDEx2 function
aldex2_da <- ALDEx2::aldex(ps2, phyloseq::sample_data(ps)$Treatment, taxa_rank = "all", gamma = 0.5, method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 1000, denom = "none")

sig_aldex2 <- aldex2_da %>%
  filter(wi.eBH < 0.05)

# setup tax table to be able to merge
taxa_info <- data.frame(tax_table(ps))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
#No differences 

#Check the other treatments - EO x Prob
# load data
set.seed(200789)
ps <- readRDS("ps/ps_taxfungi.RDS")
ps

ps <- subset_samples(ps, Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Basal Diet" | Treatment == "Essential oils")

# Inspect tables
head(otu_table(ps), n=3)
head(tax_table(ps), n=3)
head(sample_data(ps), n=3)

ps <- subset_samples(ps, Treatment == "Probiotic" | Treatment == "Essential oils")

sample_data(ps)
tail(sample_data(ps))

# aggregate at genus level for counts
ps1 <- ps %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

head(sample_data(ps1))

#transform to df
ps2 <- as.data.frame(otu_table(ps1))

# run AlDEx2 function
aldex2_da <- ALDEx2::aldex(ps2, phyloseq::sample_data(ps)$Treatment, taxa_rank = "all", gamma = 0.5, method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 1000, denom = "none")

sig_aldex2 <- aldex2_da %>%
  filter(wi.eBH < 0.05)

# setup tax table to be able to merge
taxa_info <- data.frame(tax_table(ps))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
#No differences 

#Check the other treatments - Basal x EO
# load data
set.seed(200789)
ps <- readRDS("ps/ps_taxfungi.RDS")
ps

ps <- subset_samples(ps, Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Basal Diet" | Treatment == "Essential oils")

# Inspect tables
head(otu_table(ps), n=3)
head(tax_table(ps), n=3)
head(sample_data(ps), n=3)

ps <- subset_samples(ps, Treatment == "Basal Diet" | Treatment == "Essential oils")

sample_data(ps)
tail(sample_data(ps))

# aggregate at genus level for counts
ps1 <- ps %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

head(sample_data(ps1))

#transform to df
ps2 <- as.data.frame(otu_table(ps1))

# run AlDEx2 function
aldex2_da <- ALDEx2::aldex(ps2, phyloseq::sample_data(ps)$Treatment, taxa_rank = "all", gamma = 0.5, method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 1000, denom = "none")

sig_aldex2 <- aldex2_da %>%
  filter(wi.eBH < 0.05)

# setup tax table to be able to merge
taxa_info <- data.frame(tax_table(ps))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
#No differences 

##AGE ###

set.seed(200789)
#Load data
ps <- readRDS("ps/ps_taxfungi.RDS")
ps

ps <- subset_samples(ps, Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Basal Diet" | Treatment == "Essential oils")

# Inspect tables
head(otu_table(ps), n=3)
head(tax_table(ps), n=3)
head(sample_data(ps), n=3)

#transform the ps 
#transposed_otu_table <- t(otu_table)
#otu_table(ps) <- transposed_otu_table

dim(otu_table(ps))  # Check dimensions of the OTU table
head(otu_table(ps))  # View the first few rows of the OTU table

sample_data(ps)

#Two treatments only - it is a t-test  - Basal X Activo
##select two samples
ps <- subset_samples(ps, Age_Days == "1" | Age_Days == "10")
sample_data(ps)

# aggregate at genus level for counts
ps1 <- ps %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

head(sample_data(ps1))

#transform to df
ps2 <- as.data.frame(otu_table(ps1))

# run AlDEx2 function
aldex2_da <- ALDEx2::aldex(ps2, phyloseq::sample_data(ps)$Age_Days, taxa_rank = "all", gamma = 0.5, method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 1000, denom = "none")

sig_aldex2 <- aldex2_da %>%
  filter(wi.eBH < 0.05)

# setup tax table to be able to merge
taxa_info <- data.frame(tax_table(ps))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
#no differences

#Check the other treatments - Basal x Magni
# load data
set.seed(200789)
ps <- readRDS("ps/ps_taxfungi.RDS")
ps

ps <- subset_samples(ps, Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Basal Diet" | Treatment == "Essential oils")

# Inspect tables
head(otu_table(ps), n=3)
head(tax_table(ps), n=3)
head(sample_data(ps), n=3)

ps <- subset_samples(ps, Age_Days == "1" | Age_Days == "21")

sample_data(ps)
tail(sample_data(ps))

# aggregate at genus level for counts
ps1 <- ps %>%
  tax_fix(unknowns = c("Incertae Sedis")) %>%
  tax_fix()

head(sample_data(ps1))

#transform to df
ps2 <- as.data.frame(otu_table(ps1))

# run AlDEx2 function
aldex2_da <- ALDEx2::aldex(ps2, phyloseq::sample_data(ps)$Age_Days, taxa_rank = "all", gamma = 0.5, method = "t.test", p_adjust = "BH", pvalue_cutoff = 0.05, mc_samples = 1000, denom = "none")

sig_aldex2 <- aldex2_da %>%
  filter(wi.eBH < 0.05)

# setup tax table to be able to merge
taxa_info <- data.frame(tax_table(ps))
taxa_info <- taxa_info %>% rownames_to_column(var = "OTU")

# look at plot to see if positive/negative effect size corresponds to which treatment level
ALDEx2::aldex.plot(aldex2_da, type="MW", test="wilcox", called.cex = 1, cutoff = 0.05)
#There NO differences between 1 and 10
#there are differences between 1 and 21 

# make a table of significant corrected p-values
sig_aldex2 <- aldex2_da %>%
  rownames_to_column(var = "OTU") %>%
  filter(wi.eBH < 0.05) %>%
  arrange(effect, wi.eBH) %>%
  dplyr::select(OTU, diff.btw, diff.win, effect, wi.ep, wi.eBH)

# add in previously formed taxa information to complete the table
sig_aldex2 <- left_join(sig_aldex2, taxa_info)
save(sig_aldex2, file = "tables/table_sig_taxa_age1_21.RData")
write.table(sig_aldex2, file = "tables/table_sig_taxa_age1_21.csv", sep = ",", col.names = TRUE, row.names = FALSE)
write_xlsx(sig_aldex2, "tables/table_sig_taxa_age1-21.xlsx")


##graph 1
dat1 <- readxl::read_excel("tables/table_sig_taxa_age1-21.xlsx")

ggplot(data=dat1, aes(x=Species, y=(effect), col=Age, label=Species)) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  #scale_color_brewer() +
  #geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=0, col="black")


###Lollipop graph
ggplot(dat1, aes(x = effect, y = Genus)) +
  geom_segment(aes(yend = Genus), xend = 0, colour = "grey50") +
  geom_point(size = 3, aes(colour = Age)) +
  # scale_colour_brewer(palette = "Set1", limits = c("Starter", "Grower"), guide = "none") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank()) +
  facet_grid(Age ~ ., scales = "free_y", space = "free_y")


##cleavand dot plot graph
p1 <- ggplot(dat1, aes(effect, Genus_Species)) +
  geom_line(aes(group = Genus_Species)) +
  geom_point(aes(color = factor(Age)), size = 8) +  # Ensure Age is categorical
  scale_color_manual(values = c("#996699", "#333399")) +
  labs(y = "Species", color = "Age (Days)") +     # Change axis and legend labels
  theme_bw(base_size = 15) +
  theme(
    axis.title.x = element_text(size = 20),  # X-axis title font size
    axis.title.y = element_text(size = 20),  # Y-axis title font size
    axis.text.x = element_text(size = 18),   # X-axis tick labels
    axis.text.y = element_text(size = 18),
    legend.title = element_text(size = 20))

p2 <- ggpar(p1, font.ytickslab = "italic")
p2




