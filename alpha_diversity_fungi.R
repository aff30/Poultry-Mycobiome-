## Alpha diversity 

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
library(ggpubr)
library(car)

#Load data
ps <- readRDS("ps/ps_taxfungi.RDS")
ps

ps <- subset_samples(ps, Treatment == "Probiotic" | Treatment == "BMD" | Treatment == "Basal Diet" | Treatment == "Essential oils")

## subset the treatments
sample_data(ps)

# get sample data
sampdf <- sample_data(ps) %>% 
  data.frame() %>% 
  rownames_to_column(var = "Sample_ID") 

# get diversity
alpha <- estimate_richness(ps, measures = c("Shannon", "Simpson", "Observed")) %>% 
  # make ID
  rownames_to_column(var = "Sample_ID") 
#alpha$Sample_ID <- str_remove(alpha$Sample_ID, "X")
#str_remove(alpha$Project_ID,"X")
# merge(sampdf, alpha, by = "Project_ID") 
alpha <- merge(sampdf, alpha, by = "Sample_ID") 

## ---- Shannon Kruskal Wallis  ----
#test for normality 
ggdensity(alpha$Shannon, 
          main = "Density plot of diversity",
          xlab = "diversity")

ggqqplot(alpha$Shannon)
qqPlot(alpha$Shannon)

#Not normal but just to confirm 
shapiro.test(alpha$Shannon)
#data:  alpha$Shannon
#W = 0.93039, p-value = 6.765e-05

kruskalShannonTRT <- kruskal.test(Shannon ~ Treatment, data = alpha) 
kruskalShannonTRT
# p-value = 0.5978

kruskalObservedTRT <- kruskal.test(Observed ~ Treatment, data = alpha)
kruskalObservedTRT
# p-value = 0.571


#Age 
kruskalShannonAGE <- kruskal.test(Shannon ~ Age_Days, data = alpha) 
kruskalShannonAGE
# p-value = 5.594e-09

kruskalObservedAGE <- kruskal.test(Observed ~ Age_Days, data = alpha)
kruskalObservedAGE
# p-value = 8.004e-08


#colors
trts <- c("Basal Diet" = "#996666", 
          "Probiotic" = "#24868EFF",
          "Essential oils" = "#35B779FF",
          "BMD" = "#C7E020FF")


AGE <- c("1" ="#996699", 
         "10" = "#006666", 
         "21" = "#333399")

#boxplots

shannonbox <- 
  ggplot(alpha,
         aes(x = Treatment, y = Shannon, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = trts) +
  xlab("") + geom_point(position = "jitter", color = "black", aes(shape = Treatment), size = 3) +
  theme_classic() + 
  labs(title = "",
       x = "", 
       y = "Shannon's Index",
       shape = "Treatment",
       fill = "Treatment") + 
  #scale_x_discrete(name="", labels = c("cloaca" = "Cloaca", "sock"="Sock")) +
  theme(axis.text= element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none")

shannonbox

observedbox <- 
  ggplot(alpha, aes(x = Treatment, y = Observed, fill = Treatment)) +
  geom_boxplot() +
  scale_fill_manual(values = trts) +
  xlab("") +
  geom_point(position = "jitter", color = "black", aes(shape = Treatment), size = 3) +
  theme_classic() + 
  labs(title = "",
       x = "", 
       y = "Observed ASVs",
       shape = "Treatment",
       fill = "Treatment") + 
  #scale_x_discrete(name="", labels = c("cloaca" = "Cloaca", "sock"="Sock")) +
  theme(axis.text= element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none")

observedbox


##Age Days 

shannonboxage <- 
  ggplot(alpha,
         aes(x = Age_Days, y = Shannon, fill = Age_Days)) +
  geom_boxplot() +
  scale_fill_manual(values = AGE) +
  xlab("") + geom_point(position = "jitter", color = "black", aes(shape = Treatment), size = 3) +
  theme_classic() + 
  labs(title = "",
       x = "", 
       y = "Shannon's Index",
       shape = "Age_Days",
       fill = "Age_Days") + 
  #scale_x_discrete(name="", labels = c("cloaca" = "Cloaca", "sock"="Sock")) +
  theme(axis.text= element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none")

shannonboxage


observedboxage <- 
  ggplot(alpha, aes(x = Age_Days, y = Observed, fill = Age_Days)) +
  geom_boxplot() +
  scale_fill_manual(values = AGE) +
  xlab("") +
  geom_point(position = "jitter", color = "black", aes(shape = Treatment), size = 3) +
  theme_classic() + 
  labs(title = "",
       x = "", 
       y = "Observed ASVs",
       shape = "Age_Days",
       fill = "Age_Days") + 
  #scale_x_discrete(name="", labels = c("cloaca" = "Cloaca", "sock"="Sock")) +
  theme(axis.text= element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        plot.title = element_text(hjust = 0.5)) +
  guides(fill = "none")

observedboxage


