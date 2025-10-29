# Functional Anno 
# Ana - May 2025 - adapted by SK code

# ---- Load packages ----
library(tidyverse)
library(tidylog)
#if (!require("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(microViz)
library(phyloseq)

# ---- Load tables and combine ----
path <- "ooutput_funcprofiler_AF/"

# Get file list
file_list <- list.files(path = path, pattern = "\\.tsv$", full.names = TRUE)

# Loop through and combine
combined <- lapply(file_list, function(file) {
  # Read the file
  df <- read.delim(file, sep = "\t")
  # Add a column with the filename
  df$file <- basename(file)
  return(df)
}) %>%
  bind_rows()  # Combine all dataframes into one

# quick check bc i'm curious if the abundance is a relative abundance within sample 

combined %>% 
  group_by(file) %>%
  mutate(total = sum(abundance)) %>%
  select(file, total) %>%
  unique() # abundance from the tool is relative 


# ---- Clean names and add KO information and metadata----
kocat <- read.delim("ko_ref/ko_pathwayref.txt", header=F, sep = "\t")
komap <- read.csv("ko_ref/kegg_maps_table.csv", quote = "")
koterm <- read.delim("ko_ref/ko_termtopath.txt", header = F, sep = "\t")
met<- read.csv("met_Chicken/Functional_metadata.csv")

combined2 <- combined %>%
  mutate(Sample_ID = str_remove_all(file, "_merged.non.host_temp.tsv")) %>%

  # add metadata
  inner_join(met, by = "Sample_ID") %>%
  # add ko information
  inner_join(koterm %>%
               mutate(ko_id = V2,
                      path_id = V1) %>%
               select(ko_id, path_id) %>%
               filter(str_detect(path_id, "map")) %>%
               mutate(map_ID = str_remove_all(path_id, "path:map")) %>%
               mutate(path_id = str_remove_all(path_id, "path:")), by = "ko_id") %>%
  inner_join(kocat %>%
               mutate(path_id=V1,
                      type1=V2) %>%
               select(path_id, type1), by = "path_id") %>%
  inner_join(komap %>%
               mutate(map_ID = as.character(map_ID)), by = "map_ID") %>%
  # fix some formatting things
  mutate(type = str_remove_all(type, '"')) %>%
  mutate(type = str_split_fixed(type, "\\(http", 2)[,1]) %>%
  # reorder columns for my sanity 
  select(Sample_ID, Age_Days, Treatment, Cage, Phase, ko_id, abundance, path_id, map_ID, broadgroup, subgroup, type, type1, flags)

sample_data(combined2)

# Example: saving a data frame to CSV
write.csv(combined2, "combined_functional.csv", row.names = FALSE)


# I think from here the most informative thing would be to do differentail relative abundance on the ko_ids and then 
# plot that using the available ko metadata - you could even do that first then add the metadata later. 

## Relative Abundance


