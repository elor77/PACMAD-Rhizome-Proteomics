# Script Title: Summarize Cross-Species Differentially Abundant Proteins by Orthogroups
# 
# Description: Loads differential abundance summaries for six species and orthogroup 
#              mappings. Summarizes log2 fold change and p-values at the orthogroup level, 
#              ranks proteins and orthogroups within species, and outputs long-format tables 
#              for cross-species correlation analyses.
#
# Author: Elad Oren
# Created: 2025-01-01
# Last Modified: 2025-04-27
#
# Input: 
#   - 0_Orthogroups_Ag_Mg_Pv_Sc_Sn_Td-v1.0.txt - Orthogroup mapping file
#   - 1_Ag_summary_log2FC_DAP.csv - Andropogon gerardii DAP summary
#   - 1_Mg_summary_log2FC_DAP.csv - Miscanthus giganteus DAP summary
#   - 1_Pv-DD_summary_log2FC_DAP.csv - Panicum virgatum DAP summary
#   - 1_Sn_summary_log2FC_DAP.csv - Sorghastrum nutans DAP summary
#   - 1_Trip2019_summary_log2FC_DAP.csv - Tripsacum 2019 DAP summary
#   - 1_Trip2022_summary_log2FC_DAP.csv - Tripsacum 2022 DAP summary
#
# Output:
#   - 3a_summary-by-OG_long_log2FC_DAP.csv - Long-format summary by orthogroup
#   - 3b_long_log2FC_DAP_rank_protein_id.csv - Long-format summary by protein ID
#
# Dependencies: dplyr, readr, tidyr
#
# Notes: Summarizes protein log2 fold changes by orthogroup using mean values and 
#        selects minimum adjusted p-values per group. Classifies differentially abundant 
#        proteins using thresholds of |log2FC| ≥ 1 and padj ≤ 0.05. Provides ranked 
#        lists of proteins and orthogroups within each species for downstream comparative analyses.


# Reset
rm(list = ls())

# Libraries
library(dplyr)
library(readr)
library(tidyr)

setwd("~/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/1_tables")

# Step 1: Load orthogroup mapping
orthogroups <- read_delim("0_Orthogroups_Ag_Mg_Pv_Sc_Sn_Td-v1.0.txt")

# Step 2: Load each DAP summary
Ag_raw <- read_csv("1_Ag_summary_log2FC_DAP.csv")
Mg_raw <- read_csv("1_Mg_summary_log2FC_DAP.csv")
Pv_raw <- read_csv("1_Pv-DD_summary_log2FC_DAP.csv")
Sn_raw <- read_csv("1_Sn_summary_log2FC_DAP.csv")
Trip19_raw <- read_csv("1_Trip2019_summary_log2FC_DAP.csv")
Trip22_raw <- read_csv("1_Trip2022_summary_log2FC_DAP.csv")

##Create long format summary by orthogroups for corss species correlation matrix

# Step 3: Define a function to summarize each species
summarize_species <- function(df, species_name) {
  df %>%
    left_join(orthogroups, by = "protein_id") %>%
    group_by(orthogroup) %>%
    summarise(
      log2FC = mean(log2FC, na.rm = TRUE),
      padj = min(padj, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      species = species_name,
      DAP = case_when(
        padj < 0.05 & log2FC > 1 ~ "UP",
        padj < 0.05 & log2FC < -1 ~ "DOWN",
        TRUE ~ "NS"
      )
    )
}

# Step 4: Summarize each species separately
Ag <- summarize_species(Ag_raw, "Ag")
Mg <- summarize_species(Mg_raw, "Mg")
Pv <- summarize_species(Pv_raw, "Pv")
Sn <- summarize_species(Sn_raw, "Sn")
Trip19 <- summarize_species(Trip19_raw, "Trip19")
Trip22 <- summarize_species(Trip22_raw, "Trip22")

# Step 5: Merge all into one big long table
all_long <- bind_rows(Ag, Mg, Pv, Sn, Trip19, Trip22)

# Step 6: Rank log2FC within each species
all_long <- all_long %>%
  group_by(species) %>%
  mutate(rank = rank(-log2FC, ties.method = "min", na.last = "keep")) %>%
  ungroup()

# Step 7: Save long-format table
write_csv(all_long, "3a_summary-by-OG_long_log2FC_DAP.csv")



## Create a long format with protein IDs and orthogroups
Ag_raw <- Ag_raw %>%
  left_join(orthogroups %>%
              select(protein_id, orthogroup),
            by = "protein_id")
Mg_raw <- Mg_raw %>% left_join(orthogroups %>%
                                 select(protein_id, orthogroup),
                               by = "protein_id")
Pv_raw <- Pv_raw %>% left_join(orthogroups %>%
                                 select(protein_id, orthogroup),
                               by = "protein_id")
Sn_raw <- Sn_raw %>% left_join(orthogroups %>%
                                 select(protein_id, orthogroup),
                               by = "protein_id")
Trip19_raw <- Trip19_raw %>% left_join(orthogroups %>%
                                         select(protein_id, orthogroup),
                                       by = "protein_id")
Trip22_raw <- Trip22_raw %>% left_join(orthogroups %>%
                                         select(protein_id, orthogroup),
                                       by = "protein_id")

# Combine all species into one long table
all_long_id <- bind_rows(
  Ag_raw %>% mutate(species = "Ag"),
  Mg_raw %>% mutate(species = "Mg"),
  Pv_raw %>% mutate(species = "Pv"),
  Sn_raw %>% mutate(species = "Sn"),
  Trip19_raw %>% mutate(species = "Trip19"),
  Trip22_raw %>% mutate(species = "Trip22")
)

# Rank log2FC within each species
all_long_id <- all_long_id %>%
  group_by(species) %>%
  mutate(rank = rank(-log2FC, ties.method = "min", na.last = "keep")) %>%
  ungroup()

# Save long-format table with IDs
write_csv(all_long_id, "3b_long_log2FC_DAP_rank_protein_id.csv")
