# Script Title: Prepare Protein ID Lists for Differentially Abundant Proteins
# 
# Description: Filters significant differentially abundant proteins (DAPs) from the 
#              long-format dataset, ranks them by log2 fold change, and generates 
#              species-specific protein ID lists for downstream sequence extraction 
#              and functional annotation. Tripsacum 2019 and 2022 datasets are merged 
#              into a single group.
#
# Author: Elad Oren
# Created: 2025-01-01
# Last Modified: 2025-04-27
#
# Input: 
#   - 3b_long_log2FC_DAP_rank_protein_id.csv - Long-format dataset with log2FC and DAP classifications
#
# Output:
#   - 4_all_DAPs_log2FC.csv - Filtered dataset of Up/Down regulated DAPs with ranks
#   - protein_ids_Tripsacum.txt - Protein IDs for Tripsacum (combined Trip19 and Trip22)
#   - protein_ids_Ag.txt - Protein IDs for Andropogon gerardii
#   - protein_ids_Mg.txt - Protein IDs for Miscanthus giganteus
#   - protein_ids_Pv.txt - Protein IDs for Panicum virgatum
#   - protein_ids_Sn.txt - Protein IDs for Sorghastrum nutans
#
# Dependencies: dplyr, readr, tidyverse
#
# Notes: Protein ID files are prepared separately for each species. Tripsacum IDs are 
#        combined from 2019 and 2022 datasets to form a single list. Outputs are used 
#        for subsequent protein sequence extraction and GO annotation workflows.

library(tidyverse)

# Clear Environment

rm(list = ls())

setwd("~/OneDrive/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/1_tables/")


# Load DAP Dataset

long_log2FC <- read_csv("3b_long_log2FC_DAP_rank_protein_id.csv")


# Extract Significant DAPs

all_DAPs <- long_log2FC %>%
  filter(DAP %in% c("Up", "Down")) %>%
  group_by(species) %>%
  arrange(desc(log2FC), .by_group = TRUE) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  select(rank, species, protein_id, orthogroup, log2FC, DAP)


# Save all DAPs to CSV
write_csv(all_DAPs, "4_all_DAPs_log2FC.csv")


# Summarize DAP counts per species (optional)

all_DAPs_summary <- all_DAPs %>%
  group_by(species, DAP) %>%
  summarise(count = n(), .groups = "drop")

# # Save summary if needed
# write_csv(all_DAPs_summary, "all_DAPs_summary.csv")


# Handle Trip19 and Trip22 Together


# Combine Trip19 and Trip22 with unique protein IDs
trip_DAPs <- all_DAPs %>%
  filter(species %in% c("Trip19", "Trip22")) %>%
  distinct(protein_id, .keep_all = TRUE)

# Extract other species separately
other_DAPs <- all_DAPs %>%
  filter(!species %in% c("Trip19", "Trip22"))


# Write Output Files for annotation


# Write combined Tripsacum protein IDs
write.table(trip_DAPs$protein_id,
            file = "protein_ids_Tripsacum.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

# Write other species separately
other_DAPs %>%
  group_by(species) %>%
  group_walk(~ write.table(.x$protein_id,
                           file = paste0("protein_ids_", .y$species, ".txt"),
                           quote = FALSE, row.names = FALSE, col.names = FALSE))


# Final Message

cat("Protein ID files prepared for: Tripsacum, Ag, Mg, Pv, Sn.\n")
cat("You can now proceed to sequence extraction and annotation.\n")