# Script Title: Collapse and Complement GO Biological Process Annotations for Up-Regulated DAPs
# 
# Description: Assigns the best enriched GO biological process (BP) term to each 
#              up-regulated differentially abundant protein (DAP). Complements missing 
#              annotations based on orthogroup relationships and rescues remaining 
#              missing proteins using available GO annotations from the full dataset.
#
# Author: Elad Oren
# Created: 2025-01-01
# Last Modified: 2025-04-27
#
# Input: 
#   - 4b_all_DAPs_with_DE_BP.csv - Full DAP + functional annotation dataset
#   - 4c_GO_enrichment_UP.csv - GO enrichment results for up-regulated proteins
#
# Output:
#   - 4d_collapsed_BP_UP_only.csv - Best enriched BP term per up-regulated protein (collapsed version)
#   - 4e_daps_collapsed_BP_UP_only.csv - Collapsed dataset merged with DAP information
#   - 4f_daps_final_collapsed_complemented_rescued_BP_UP_only.csv - Final complemented and rescued BP assignments
#
# Dependencies: tidyverse (dplyr, readr, stringr)
#
# Notes: GO assignment prioritizes enriched BP terms with the lowest elimFisher p-value. 
#        Missing annotations are first complemented based on orthogroup relationships, 
#        then rescued by direct protein ID lookup if available. Final dataset is used 
#        for biological process categorization and downstream visualization (e.g., mosaic plots).

library(tidyverse)

# Clear Environment
rm(list = ls())

# Set working directory
setwd("~/OneDrive/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/1_tables/")

# Load full DAP + BP annotation table
daps_full <- read_csv("4b_all_DAPs_with_DE_BP.csv")

# Filter for Up-regulated proteins only
daps_up <- daps_full %>%
  filter(DAP == "Up")

# Load GO enrichment for UP only
table_up <- read_csv("4c_GO_enrichment_UP.csv")

# Prepare list of enriched BP GO terms (elimFisher < 0.05)
enriched_BP <- table_up %>%
  mutate(elimFisher = as.numeric(elimFisher)) %>%
  filter(elimFisher < 0.05) %>%
  select(GO_ID = GO.ID, Term, elimFisher)

# Format GO ID properly
daps_up <- daps_up %>%
  mutate(go_id = if_else(str_detect(as.character(go_id), "^GO:"), 
                         as.character(go_id), 
                         paste0("GO:", str_pad(as.character(go_id), 7, pad = "0"))))

# Collapse: select best enriched BP for each protein
collapsed_bp <- daps_up %>%
  filter(!is.na(go_id)) %>%
  left_join(enriched_BP, by = c("go_id" = "GO_ID")) %>%
  mutate(elimFisher = as.numeric(elimFisher)) %>%
  filter(!is.na(elimFisher)) %>%
  group_by(protein_id, go_id, Term, PPV) %>%
  summarise(min_elimFisher = min(elimFisher, na.rm = TRUE), .groups = "drop") %>%
  group_by(protein_id) %>%
  slice_min(order_by = min_elimFisher, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(protein_id, go_id, Term, PPV)

# Save collapsed BP table
write_csv(collapsed_bp, "4d_collapsed_BP_UP_only.csv")

# Merge collapsed BP back onto DAPs
daps_collapsed <- daps_up %>%
  distinct(protein_id, species, orthogroup, rank, log2FC, DAP, description) %>%
  left_join(collapsed_bp, by = "protein_id")

write_csv(daps_collapsed, "4e_daps_collapsed_BP_UP_only.csv")

cat("=== Summary after initial collapse ===\n")
daps_collapsed %>%
  summarise(
    total_proteins = n(),
    proteins_with_GO = sum(!is.na(go_id)),
    proteins_without_GO = sum(is.na(go_id))
  ) %>%
  print()



# Impute missing GO terms based on orthogroup

# 1. Separate proteins with and without BP
daps_with_bp <- daps_collapsed %>% filter(!is.na(go_id))
daps_without_bp <- daps_collapsed %>% filter(is.na(go_id))

# 2. Best BP per Orthogroup (highest PPV)
best_bp_per_OG <- daps_with_bp %>%
  group_by(orthogroup) %>%
  arrange(desc(PPV)) %>%
  slice(1) %>%
  ungroup() %>%
  select(orthogroup, go_id, Term, PPV)

# 3. Complement missing BPs
daps_complemented <- daps_without_bp %>%
  filter(!is.na(orthogroup)) %>%
  left_join(best_bp_per_OG, by = "orthogroup") %>%
  mutate(
    go_id = go_id.y,
    Term = Term.y,
    PPV = PPV.y,
    complemented = if_else(!is.na(go_id), TRUE, FALSE)
  ) %>%
  select(-go_id.y, -Term.y, -PPV.y, -go_id.x, -Term.x, -PPV.x)

# 4. Handle proteins without orthogroup
daps_without_bp_no_og <- daps_without_bp %>%
  filter(is.na(orthogroup)) %>%
  mutate(complemented = FALSE)

# 5. Recombine everything
daps_final <- bind_rows(
  daps_with_bp %>% mutate(complemented = FALSE),
  daps_complemented,
  daps_without_bp_no_og
)

# # Save final table
# write_csv(daps_final, "daps_final_collapsed_complemented_BP_UP_only.csv")

cat("\n=== Summary after orthogroup complementing ===\n")
daps_final %>%
  summarise(
    total_proteins = n(),
    proteins_with_GO = sum(!is.na(go_id)),
    proteins_without_GO = sum(is.na(go_id))
  ) %>%
  print()

###Fill from table
# Find proteins still missing GO info
still_missing <- daps_final %>%
  filter(is.na(go_id) | is.na(Term)) %>%
  select(protein_id)

# Lookup in full daps_full table
rescue_lookup <- daps_full %>%
  filter(protein_id %in% still_missing$protein_id) %>%
  filter(!is.na(go_id)) %>%
  mutate(go_id = if_else(str_detect(as.character(go_id), "^GO:"), 
                         as.character(go_id), 
                         paste0("GO:", str_pad(as.character(go_id), 7, pad = "0")))) %>%
  group_by(protein_id) %>%
  arrange(desc(PPV)) %>%
  slice(1) %>%
  ungroup() %>%
  select(protein_id, go_id, go_term, PPV)  # <-- pay attention: go_term, not Term

# Merge rescued annotations
daps_final_rescued <- daps_final %>%
  left_join(rescue_lookup, by = "protein_id", suffix = c("", "_rescued")) %>%
  mutate(
    go_id = if_else(is.na(go_id), go_id_rescued, go_id),
    Term = if_else(is.na(Term), go_term, Term),  # use go_term here!
    PPV = if_else(is.na(PPV), PPV_rescued, PPV),
    rescued_by_protein_id = if_else(!is.na(go_id_rescued), TRUE, FALSE)
  ) %>%
  select(-go_id_rescued, -go_term, -PPV_rescued)  # drop correctly

# Save final rescued version
write_csv(daps_final_rescued, "4f_daps_final_collapsed_complemented_rescued_BP_UP_only.csv")


cat("\n=== Final summary after protein_id rescue ===\n")
daps_final_rescued %>%
  summarise(
    total_proteins = n(),
    proteins_with_GO = sum(!is.na(go_id)),
    proteins_without_GO = sum(is.na(go_id))
  ) %>%
  print()
