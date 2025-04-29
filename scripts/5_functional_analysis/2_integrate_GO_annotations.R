# Script Title: GO Annotation, Enrichment, and Visualization for Differentially Abundant Proteins
# 
# Description: Integrates PANNZER2 functional annotations with differentially abundant 
#              protein data, performs GO biological process enrichment analysis using topGO, 
#              and visualizes the most significantly enriched GO terms for up- and 
#              down-regulated proteins across species.
#
# Author: Elad Oren
# Created: 2025-01-01
# Last Modified: 2025-04-27
#
# Input: 
#   - 4a_all_DAPs_log2FC.csv - Significant DAPs with log2 fold changes and species info
#   - PANNZER2 annotation files for each species (Ag_anno.out, Mg_anno.out, Pv_anno.out, Sn_anno.out, Trip_anno.out)
#
# Output:
#   - 4b_all_DAPs_with_DE_BP.csv - Merged DAPs with protein descriptions and GO annotations
#   - 4c_GO_enrichment_UP.csv - GO enrichment results for up-regulated proteins
#   - 4c_GO_enrichment_DOWN.csv - GO enrichment results for down-regulated proteins
#   - (Plots generated interactively for visualization; can be saved manually)
#
# Dependencies: tidyverse (dplyr, readr, ggplot2, tidyr), cowplot, topGO
#
# Notes: GO enrichment analysis is performed separately for up- and down-regulated proteins 
#        using Fisher's exact test with the "elim" algorithm in topGO. Bubble plots show 
#        -log10(p-value) versus GO biological process term with color-coded gene ratios 
#        and point size indicating gene counts.

library(tidyverse)
library(ggplot2)
library(cowplot)
library(topGO)


# Clear Environment
rm(list = ls())

setwd("~/OneDrive/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/1_tables/")


# Load your master DAP table
daps <- read_csv("4a_all_DAPs_log2FC.csv")

# Define species list
species_list <- c("Ag", "Mg", "Pv", "Sn", "Trip")

# Path to annotation files
anno_path <- "/Users/eo235/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/PANNZER_annotations/input/"

# Read and combine all annotation files
annotations <- map_dfr(species_list, function(sp) {
  file_path <- file.path(anno_path, paste0(sp, "_anno.out"))
  read_tsv(file_path, col_types = cols()) %>%
    mutate(species = sp)
})

# Keep only DE (description) and BP (Biological Process) annotations
annotations_filtered <- annotations %>%
  filter(type == "DE" | str_detect(type, "BP"))

DE_annotations <- annotations_filtered %>%
  filter(type == "DE") %>%
  select(protein_id = qpid, description = desc)

BP_annotations <- annotations_filtered %>%
  filter(str_detect(type, "BP")) %>%
  select(protein_id = qpid, go_id = id, go_term = desc, PPV)

# Merge DE annotations
daps_anno <- daps %>%
  left_join(DE_annotations, by = "protein_id")

# Merge BP annotations
daps_full <- daps_anno %>%
  left_join(BP_annotations, by = "protein_id")

# Save merged table
write_csv(daps_full, "4b_all_DAPs_with_DE_BP.csv")

######TopGO analysis######

# Read your merged DAP + annotation file
daps_full <- read_csv("4b_all_DAPs_with_DE_BP.csv")

# Prepare gene-to-GO mapping
gene2GO <- daps_full %>%
  filter(!is.na(go_id)) %>%
  distinct(protein_id, go_id) %>%
  mutate(go_id = paste0("GO:", go_id)) %>%
  group_by(protein_id) %>%
  summarise(go_terms = list(go_id), .groups = "drop") %>%
  deframe()

# Background gene universe: all unique protein_ids in your dataset
gene_universe <- unique(daps_full$protein_id)

# Subset DAPs into Up and Down
up_genes <- daps_full %>%
  filter(DAP == "Up") %>%
  pull(protein_id) %>%
  unique()

down_genes <- daps_full %>%
  filter(DAP == "Down") %>%
  pull(protein_id) %>%
  unique()

# Create gene lists for topGO (all genes = 0 or 1)
gene_list_up <- factor(as.integer(gene_universe %in% up_genes))
names(gene_list_up) <- gene_universe

gene_list_down <- factor(as.integer(gene_universe %in% down_genes))
names(gene_list_down) <- gene_universe

# Create topGOdata objects
GOdata_up <- new("topGOdata",
                 ontology = "BP",
                 allGenes = gene_list_up,
                 annot = annFUN.gene2GO,
                 gene2GO = gene2GO,
                 nodeSize = 10)

GOdata_down <- new("topGOdata",
                   ontology = "BP",
                   allGenes = gene_list_down,
                   annot = annFUN.gene2GO,
                   gene2GO = gene2GO,
                   nodeSize = 10)

# Run enrichment tests
result_up <- runTest(GOdata_up, algorithm = "elim", statistic = "Fisher")
result_down <- runTest(GOdata_down, algorithm = "elim", statistic = "Fisher")

# Create result tables
table_up <- GenTable(GOdata_up,
                     elimFisher = result_up,
                     topNodes = 50)

table_down <- GenTable(GOdata_down,
                       elimFisher = result_down,
                       topNodes = 50)

# Save results to CSV
write_csv(table_up, "4c_GO_enrichment_UP.csv")
write_csv(table_down, "4c_GO_enrichment_DOWN.csv")


####Visualization####
# Process UPregulated
plot_up <- table_up %>%
  mutate(
    elimFisher = as.numeric(elimFisher),  # important: topGO outputs it as character
    GeneRatio = Significant / Annotated,
    logP = -log10(elimFisher),
    Regulation = "Upregulated"
  ) %>%
  filter(!is.na(elimFisher), Annotated >= 5) %>%
  slice_max(order_by = logP, n = 15, with_ties = FALSE) %>%
  arrange(logP) %>%
  mutate(Term = factor(Term, levels = Term))

# Process DOWNregulated
plot_down <- table_down %>%
  mutate(
    elimFisher = as.numeric(elimFisher),
    GeneRatio = Significant / Annotated,
    logP = -log10(elimFisher),
    Regulation = "Downregulated"
  ) %>%
  filter(!is.na(elimFisher), Annotated >= 5) %>%
  slice_max(order_by = logP, n = 15, with_ties = FALSE) %>%
  arrange(logP) %>%
  mutate(Term = factor(Term, levels = Term))

#Plot UP
p_up <- ggplot(plot_up, aes(x = logP, y = Term, size = Significant, color = GeneRatio)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() +
  labs(
    title = "Upregulated",
    x = "-log10(p-value)",
    y = NULL,
    color = "Gene Ratio",
    size = "Gene Count"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 14)
  )

#Plot Down
p_down <- ggplot(plot_down, aes(x = logP, y = Term, size = Significant, color = GeneRatio)) +
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() +
  labs(
    title = "Downregulated",
    x = "-log10(p-value)",
    y = NULL,
    color = "Gene Ratio",
    size = "Gene Count"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 14)
  )

plot_grid(p_down, p_up, ncol = 2, align = "h")


