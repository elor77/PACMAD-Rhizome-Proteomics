# Script Title: Visualize HSP Phylogenetic Tree with Expression Data and Clade Labels
# 
# Description: Constructs and visualizes a phylogenetic tree of heat shock proteins (HSPs) 
#              across multiple grass species. Integrates differential protein abundance data 
#              (log2 fold change), overlays expression as a heatmap, highlights significant 
#              changes, and labels orthogroup clades based on specific node positions.
#
# Author: Elad Oren (adapted from Sheng-Kai Hsu)
# Created: 2025-03-20
# Last Modified: 2025-04-20
#
# Input: 
#   - ../1_tables/6a_HSP_protein_og_id_desc.txt - HSP protein-to-orthogroup and description mapping
#   - ../1_tables/3b_long_log2FC_DAP_rank_protein_id.csv - Differential abundance data
#   - ../0_raw_data/Proteome_reference_table.txt - Species-to-proteome reference
#   - ../RAxML/RAxML_bestTree.212_hsp_proteins.tree - Phylogenetic tree from RAxML
#
# Output:
#   - HSP_proteins_for_extraction.txt - Protein IDs and corresponding proteomes for sequence retrieval
#   - HSP_protein_name_mapping.txt - Mapping of full protein IDs to short names
#   - HSP_expression_with_descriptions.txt - Final expression dataset with descriptions
#   - Fig5a_HSP_tree_basic.pdf - Basic phylogenetic tree with expression heatmap
#   - Fig5a_HSP_tree_basic_tips_colored.pdf - Tree with tip colors by orthogroup
#   - Fig5a_HSP_tree_p_nodes.pdf - Tree with node IDs labeled
#   - Fig5a_HSP_tree_bracketed.pdf - Final tree with clade brackets and significance markers
#
# Dependencies: ape, ggtree, ggplot2, ggnewscale, phytools, dplyr, tidyr, readr, stringr
#
# Notes: 
#   - Tripsacum 2019 and 2022 samples are merged into a single "Trip" group.
#   - Significant proteins (padj < 0.05) are marked with asterisks on the tree.
#   - Orthogroup clades are labeled manually based on node positions.
#   - Rooting and tip labels are adjusted for clarity.

# Clear environment and set working directory
rm(list=ls())
setwd("/Users/eo235/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/temp")

# Load required libraries
library(ape)          # Phylogenetic analysis
library(ggplot2)      # Plotting
library(ggtree)       # Tree visualization
library(ggnewscale)   # Multiple scales in plots
library(phytools)     # Phylogenetic tools
library(dplyr)        # Data manipulation
library(tidyr)        # Data tidying
library(readr)        # Reading files
library(stringr)      # String manipulation

# Read and prepare data
data <- read_delim("../1_tables/6a_HSP_protein_og_id_desc.txt", 
                   col_names = c("orthogroup", "protein_id", "description"),
                   delim = "\t")

# Fix duplicate proteins by keeping the longest description
data_deduped <- data %>%
  group_by(protein_id, orthogroup) %>%
  summarize(description = description[which.max(nchar(description))], .groups = "drop")

# Replace original data with deduplicated version
data <- data_deduped

# Read expression data and proteome reference
expression_data <- read_delim("../1_tables/3b_long_log2FC_DAP_rank_protein_id.csv")
proteome_ref <- read_delim("../0_raw_data/Proteome_reference_table.txt")

# Extract all unique orthogroups
orthogroups <- unique(data$orthogroup)

# Match protein IDs with expression data
filtered_expression_data <- expression_data %>%
  filter(protein_id %in% data$protein_id | orthogroup %in% orthogroups)

# Create short names for proteins
expression_with_short_names <- filtered_expression_data %>%
  mutate(
    protein_id_base = gsub("_T\\d+$", "", protein_id),
    short_name = case_when(
      str_detect(protein_id, "^Td") ~ str_extract(protein_id, "^Td\\d+aa\\d+"),
      str_detect(protein_id, "^Ag") ~ str_extract(protein_id, "^Ag\\d+aa\\d+"),
      str_detect(protein_id, "^Pv") ~ str_replace(str_extract(protein_id, "^Pv[^_]+"), 
                                                  "\\.\\d+\\.v\\d+\\.\\d+$", ""),
      str_detect(protein_id, "^Sn") ~ str_extract(protein_id, "^Sn\\d+aa\\d+"),
      str_detect(protein_id, "^Misin") ~ str_replace(str_extract(protein_id, "^Misin[^_]+"), 
                                                     "\\.\\d+\\.v\\d+\\.\\d+$", ""),
      TRUE ~ str_extract(protein_id, "^[^_]+")
    )
  )

# Consolidate Trip19 and Trip22 data and handle short name duplicates
consolidated_expression_data <- expression_with_short_names %>%
  mutate(species = ifelse(species %in% c("Trip19", "Trip22"), "Trip", species)) %>%
  group_by(protein_id, protein_id_base, orthogroup, short_name) %>%
  summarize(
    species = first(species),
    log2FC = mean(log2FC, na.rm = TRUE),
    padj = min(padj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(short_name) %>%
  mutate(has_duplicates = n() > 1) %>%
  ungroup() %>%
  mutate(short_name = if_else(has_duplicates, paste0(short_name, "_", row_number()), short_name)) %>%
  select(-has_duplicates)

# Create consolidated description table by orthogroup
orthogroup_descriptions <- data %>%
  group_by(orthogroup) %>%
  summarize(description = names(table(description))[which.max(table(description))], .groups = "drop")

# Join descriptions with expression data
expression_with_descriptions <- consolidated_expression_data %>%
  left_join(orthogroup_descriptions, by = "orthogroup") %>%
  left_join(select(data, protein_id, specific_description = description), by = "protein_id") %>%
  mutate(description = ifelse(is.na(specific_description), description, specific_description)) %>%
  select(-specific_description)

# Create files for sequence extraction and protein mapping
proteins_for_extraction <- expression_with_descriptions %>%
  select(protein_id, species) %>%
  left_join(proteome_ref, by = "species") %>%
  select(protein_id, proteome)
write_tsv(proteins_for_extraction, "HSP_proteins_for_extraction.txt")

protein_mapping <- expression_with_descriptions %>%
  select(protein_id, short_name) %>%
  distinct()
write_tsv(protein_mapping, "HSP_protein_name_mapping.txt")

write_tsv(expression_with_descriptions, "HSP_expression_with_descriptions.txt")

# Extract the expression data and prepare for visualization
pheno <- as.data.frame(expression_with_descriptions)
pheno_noNA <- na.omit(pheno)
pheno_noNA <- pheno_noNA[pheno_noNA$log2FC != 0, ]
log2fc <- setNames(pheno_noNA$log2FC, pheno_noNA$short_name)

# Load and prepare phylogenetic tree
tre <- read.tree("../RAxML/RAxML_bestTree.212_hsp_proteins.tree")
rooted_tre <- root(tre, outgroup = "Td00001aa001478", resolve.root = TRUE)

# Filter expression data and create subset tree
pheno_filtered <- pheno_noNA %>%
  filter(short_name %in% rooted_tre$tip.label)
subTre <- keep.tip(rooted_tre, pheno_filtered$short_name)

# Function to create tree with heatmap
create_tree_heatmap <- function(tree, data, layout = "rectangular",
                                color_tips = FALSE, min_color = -5.1, max_color = 5.1) {
  # Create base tree
  p <- ggtree(tree, layout = layout, branch.length = "none",
              ladderize = TRUE, size = 0.3) + xlim(0, 40)
  
  # Add tip labels with optional coloring
  if (color_tips) {
    # Match tree tip labels to pheno data
    tip_data <- pheno_noNA[match(tree$tip.label, pheno_noNA$short_name), ]
    n_groups <- length(unique(tip_data$orthogroup))
    
    # Create color palette
    if (n_groups <= 8) {
      colors <- RColorBrewer::brewer.pal(max(8, n_groups), "Set1")
    } else {
      colors <- c(
        RColorBrewer::brewer.pal(8, "Set1"),
        RColorBrewer::brewer.pal(8, "Set2"),
        RColorBrewer::brewer.pal(8, "Set3"),
        RColorBrewer::brewer.pal(8, "Dark2"),
        RColorBrewer::brewer.pal(8, "Paired")
      )
      if (n_groups > length(colors)) {
        more_colors <- viridis::viridis_pal()(n_groups - length(colors))
        colors <- c(colors, more_colors)
      }
    }
    
    orthogroup_factor <- factor(tip_data$orthogroup)
    tip_colors <- colors[as.numeric(orthogroup_factor)]
    p <- p + geom_tiplab(size = 2, color = tip_colors)
  } else {
    p <- p + geom_tiplab(size = 2)
  }
  
  # Add heatmap
  p1 <- gheatmap(p, as.data.frame(data), offset = 2, width = 0.1,
                 colnames = FALSE, color = NA) +
    scale_x_ggtree() +
    scale_y_continuous(expand = c(0, 0.3)) +
    scale_fill_gradient2(
      low = "red",
      mid = "white",
      high = "blue",
      midpoint = 0,
      limits = c(min_color, max_color),
      name = expression(log[2](FC))
    ) +
    theme(legend.position = "left")
  
  return(p1)
}

# Create and save basic heatmap
plot_1 <- create_tree_heatmap(subTre, log2fc, color_tips = FALSE)
print(plot_1)
ggsave("Fig5a_HSP_tree_basic.pdf", plot_1, width = 10, height = 12)

# Create and save heatmap with colored tips
plot_2 <- create_tree_heatmap(subTre, log2fc, color_tips = TRUE)
print(plot_2)
ggsave("Fig5a_HSP_tree_basic_tips_colored.pdf", plot_2, width = 10, height = 12)

# Create tree with node IDs for reference
p_nodes <- ggtree(subTre, layout = "rectangular", branch.length = "none", 
                  ladderize = TRUE, size = 0.3) + 
  xlim(0, 40) +
  geom_tiplab(
    offset = 3, size = 3,
    aes(color = factor(pheno_noNA$orthogroup[match(label, pheno_noNA$short_name)]))
  ) +
  geom_text(aes(label = node), hjust = -0.3, size = 3) +
  #scale_color_viridis_d(option = "D", guide = "none")
  scale_alpha_discrete(guide = "none")

print(p_nodes)
ggsave("Fig5a_HSP_tree_p_nodes.pdf", p_nodes, width = 10, height = 25)

# Define specific nodes for labeling
og_nodes <- c(291, 287, 67, 309, 314, 280, 61, 60, 277, 276, 261, 268, 323, 106, 318, 100, 251, 38, 256, 247, 245, 356, 388, 350, 346, 341, 401, 395, 333, 412, 195, 202, 224, 231, 238, 237, 19)
Ntip <- length(subTre$tip.label)
first_tips <- sapply(og_nodes, function(nd) {
  if (nd <= Ntip) {
    subTre$tip.label[nd]
  } else {
    extract.clade(subTre, nd)$tip.label[1]
  }
})
node_labels <- pheno_noNA$orthogroup[match(first_tips, pheno_noNA$short_name)]

# Add clade labels and significance markers
plot_bracketed <- plot_1
for(i in seq_along(og_nodes)) {
  plot_bracketed <- plot_bracketed +
    geom_cladelabel(
      node = og_nodes[i],
      label = node_labels[i],
      align = TRUE,
      angle = 45,
      offset = 5,
      barsize = 0.5,
      fontsize = 3,
      offset.text = 0.5
    )
}

# Function to determine if a tip is significant
is_significant <- function(tip_label) {
  padj_value <- pheno_noNA$padj[match(tip_label, pheno_noNA$short_name)]
  return(!is.na(padj_value) && padj_value < 0.05)
}

# Add asterisks for significant tips
plot_bracketed <- plot_bracketed +
  geom_tiplab(
    aes(label = ifelse(sapply(label, is_significant), "*", "")),
    offset = 2.5,
    color = "black",
    size = 3,
    fontface = "bold"
  ) +
  theme(plot.margin = margin(5, 80, 5, 10, "pt"))

print(plot_bracketed)

# Save final plot
ggsave("Fig5a_HSP_tree_bracketed.pdf", plot_bracketed, width = 14, height = 15)
