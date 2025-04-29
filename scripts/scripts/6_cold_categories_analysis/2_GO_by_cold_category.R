# Script Title: Visualize Biological Process Contributions Within Cold-Related Categories
# 
# Description: Analyzes the top 50 up-regulated differentially abundant proteins (DAPs) 
#              assigned to cold-related functional categories. Maps biological process (BP) 
#              annotations to consolidated BP categories and generates a stacked bar plot 
#              showing the distribution of biological processes within each cold-related category.
#
# Author: Elad Oren
# Created: 2025-01-01
# Last Modified: 2025-04-27
#
# Input: 
#   - 5b_top50_proteins_manual_curation.csv - Manually curated top 50 up-regulated proteins
#   - 4g_daps_UP_consolidated_BPs.csv - Mapping of GO terms to consolidated BP categories
#   - BP_Category_Colors.csv - Color mapping for BP categories
#
# Output:
#   - ../Figures/FigS1_cold_category_bp_counts.pdf - Stacked bar plot of BP category distributions by cold-related category
#
# Dependencies: tidyverse (dplyr, readr, ggplot2), scales
#
# Notes: Proteins with missing or unmapped BP terms are assigned to the "Other" category. 
#        BP categories and cold-related categories are manually ordered for clarity in visualization. 
#        Plot colors are defined according to a pre-specified palette.


# Load necessary libraries
library(ggplot2)
library(dplyr)

rm(list = ls())
# Read the data files
# Read the BP mappings file with proper encoding
bp_mappings <- read.csv("4g_daps_UP_consolidated_BPs.csv", stringsAsFactors = FALSE, check.names = FALSE)

# Read the protein data file with proper encoding
protein_data <- read.csv("5b_top50_proteins_manual_curation.csv", stringsAsFactors = FALSE, check.names = FALSE)

# Create a mapping from BP terms to consolidated terms
bp_to_consolidated <- setNames(bp_mappings$Assigned_Category, bp_mappings$Biological_Process)

# Count proteins by cold-related category and biological process
# FIXED: Don't filter out NA BPs, instead map them to "Other"
filtered_data <- protein_data %>%
  filter(!is.na(cold_related_category) & 
           cold_related_category != "Other")

# Debug: Check how many rows were kept 
print("Rows after filtering 'Other' cold categories:")
print(nrow(filtered_data))

# Debug: Check how many NA BPs we have (these will be mapped to "Other")
print("Number of NA BPs in filtered data:")
print(sum(is.na(filtered_data$BP)))

# Add the consolidated BP column - correctly handle NA values
consolidated_data <- filtered_data %>%
  mutate(consolidated_bp = case_when(
    is.na(BP) ~ "Other",
    BP %in% names(bp_to_consolidated) ~ bp_to_consolidated[BP],
    TRUE ~ "Other"
  ))

# Calculate counts
bp_counts <- consolidated_data %>%
  group_by(cold_related_category, consolidated_bp) %>%
  summarize(count = n(), .groups = "drop")

# Calculate total proteins for subtitle
total_proteins <- sum(bp_counts$count)

# Order the cold categories by total protein count
category_order <- bp_counts %>%
  group_by(cold_related_category) %>%
  summarize(total = sum(count)) %>%
  arrange(desc(total)) %>%
  pull(cold_related_category)


bp_counts$cold_related_category <- factor(bp_counts$cold_related_category, 
                                          levels = category_order)

# Step 1: Order BPs by total count, excluding "Other"
bp_order <- bp_counts %>%
  filter(consolidated_bp != "Other") %>%
  group_by(consolidated_bp) %>%
  summarise(total_n = sum(count), .groups = "drop") %>%
  arrange(desc(total_n)) %>%
  pull(consolidated_bp)

# Step 2: Add "Other" at the end
if ("Other" %in% bp_counts$consolidated_bp) {
  bp_order <- c(bp_order, "Other")
}

# Step 3: Re-factor for correct fill stacking and legend order
bp_counts$consolidated_bp <- factor(
  bp_counts$consolidated_bp,
  levels = bp_order
)

# Load the color mapping table
bp_colors_df <- read.csv("BP_Category_Colors.csv", stringsAsFactors = FALSE)

# Create a named character vector for ggplot2
bp_colors <- setNames(bp_colors_df$Color, bp_colors_df$Name)

# Create the stacked bar chart
p <- ggplot(bp_counts, aes(x = cold_related_category, y = count, fill = consolidated_bp)) +
  geom_bar(stat = "identity", position = "stack", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = bp_colors) +
  labs(
    x = "Cold-related category",
    y = "Protein count",
    fill = "Biological process"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

# Print the plot to screen
print(p)

# Save the plots
#ggsave("cold_category_bp_counts.png", p, width = 14, height = 8, dpi = 300)
ggsave("../Figures/FigS1_cold_category_bp_counts.pdf", p, width = 10, height = 8)
 