# Script Title: Visualize Consolidated Biological Processes in Up-Regulated DAPs
# 
# Description: Joins up-regulated differentially abundant proteins (DAPs) to manually 
#              assigned consolidated biological process (BP) categories, generates 
#              a mosaic-style stacked proportional bar plot showing BP category 
#              distributions across six grass species, and performs a Chi-squared 
#              test to evaluate differences among species.
#
# Author: Elad Oren
# Created: 2025-01-01
# Last Modified: 2025-04-27
#
# Input: 
#   - 4f_daps_final_collapsed_complemented_rescued_BP_UP_only.csv - Final DAP dataset with BP terms
#   - 4g_daps_UP_consolidated_BPs.csv - Mapping of GO terms to consolidated BP categories
#   - BP_Category_Colors.csv - Color mapping for BP categories
#
# Output:
#   - ../Figures/FigS3_DAPs_mosaic_plot_consolidated_BP_final.pdf - Mosaic plot of BP categories by species
#
# Dependencies: tidyverse (dplyr, ggplot2, readr, stringr), scales, RColorBrewer
#
# Notes: Species and BP categories are manually ordered. The Chi-squared test is 
#        performed to test for independence between species and BP category distributions.
#        Optionally, "Other" BP terms can be excluded from the test.

# Load packages
library(tidyverse)
library(scales)
library(RColorBrewer)

# Clear Environment
rm(list = ls())

setwd("~/OneDrive/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/1_tables/")


# Step 0: Load rescued DAPs file and mapping table
daps_final_rescued <- read_csv("4f_daps_final_collapsed_complemented_rescued_BP_UP_only.csv")
bp_mapping_table <- read_csv("4g_daps_UP_consolidated_BPs.csv")

# Step 1: Join mapping table to rescued data and prepare the data
plot_data <- daps_final_rescued %>%
  filter(!is.na(Term)) %>%
  # join mapping based on matching the Term with Biological_Process
  left_join(bp_mapping_table, by = c("Term" = "Biological_Process")) %>%
  # replace NA Assigned_Category with "Other"
  mutate(Assigned_Category = coalesce(Assigned_Category, "Other")) %>%
  count(species, Assigned_Category)

# Step 3: Define species order
species_order <- c("Trip19", "Trip22", "Ag", "Mg", "Pv", "Sn")

# Step 4: Sort categories descending by total count, "Other" at the bottom
BP_counts <- plot_data %>%
  group_by(Assigned_Category) %>%
  summarise(total_n = sum(n), .groups = "drop") %>%
  mutate(is_other = if_else(tolower(Assigned_Category) == "other", 1, 0)) %>%
  arrange(is_other, desc(total_n))

BP_order <- BP_counts$Assigned_Category

# Load the color mapping table
bp_colors_df <- read.csv("BP_Category_Colors.csv", stringsAsFactors = FALSE)

# Create a named character vector for ggplot2
bp_colors <- setNames(bp_colors_df$Color, bp_colors_df$Name)

# Step 6: Calculate total N
N_total <- sum(plot_data$n)

# Step 7: Plot
ggplot(plot_data, aes(x = factor(species, levels = species_order),
                      y = n,
                      fill = factor(Assigned_Category, levels = BP_order))) +
  geom_bar(stat = "identity", position = "fill", color = "black") +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.1))) +
  scale_fill_manual(values = bp_colors, breaks = BP_order) +
  annotate("text", x = 1, y = 1.05, label = paste0("N = ", N_total), hjust = 0, size = 3.5) +
  labs(
    x = "Species",
    y = "Cumulative relative abundance",
    fill = "Consolidated Biological Process"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

ggsave("../Figures/FigS3_DAPs_mosaic_plot_consolidated_BP_final.pdf", width = 10, height = 6, dpi = 300)

# Step 8: Chi-square test for biological processes
contingency_table <- xtabs(n ~ species + Assigned_Category, data = plot_data)

# Optionally drop "other" before chi-square test
#contingency_table <- contingency_table[, colnames(contingency_table) != "Other"]

chisq_test_result <- chisq.test(contingency_table)
#chisq_test_result <- chisq.test(contingency_table, simulate.p.value = TRUE, B = 10000)

print(chisq_test_result)
