# Script Title: Visualize Cold-Related Functional Categories in Top 50 DAPs
# 
# Description: Counts and visualizes cold-related functional categories among the 
#              manually curated top 50 up-regulated differentially abundant proteins (DAPs) 
#              across six grass species. Generates a mosaic-style stacked bar plot and 
#              performs a Chi-squared test to evaluate differences in category distributions.
#
# Author: Elad Oren
# Created: 2025-01-01
# Last Modified: 2025-04-27
#
# Input: 
#   - 5b_top50_proteins_manual_curation.csv - Manually curated top 50 up-regulated proteins dataset
#
# Output:
#   - ../Figures/Fig4a_Top50_mosaic_plot_cold-categories.pdf - Mosaic plot of cold-related functional categories
#
# Dependencies: tidyverse (dplyr, readr, ggplot2), scales, viridis
#
# Notes: Cold-related categories are manually defined and colored. The Chi-squared test assesses 
#        whether the distribution of functional categories differs significantly among species. 
#        The "Other" category can be optionally excluded before statistical testing.


library(readr)
library(dplyr)
library(ggplot2)
library(viridis)

rm(list = ls())
setwd("~/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/1_tables")


top50_proteins <- read_csv("5b_top50_proteins_manual_curation.csv")


###Now for cold-related categories

# Step 1: Prepare cold data
cold_plot_data <- top50_proteins %>%
  count(species, cold_related_category)

# Step 2: Species ordering
species_order <- c("Trip19", "Trip22", "Ag", "Mg", "Pv", "Sn")

# Summarize counts
cold_counts <- cold_plot_data %>%
  group_by(cold_related_category) %>%
  summarise(total_n = sum(n), .groups = "drop")

# Sort ONLY non-"Other" categories by descending abundance
cold_order <- cold_counts %>%
  filter(tolower(cold_related_category) != "other") %>%
  arrange(desc(total_n)) %>%
  pull(cold_related_category)

# Finally, add "Other" at the end
cold_order <- c(cold_order, "Other")

# Step 4: Set colors manually
colors_manual_cold <- c(
  "Other" = "white",
  "Metabolism / osmoregulation" = "#E69F00",
  "Protein folding / stability" = "#56B4E9",
  "Protein aggregation / membrane stability" = "#92CDDD",
  "Detoxification / ROS scavenging" = "#F0E442",
  "Cold signal transduction" = "#009E73",
  "Cytoskeleton stability" = "#D55E00",
  "Antifreeze" = "#999999",
  "Lipid metabolism" = "#CC79A7"
)

# Step 5: Plot
ggplot(cold_plot_data, aes(x = factor(species, levels = species_order),
                           y = n,
                           fill = factor(cold_related_category, levels = cold_order))) +
  geom_bar(stat = "identity", position = "fill", color = "black", size = 0.25) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(
    values = colors_manual_cold,
    breaks = cold_order
  ) +
  labs(
    x = "Species",
    y = "Cumulative relative abundance",
    fill = "Cold related category"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )
# Save the plot
ggsave("../Figures/Fig4a_Top50_mosaic_plot_cold-categories.pdf", width = 8, height = 5, dpi = 300)


# Step 6: Chi-square test for cold-related categories
# Correct contingency table using the real counts
contingency_table <- xtabs(n ~ species + cold_related_category, data = cold_plot_data)

# Drop the "other" column
#contingency_table <- contingency_table[, colnames(contingency_table) != "Other"]

# Now run the chi-square test
chisq_test_result <- chisq.test(contingency_table)

# View results
chisq_test_result
