# Script Title: PCA Analysis of Multispecies and Tripsacum 2022 Proteomics
# 
# Description: Performs Principal Component Analysis (PCA) on multispecies 
#              rhizome proteomics data and separately on Tripsacum 2022 batch-normalized 
#              proteomics data. Cleans sample metadata, computes PCA, visualizes results, 
#              and outputs PCA coordinate tables for downstream analyses.
#
# Author: Elad Oren
# Created: 2025-01-01
# Last Modified: 2025-04-27
#
# Input: 
#   - 1_rhizomes_species_filtered.csv - Multispecies proteomics dataset
#   - 0_Trip2022_clean-data_removed-NA.csv - Tripsacum 2022 batch-normalized dataset
#
# Output:
#   - 1A_PCA_coordinates_species.csv - PCA coordinates for multispecies data
#   - 1B_PCA_coordinates_trip2022.csv - PCA coordinates for Tripsacum 2022 data
#
# Dependencies: dplyr, readr, ggfortify, ggplot2
#
# Notes: Prepares sample metadata for each dataset to align with PCA structure. 
#        Visualizes PCA with `autoplot` from ggfortify package. PCA is performed 
#        after scaling (centered and scaled variables). Final PCA coordinate tables 
#        are saved for use in further analyses and plotting.

# Clean everything first
rm(list = ls())

# Load libraries
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(purrr)
library(stringr)

# Step 1: Set working directory
setwd("~/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/1_tables/")

# Step 2: Define file paths
files <- c(
  "0_Ag_log-normalized-abundances_long-format.csv",
  "0_Mg_3_log-normalized-abundances_long-format.csv",
  "0_Pv-DD_log-normalized-abundances_long-format.csv",
  "0_Sn_log-normalized-abundances_long-format.csv",
  "0_Trip2019_log-normalized-abundances_long-format.csv",
  "0_Trip2022_log-normalized-abundances_long-format.csv"
)

# Step 3: Read each file and assign species based on filename
list_of_dfs <- map(files, function(file) {
  species <- case_when(
    str_detect(file, "Ag") ~ "Ag",
    str_detect(file, "Mg") ~ "Mg",
    str_detect(file, "Pv-DD") ~ "Pv",
    str_detect(file, "Sn") ~ "Sn",
    str_detect(file, "Trip2019") ~ "Trip19",
    str_detect(file, "Trip2022") ~ "Trip22",
    TRUE ~ NA_character_
  )
  
  read_csv(file) %>%
    mutate(species = species)
})

# Step 4: Combine everything
all_data <- bind_rows(list_of_dfs)

# Step 5: Read orthogroup mapping
orthogroups <- read_delim("0_Orthogroups_Ag_Mg_Pv_Sc_Sn_Td-v1.0.txt")

# Step 6: Add orthogroup to all_data
all_data <- all_data %>%
  left_join(
    orthogroups %>% select(protein_id, orthogroup),
    by = "protein_id"
  )

# Step 7: Prepare PCA input (pivot wider)
pca_data <- all_data %>%
  filter(!is.na(orthogroup)) %>%
  group_by(sample, species, season, orthogroup) %>%
  summarise(log_abundance = mean(log_abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = orthogroup,
    values_from = log_abundance
  )

threshold <- 0.5

# Step 1: Pivot longer
pca_long <- pca_data %>%
  pivot_longer(cols = -c(sample, species, season), names_to = "OG", values_to = "log_abundance")

# Step 2: Calculate per species per OG fraction detected
og_fraction <- pca_long %>%
  group_by(species, OG) %>%
  summarise(
    n_present = sum(!is.na(log_abundance)),
    n_total = n(), 
    fraction_present = n_present / n_total,
    .groups = "drop"
  )

# Step 3: For each OG, check if ALL species pass threshold
ogs_to_keep <- og_fraction %>%
  group_by(OG) %>%
  summarise(all_pass = all(fraction_present > threshold), .groups = "drop") %>%
  filter(all_pass) %>%
  pull(OG)


# Step 10: Filter PCA data
pca_filtered <- pca_data %>%
  select(sample, species, season, all_of(ogs_to_keep))

# Step 11: Prepare input for PCA
pca_input <- pca_filtered %>% select(-sample, -species, -season)

# Step 12: Remove columns with missing values
pca_input_clean <- pca_input %>%
  select(where(~ all(!is.na(.))))

# Step 13: Run PCA
pca_result_all <- PCA(pca_input_clean, graph = FALSE)

# Step 14: Prepare sample_season for plotting
pca_filtered_season <- pca_filtered %>%
  mutate(sample_season = paste0(sample, "_", season))

pca_coords <- as.data.frame(pca_result_all$ind$coord) %>%
  mutate(sample_season = pca_filtered_season$sample_season,
         species = pca_filtered_season$species)

# Step 15: Prepare metadata
pca_metadata <- all_data %>%
  mutate(sample_season = paste0(sample, "_", season)) %>%
  select(sample_season, species, season) %>%
  distinct(sample_season, species, .keep_all = TRUE)

# Step 16: Safe join
pca_plot_data <- left_join(pca_coords, pca_metadata, by = c("sample_season", "species"))

# Step 17: Plot PCA
pca_all <- ggplot(pca_plot_data, aes(x = Dim.1, y = Dim.2, color = species, shape = season)) +
  geom_point(size = 3) +
  labs(
    title = paste0("PCA of Log Protein Abundance (Filtered OGs: ", length(ogs_to_keep), ")"),
    x = paste0("PC1 (", round(pca_result_all$eig[1, 2], 1), "% variance)"),
    y = paste0("PC2 (", round(pca_result_all$eig[2, 2], 1), "% variance)")
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom"
  )


ggsave("../Figures/FigS2A_PCA_all.pdf", plot = pca_all, width = 8, height = 8, units = "in", dpi = 300)


print(pca_all)

# Step 1: Filter only Trip19 and Trip22
trip_data <- all_data %>%
  filter(species %in% c("Trip19", "Trip22"))

# Step 2: Assign tmt_batch = 5 to Trip19 samples
trip_data <- trip_data %>%
  mutate(tmt_batch = ifelse(species == "Trip19", 5, tmt_batch))


# Read Hybrid class
hybrid_class <- read_delim("../0_raw_data/Tripsacum_samples_Tf-vs-Td.txt")

# Remove replicate from hybrid_class table
hybrid_class <- hybrid_class %>%
  mutate(sample_core = str_remove(sample_name, "_[0-9]+$"))

# Now deduplicate (keep only one Hybrid_class per core sample)
hybrid_class <- hybrid_class %>%
  distinct(sample_core, Hybrid_class)

# Prepare trip_data
trip_data <- trip_data %>%
  mutate(
    sample_name = str_replace(sample, "-", "_"),        # Replace "-" with "_"
    sample_core = str_remove(sample_name, "_[0-9]+$")    # Remove _1, _2, _3
  )

# Join by sample_core
trip_data <- trip_data %>%
  left_join(hybrid_class, by = "sample_core")

# Optional: after join, you can remove sample_core/sample_name if you want
trip_data <- trip_data %>%
  select(-sample_name, -sample_core)

# Step 5: Batch correction per protein_id
trip_data_corrected <- trip_data %>%
  group_by(protein_id) %>%
  mutate(
    residual = {
      batches <- factor(tmt_batch)
      if (nlevels(batches) > 1) {
        model <- lm(log_abundance ~ batches)
        resid(model)
      } else {
        log_abundance  # No correction possible, keep raw
      }
    }
  ) %>%
  ungroup()

# Step 6: Prepare wide format (for PCA)
trip_data_wide <- trip_data_corrected %>%
  select(sample, species, season, Hybrid_class, protein_id, residual) %>%
  pivot_wider(
    names_from = protein_id,
    values_from = residual
  )

# Step 7: Prepare PCA input
pca_input <- trip_data_wide %>%
  select(-sample, -species, -season, -Hybrid_class)

# Remove constant or NA columns
pca_input_clean <- pca_input %>%
  select(where(~ is.numeric(.) && var(., na.rm = TRUE) != 0)) %>%   # Remove zero-variance
  select(where(~ !any(is.na(.))))                                    # Remove any column with NA

# Step 8: Run PCA
pca_result_trip <- PCA(pca_input_clean, graph = FALSE)

# Step 9: Prepare data for ggplot
pca_scores <- as.data.frame(pca_result_trip$ind$coord) %>%
  mutate(
    sample = trip_data_wide$sample,
    species = trip_data_wide$species,
    season = trip_data_wide$season,
    Hybrid_class = trip_data_wide$Hybrid_class
  )

# Step 10: Plot PCA
pca_trip <- ggplot(pca_scores, aes(x = Dim.1, y = Dim.2, color = Hybrid_class, shape = season)) +
  geom_point(size = 3) +
  labs(
    title = "Trip Batch-Corrected Protein Abundances (n=1,889)",
    x = paste0("PC1 (", round(pca_result_trip$eig[1, 2], 1), "% variance)"),
    y = paste0("PC2 (", round(pca_result_trip$eig[2, 2], 1), "% variance)")
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave("../Figures/FigS2B_PCA_trip.pdf", plot = pca_trip, width = 8, height = 8, units = "in", dpi = 300)

print(pca_trip)
