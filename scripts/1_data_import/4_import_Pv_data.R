# Script Title: Import and Clean Panicum (Pv) 2022 TMT Proteomics Data
# 
# Description: Imports TMT proteomics data from Excel for Panicum 2022 samples, 
#              cleans column names, maps TMT labels to sample names, standardizes
#              naming conventions (including consistent use of dashes in species
#              names), and reformats statistical column names for analysis.
#
# Author: Elad Oren
# Created: 2025-01-01
# Last Modified: 2025-04-01
#
# Input: 
#   - CintaRomay10474963_Pv_Set2_Fr-2_with2peptides_with-adj-p-values.xlsx
#   - Pv_tmt_labels.txt - Sample metadata and TMT label mapping
#
# Output:
#   - 0_Pv_clean-data.csv - Full cleaned dataset
#   - 0_Pv_clean-data_removed-NA.csv - Filtered dataset without missing p-values
#
# Dependencies: readxl, stringr, dplyr, readr
#
# Notes: Removes Fall/Autumn samples, converts underscore to dash in species 
#        names (e.g., "Pani_Dust" to "Pani-Dust"), and calculates protein 
#        detection statistics by species and season.

setwd("~/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/0_raw_data/Pv_2022/")

# Simplified script to import TMT proteomics data, relabel, and extract essential columns
library(readxl)    # For reading Excel files
library(stringr)   # For string manipulation
library(dplyr)     # For data manipulation
library(readr)

rm(list=ls())

# Import the protein data from Excel
protein_data <- read_excel("CintaRomay10474963_Pv_Set2_Fr-2_with2peptides_with-adj-p-values.xlsx", 
                           sheet = "Proteins")

## Clean column names
# Inspect the column names
original_colnames <- colnames(protein_data)
print("Original Column Names:")
print(original_colnames)

# Function to clean column names
clean_colnames <- function(colnames) {
  # Remove ' n/a'
  colnames <- gsub(" n/a", "", colnames)
  # Replace spaces and special characters with underscores
  colnames <- gsub("[[:space:][:punct:]]+", "_", colnames)
  # Remove multiple consecutive underscores
  colnames <- gsub("_+", "_", colnames)
  # Remove leading and trailing underscores
  colnames <- gsub("^_|_$", "", colnames)
  return(colnames)
}

# Apply the cleaning function to the column names
cleaned_colnames <- clean_colnames(original_colnames)
colnames(protein_data) <- cleaned_colnames

# Verify the column name changes
print(cleaned_colnames)


# Load your key
key <- read.delim("Pv_tmt_labels.txt", stringsAsFactors = FALSE)

# Get current colnames
all_cols <- colnames(protein_data)

# Match abundance-related columns with F1_<tmt> and season
pattern <- "^(.*)_F2_(\\d+[A-Z]?)_Sample_.*_(Autumn|Fall|Winter|Summer)$"
matched <- str_match(all_cols, pattern)

# Only keep rows with a match
col_info <- data.frame(
  old_name = all_cols[!is.na(matched[, 1])],
  prefix = matched[!is.na(matched[, 1]), 2],
  tmt_label = matched[!is.na(matched[, 1]), 3],
  season = matched[!is.na(matched[, 1]), 4],
  stringsAsFactors = FALSE
)

# Normalize "Autumn" to "Fall" for matching
col_info$season[col_info$season == "Autumn"] <- "Fall"

# Join with key to get sample_name (not botanical_accession!)
col_info <- left_join(col_info, key, by = c("tmt_label", "season"))

# Check for missing matches
if (any(is.na(col_info$sample_name))) {
  warning("Some columns did not match any entry in the key file.")
  print(col_info[is.na(col_info$sample_name), ])
}

# Build new names using sample_name
col_info <- col_info %>%
  mutate(new_name = paste0(prefix, "_", sample_name, "_", season))

# Apply renaming
rename_vector <- setNames(col_info$new_name, col_info$old_name)
colnames(protein_data)[colnames(protein_data) %in% names(rename_vector)] <-
  rename_vector[colnames(protein_data)[colnames(protein_data) %in% names(rename_vector)]]

colnames(protein_data)

#rename_protein_id
protein_data <- protein_data %>%
  rename(protein_id = Accession)

# Step 2: Drop columns with _Fall or _Autumn suffix
protein_data <- protein_data %>%
  select(-matches("(_Fall|_Autumn)"))

# Step 2: First fix the species names to use dashes
colnames(protein_data) <- colnames(protein_data) %>%
  # Convert known species names with underscore to dash
  str_replace_all("Pani_Dust", "Pani-Dust")
  #str_replace_all("Pv_DD", "Pv-DD") # (if you have more, add more here)

protein_data_renamed <- protein_data

# Step 3: Now clean stat columns properly
colnames(protein_data_renamed) <- colnames(protein_data) %>%
  # Clean P-Value columns
  str_replace_all("Abundance_Ratio_P_Value_([A-Za-z0-9-]+)_([A-Za-z]+)_\\1_([A-Za-z]+)", "pval_\\1_\\2_\\3") %>%
  str_replace_all("Abundance_Ratio_Adj_P_Value_([A-Za-z0-9-]+)_([A-Za-z]+)_\\1_([A-Za-z]+)", "adj_pval_\\1_\\2_\\3") %>%
  # Clean Abundance Ratio columns
  str_replace_all("Abundance_Ratio_([A-Za-z0-9-]+)_([A-Za-z]+)_\\1_([A-Za-z]+)", "Abundance_Ratio_\\1_\\2_\\3") %>%
  # Clean grouped CV
  str_replace_all("Abundances_Grouped_CV_([A-Za-z0-9-]+)_([A-Za-z]+)", "Abundance_Grouped_CV_\\1_\\2") %>%
  # Clean grouped abundances
  str_replace_all("Abundances_Grouped_([A-Za-z0-9-]+)_([A-Za-z]+)", "Abundance_Grouped_\\1_\\2")


colnames(protein_data_renamed)


# Step 2: reorder columns — stats first, then samples
protein_data_renamed <- protein_data_renamed %>%
  select(
    protein_id, Checked, Master, Description, Coverage, Peptides, PSMs, AAs, Modifications,
    matches("^Abundance_Ratio_"),
    matches("^pval_"),
    matches("^adj_pval_"),
    matches("^Abundance_Grouped"),
    everything()  # puts all remaining columns (samples etc.) after
  )
  



# Save the original dataset with renamed columns
write_csv(protein_data_renamed, "0_Pv_clean-data.csv")

#####statisctics for proptein counts

# Step 1: Select protein_id + Abundances_Normalized columns
df_norm <- protein_data_renamed %>%
  select(protein_id, starts_with("Abundances_Normalized_"))

# Step 2: Count proteins > 0 per sample
counts_per_sample <- df_norm %>%
  select(-protein_id) %>%
  summarise(across(
    everything(),
    ~ sum(.x > 0, na.rm = TRUE)
  )) %>%
  pivot_longer(
    cols = everything(),
    names_to = "sample",
    values_to = "num_detected_proteins"
  )

# Step 3: Extract species and season
counts_per_sample <- counts_per_sample %>%
  mutate(
    species = str_extract(sample, "Pv-DD|Pv-RB|Pv|Sn"),  # Extract species (add more if needed)
    season = str_extract(sample, "Summer|Winter")        # Extract season
  )

# Step 4: Summarize by species + season
summary_per_species_season <- counts_per_sample %>%
  group_by(species, season) %>%
  summarise(
    mean_detected = mean(num_detected_proteins),
    sd_detected = sd(num_detected_proteins),
    .groups = "drop"
  )

#####


# Step 3: Drop rows with any NA in adj_pval columns
protein_data_renamed <- protein_data_renamed %>%
  filter(if_all(starts_with("adj_pval_"), ~ !is.na(.)))

write.csv(protein_data_renamed, "0_Pv_clean-data_removed-NA.csv", row.names = FALSE)
