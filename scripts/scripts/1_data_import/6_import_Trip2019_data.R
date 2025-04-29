# Script Title: Import and Clean Tripsacum 2019 TMT Proteomics Data
# 
# Description: Imports TMT proteomics data from Excel for Tripsacum 2019 
#              samples, cleans column names, renames samples based on TMT 
#              label mapping, and updates protein IDs to Td-FL_v1.0 annotation.
#              Includes deduplication of proteins due to previous bad annotation
#              based on unique peptide count and coverage percentage.
#
# Author: Elad Oren
# Created: 2025-01-01
# Last Modified: 2025-04-01
#
# Input: 
#   - EladOren10455277_research-for-ThuyLa10439647_TMT16plex_with2peptides.xlsx
#   - Trip2019_protein_sample_labels.txt - Sample metadata and TMT label mapping
#   - annotation_mapping_Td-FL_V0.1_vs_Td-FL_V1.0.txt - ID mapping between versions
#
# Output:
#   - 0_Trip2019_clean-data.csv - Full cleaned dataset
#   - 0_Trip2019_clean-data_removed-NA.csv - Filtered dataset without NA ratios
#
# Dependencies: readxl, dplyr, stringr
#
# Notes: This script processes multiple TMT18plex sets and combines them into
#        a single dataset. Winter/Summer ratio is calculated and p-values are 
#        adjusted using Benjamini-Hochberg method.


setwd("~/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/0_raw_data/Tripsacum2019")

# Simplified script to import TMT proteomics data, relabel, and extract essential columns
library(readxl)    # For reading Excel files
library(stringr)   # For string manipulation
library(dplyr)     # For data manipulation

rm(list=ls())

# Import the protein data from Excel
protein_data <- read_excel("EladOren10455277_research-for-ThuyLa10439647_TMT16plex_with2peptides.xlsx", 
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
print("Cleaned Column Names:")
print(cleaned_colnames)


# Load the sample label key
sample_labels <- read.delim("Trip2019_protein_sample_labels.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Create a working copy of the protein data
renamed_data <- protein_data

# Loop through and rename columns by replacing F24_<TMT> with full Sample_name (hyphens preserved)
for (i in 1:nrow(sample_labels)) {
  tmt_code <- sample_labels$TMT16pro_assignment[i]
  sample_name <- sample_labels$Sample_name[i]
  
  # Replace "F24_<TMT>" with the Sample_name
  colnames(renamed_data) <- gsub(
    paste0("F24_", tmt_code),
    sample_name,
    colnames(renamed_data)
  )
}

#clean redundant suffixes like _Sample_Summer or _Sample_Winter
colnames(renamed_data) <- gsub("_Sample_(Summer|Winter)", "", colnames(renamed_data))

# remove leftover F24_ if any
colnames(renamed_data) <- gsub("F24_", "", colnames(renamed_data))

write_csv(renamed_data, "0_Trip2019_clean-data.csv")


# Filter out NA abundance ratios
renamed_protein_data <- renamed_data %>%
  filter(!is.na(Abundance_Ratio_Winter_Summer))

## Add the Td-FL_v1.0 annotation column
# First create a mapping function that returns the new ID if available, otherwise keeps the old ID
map_annotation <- function(accession) {
  idx <- match(accession, annotations$old_id)
  if (!is.na(idx)) {
    return(annotations$new_id[idx])
  } else {
    return(accession)  # Return the original ID if no mapping exists
  }
}


#Load the annotation mapping file
annotations <- read.delim("annotation_mapping_Td-FL_V0.1_vs_Td-FL_V1.0.txt", 
                          header = FALSE, 
                          sep = "\t",
                          col.names = c("old_id", "new_id"),
                          stringsAsFactors = FALSE)


# Apply the mapping function to create a new column
renamed_protein_data$Protein_id <- sapply(renamed_protein_data$Accession, map_annotation)


#rearrange protein_id column first
renamed_protein_data <- renamed_protein_data %>%
  select(
    Protein_id, 
    everything()
  )


# Check for duplicate protein IDs
duplicate_check <- renamed_protein_data %>%
  group_by(Protein_id) %>%
  summarize(count = n())

print(paste("Number of proteins with duplicated protein_id:", 
            sum(duplicate_check$count > 1)))
print(paste("Number of duplicate entries to handle:", 
            sum(duplicate_check$count[duplicate_check$count > 1]) - length(duplicate_check$count[duplicate_check$count > 1])))

# Remove duplicates by keeping entries with the highest number of unique peptides
# and in case of ties, the highest coverage percentage
protein_data_unique <- renamed_protein_data %>%
  group_by(Protein_id) %>%
  arrange(desc(Unique_Peptides), desc(Coverage), .by_group = TRUE) %>%
  slice(1) %>%
  ungroup()

print(paste("Removed", nrow(renamed_protein_data) - nrow(protein_data_unique), 
            "duplicate entries based on protein_id"))

# Replace the original data with the de-duplicated data
renamed_protein_data <- protein_data_unique


# Save the original dataset with renamed columns
write.csv(renamed_protein_data, "0_Trip2019_clean-data_removed-NA.csv", row.names = FALSE)

