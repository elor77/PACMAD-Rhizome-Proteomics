# Script Title: Import and Clean Tripsacum 2022 TMT Proteomics Data
# 
# Description: Imports TMT proteomics data from Excel files for Tripsacum 2022 
#              samples, cleans column names, maps TMT labels to sample names,
#              combines data from multiple sets, and filters rows with missing
#              values. Prepares data for downstream differential abundance analysis.
#
# Author: Elad Oren
# Created: 2025-01-01
# Last Modified: 2025-04-01
#
# Input: 
#   - SaraMiller10462499_TMT18plex_Set1_Fr_with2peptides.xlsx
#   - SaraMiller10462499_TMT18plex_Set2_Fr_with2peptides.xlsx
#   - SaraMiller10462499_TMT18plex_Set3_Fr_with2peptides.xlsx
#   - SaraMiller10462499_TMT18plex_Set4_Fr_with2peptides.xlsx
#   - Trip2022_TMT_sample_name_key.txt - Sample metadata and TMT label mapping
#
# Output:
#   - 0_Trip2022_clean-data_removed-NA.csv - Cleaned dataset ready for analysis
#
# Dependencies: readxl, dplyr
#
# Notes: This script processes multiple TMT18plex sets and combines them into
#        a single dataset. Winter/Summer ratio is calculated and p-values are 
#        adjusted using Benjamini-Hochberg method.

setwd("~/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/0_raw_data/Tripsacum2022/")

library(readxl)
library(dplyr)

# --- 1. File list ---
file_names <- c(
  "SaraMiller10462499_TMT18plex_Set1_Fr_with2peptides.xlsx",
  "SaraMiller10462499_TMT18plex_Set2_Fr_with2peptides.xlsx",
  "SaraMiller10462499_TMT18plex_Set3_Fr_with2peptides.xlsx",
  "SaraMiller10462499_TMT18plex_Set4_Fr_with2peptides.xlsx"
)

# --- 2. Load sample key ---
sample_key <- read.delim("Trip2022_TMT_sample_name_key.txt", stringsAsFactors = FALSE)

# --- 3. Function to rename columns by sample ---
rename_abundance_columns <- function(df, set_num) {
  set_key <- sample_key[sample_key$protein_set == set_num, ]
  
  for (i in seq_len(nrow(set_key))) {
    tmt <- set_key$TMT18plex_label[i]
    sample <- set_key$sample_name[i]
    
    # Match _<TMT>_Sample_ and replace with _<sample>_Sample_
    pattern <- paste0("_", tmt, "_Sample_")
    replacement <- paste0("_", sample, "_Sample_")
    
    colnames(df) <- gsub(pattern, replacement, colnames(df), fixed = TRUE)
  }
  
  return(df)
}

# --- 4. Read, clean, and process each file ---
protein_data_list <- lapply(seq_along(file_names), function(i) {
  df <- read_excel(file_names[i], sheet = "Proteins")
  #df <- read_excel(file_names[1], sheet = "Proteins")

    # Clean column names
  colnames(df) <- gsub(" n/a", "", colnames(df))
  colnames(df) <- gsub("[[:space:][:punct:]]+", "_", colnames(df))
  colnames(df) <- gsub("_+", "_", colnames(df))
  colnames(df) <- gsub("^_|_$", "", colnames(df))
  colnames(df) <- gsub("_F[1-9]_", "_", colnames(df))
  
  # Add protein set identifier
  df$tmt_batch <- paste0(i)
  #df$tmt_batch <- paste0(1)
  
  # Rename abundance-related columns using sample key
  df <- rename_abundance_columns(df, i)
  #df <- rename_abundance_columns(df, 1)
  colnames(df) <- gsub("_Sample_", "_", colnames(df), fixed = TRUE)
  
  # Calculate the Winter/Summer ratio
  df$Abundance_Ratio_Winter_Summer <- ifelse(
    !is.na(df$Abundance_Ratio_Summer_Winter) & df$Abundance_Ratio_Summer_Winter != 0, 
    1 / df$Abundance_Ratio_Summer_Winter, 
    NA_real_
  )
  
  # Remove all "Fall" columns
  fall_columns <- grep("Fall", colnames(df), value = TRUE)
  if(length(fall_columns) > 0) {
    df <- df[, !colnames(df) %in% fall_columns]
  }
  
  # Add protein set identifier
  df$tmt_batch <- paste0(i)
  
  # Rearrange columns
  razor_pos <- which(colnames(df) == "Razor_Peptides")
  if(length(razor_pos) > 0) {
    cols <- colnames(df)
    new_order <- c(
      "tmt_batch",
      cols[!cols %in% c("tmt_batch", "Abundance_Ratio_Winter_Summer")]
    )
    razor_pos_new <- which(new_order == "Razor_Peptides")
    final_order <- c(
      new_order[1:razor_pos_new],
      "Abundance_Ratio_Winter_Summer",
      new_order[(razor_pos_new+1):length(new_order)]
    )
    df <- df[, final_order]
  } else {
    cols <- colnames(df)
    cols <- cols[!cols %in% c("tmt_batch", "Abundance_Ratio_Winter_Summer")]
    df <- df[, c("tmt_batch", cols[1], "Abundance_Ratio_Winter_Summer", cols[-1])]
  }
  
  return(df)
})

# --- 5. Combine all sets ---
combined_protein_data <- bind_rows(protein_data_list)

# --- 6. Rename key columns ---
combined_protein_data <- combined_protein_data %>%
  rename(
    protein_id = Accession,
    pval = Abundance_Ratio_P_Value_Summer_Winter
  )


write_csv(combined_protein_data, "0_Trip2022_clean-data_removed-NA.csv")
