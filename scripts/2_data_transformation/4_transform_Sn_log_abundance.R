# Script Title: Transform Sorghastrum Abundance Data to Long Format
# 
# Description: Transforms Sorghastrum normalized abundance data from wide to long
#              format for visualization and analysis. Extracts sample and season
#              information from column names and adds a log-transformed abundance
#              column suitable for downstream analyses like PCA.
#
# Author: Elad Oren
# Created: 2025-01-01
# Last Modified: 2025-04-01
#
# Input: 
#   - 0_Sn_clean-data_removed-NA.csv - Filtered Sorghastrum dataset
#
# Output:
#   - 0_Sn_log-normalized-abundances_long-format.csv - Long-format data
#     with log-transformed abundances
#
# Dependencies: readr, dplyr, tidyr
#
# Notes: Adds a small constant (1e-6) before log transformation to avoid 
#        log(0) issues. Extracts sample names and season information using
#        regular expression pattern matching.

setwd("~/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/0_raw_data/Sn_2022/")


library(readr)
library(dplyr)
library(tidyr)

rm(list=ls())



# Step 1: Load your data
df <- read_csv("0_Sn_clean-data_removed-NA.csv")

# Step 2: Select protein_id + normalized abundance columns 
df_norm <- df %>%
  select(
    protein_id,
    matches("^Abundances_Normalized_Sn-")
  )

# Step 3: Pivot to long format and drop NAs
df_norm_long <- df_norm %>%
  pivot_longer(
    cols = -protein_id,
    names_to = "raw_name",
    values_to = "normalized_abundance"
  ) %>%
  filter(!is.na(normalized_abundance)) %>%
  extract(
    col = raw_name,
    into = c("sample", "season"),
    regex = "Abundances_Normalized_(.+)_(Summer|Winter)"
  )

# Step 4: Add log-transformed abundance column
df_norm_long <- df_norm_long %>%
  mutate(log_abundance = log(normalized_abundance + 1e-6))

# Step 5: Write the processed data to CSV
write_csv(df_norm_long, "0_Sn_log-normalized-abundances_long-format.csv")

