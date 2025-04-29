# Script Title: Merge and Analyze Tripsacum 2019–2022 Differentially Abundant Proteins
# 
# Description: Merges Tripsacum 2019 and 2022 proteomics datasets, calculates 
#              Spearman correlation between years, combines p-values using Fisher's 
#              method, and classifies proteins as up-regulated, down-regulated, or 
#              non-significant based on threshold criteria.
#
# Author: Elad Oren
# Created: 2025-01-01
# Last Modified: 2025-04-27
#
# Input: 
#   - 1_Trip2019_summary_log2FC_DAP.csv - Tripsacum 2019 DAP summary
#   - 1_Trip2022_summary_log2FC_DAP.csv - Tripsacum 2022 DAP summary
#
# Output:
#   - 2_Merged_Trip2019_Trip2022_log2FC_fisher.csv - Combined DAP summary across years
#
# Dependencies: dplyr, readr
#
# Notes: Calculates mean log2 fold change across years and adjusts combined p-values 
#        using the Benjamini-Hochberg method after Fisher's method. Classifies proteins 
#        as differentially abundant using thresholds of |mean log2FC| ≥ 1 and 
#        combined padj ≤ 0.05. Includes Spearman correlation analysis of log2FC values.


setwd("~/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/1_tables/")

# Load libraries
library(dplyr)
library(readr)

rm(list=ls())


# Load both datasets
df_2019 <- read_csv("1_Trip2019_summary_log2FC_DAP.csv")
df_2022 <- read_csv("1_Trip2022_summary_log2FC_DAP.csv")

# Rename columns for clarity
df_2019 <- df_2019 %>%
  rename_with(~ paste0(., "_2019"), -protein_id)

df_2022 <- df_2022 %>%
  rename_with(~ paste0(., "_2022"), -protein_id)

# Merge by protein_id
df_merged <- inner_join(df_2019, df_2022, by = "protein_id")

# --- A. Spearman correlation on log2FC ---
spearman_result <- cor.test(df_merged$log2FC_2019, df_merged$log2FC_2022, method = "spearman")
print(spearman_result)

# --- B. Fisher's method for combined adjusted p-value ---
# Function for Fisher's method
fisher_pval <- function(pvals) {
  pvals <- pmin(pmax(pvals, 1e-300), 1 - 1e-16)  # sanitize p-values
  X2 <- -2 * sum(log(pvals))
  df <- 2 * length(pvals)
  pchisq(X2, df = df, lower.tail = FALSE)
}

# Apply Fisher per row
df_merged <- df_merged %>%
  rowwise() %>%
  mutate(
    mean_log2FC = mean(c(log2FC_2019, log2FC_2022), na.rm = TRUE),
    fisher_combined_padj = fisher_pval(c(padj_2019, padj_2022))
  ) %>%
  ungroup()

df_merged <- df_merged %>%
  mutate(
    fisher_combined_padj = p.adjust(fisher_combined_padj, method = "BH")
  )

# Optional: add -log10(Fisher p) and classify DAP
df_merged <- df_merged %>%
  mutate(
    neg_log10padj = -log10(fisher_combined_padj),
    DAP = case_when(
      mean_log2FC >= 1 & fisher_combined_padj <= 0.05 ~ "Up",
      mean_log2FC <= -1 & fisher_combined_padj <= 0.05 ~ "Down",
      TRUE ~ "NS"
    )
  )



# View results
head(df_merged)

# --- Count DAP categories ---
dap_counts <- df_merged %>%
  count(DAP)

print(dap_counts)


# Optional: write to file
write_csv(df_merged, "2_Merged_Trip2019_Trip2022_log2FC_fisher.csv")
