# Script Title: Analyze Tripsacum 2022 Differentially Abundant Proteins
# 
# Description: Performs differential abundance analysis on Tripsacum 2022 
#              proteomics data across multiple TMT batches. Combines p-values 
#              across batches using Fisher's method, calculates log2 fold change,
#              and classifies proteins as up-regulated, down-regulated, or
#              non-significant based on threshold criteria.
#
# Author: Elad Oren
# Created: 2025-01-01
# Last Modified: 2025-04-01
#
# Input: 
#   - 0_Trip2022_clean-data_removed-NA.csv - Filtered Tripsacum 2022 dataset
#
# Output:
#   - 0_Trip2022_combined_batch_ratio_filtered_NA.csv - Intermediate file with combined batches
#   - 1_Trip2022_summary_log2FC_DAP.csv - Final DAP analysis results
#
# Dependencies: readxl, dplyr, readr, tidyr
#
# Notes: Performs batch-specific p-value adjustment using Benjamini-Hochberg
#        method, then combines adjusted p-values across batches using Fisher's
#        method. Only includes proteins detected in at least 3 TMT batches
#        to ensure reliability. Classifies proteins as differentially abundant
#        using thresholds of |log2FC| ≥ 1 and padj ≤ 0.05

setwd("~/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/0_raw_data/Tripsacum2022/")

library(readxl)
library(dplyr)

rm(list=ls())


combined_protein_data <- read_csv("0_Trip2022_clean-data_removed-NA.csv")


# Extract essential columns and save the smaller dataset
combined_protein_data <- combined_protein_data %>%
  dplyr::select(
    tmt_batch,
    protein_id,
    Abundance_Ratio_Winter_Summer,
    pval
  )

combined_protein_data <- combined_protein_data %>%
  group_by(tmt_batch) %>%
  mutate(
    padj = p.adjust(pval, method = "BH")
  ) %>%
  ungroup()

combined_protein_data <- combined_protein_data %>%
  filter(!is.na(pval))


write_csv(combined_protein_data, "0_Trip2022_combined_batch_ratio_filtered_NA.csv")


# Reshape to wide format
df_wide <- combined_protein_data %>%
  select(protein_id, tmt_batch, padj) %>%
  pivot_wider(names_from = tmt_batch, values_from = padj, names_prefix = "batch_")

# Compute Fisher's combined p-value
fisher_pval <- function(pvals) {
  pvals <- pvals[!is.na(pvals)]
  if (length(pvals) == 0) return(NA)
  X2 <- -2 * sum(log(pvals))
  df <- 2 * length(pvals)
  pchisq(X2, df = df, lower.tail = FALSE)
}

df_wide <- df_wide %>%
  rowwise() %>%
  mutate(fisher_combined_padj = fisher_pval(c_across(starts_with("batch_")))) %>%
  ungroup()



# --- 4. Identify proteins found in ≥ 3 sets ---
protein_set_counts <- combined_protein_data %>%
  group_by(protein_id) %>%
  summarize(n_sets = n_distinct(tmt_batch), .groups = "drop") %>%
  filter(n_sets >= 3)


# Filter data to keep only those proteins
filtered_df <- combined_protein_data %>%
  filter(protein_id %in% protein_set_counts$protein_id)


# --- 6. Summarize: mean ratio + combined p-value ---
summary_df <- filtered_df %>%
  group_by(protein_id) %>%
  summarize(
    mean_ratio = mean(Abundance_Ratio_Winter_Summer, na.rm = TRUE),
    .groups = "drop"
  )

# Join combined p-values into summary
summary_df <- summary_df %>%
  left_join(df_wide %>% select(protein_id, fisher_combined_padj), by = "protein_id")


# --- 8. Compute log2FC, -log10(padj), and DAP status ---
summary_df <- summary_df %>%
  mutate(
    log2FC = log2(mean_ratio),
    neg_log10padj = -log10(fisher_combined_padj),
    DAP = case_when(
      log2FC >= 1 & fisher_combined_padj <= 0.05 ~ "Up",
      log2FC <= -1 & fisher_combined_padj <= 0.05 ~ "Down",
      TRUE ~ "NS"
    )
  )

# --- 9. Count DAP categories ---
up_regulated <- sum(summary_df$DAP == "Up", na.rm = TRUE)
down_regulated <- sum(summary_df$DAP == "Down", na.rm = TRUE)
not_significant <- sum(summary_df$DAP == "NS", na.rm = TRUE)

cat("Up-regulated proteins (log2FC ≥ 1, padj ≤ 0.05):", up_regulated, "\n")
cat("Down-regulated proteins (log2FC ≤ -1, padj ≤ 0.05):", down_regulated, "\n")
cat("Non-significant proteins:", not_significant, "\n")

# --- 10. Save simplified output ---
summary_df_short <- summary_df %>%
  select(protein_id, log2FC, padj = fisher_combined_padj, neg_log10padj, DAP) %>%
  filter(!is.na(padj))


write.csv(summary_df_short, "1_Trip2022_summary_log2FC_DAP.csv", row.names = FALSE)


