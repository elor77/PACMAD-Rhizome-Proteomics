# Script Title: Evolutionary Conservation Analysis of Cold-Responsive Orthogroups
# 
# Description: Analyzes evolutionary conservation of cold response across grass species
#              by comparing pairwise Spearman correlations between "Conserved Cold" 
#              orthogroups (DAP in ≥2 species) and background orthogroups. Uses relaxed
#              filtering (present in 4 of 5 lineages) with strict harmonization 
#              (same N for all metrics).
#
# Author: Elad Oren
# Created: 2025-XX-XX
# Last Modified: 2025-XX-XX
#
# Input: 
#   - 0_Ag_log-normalized-abundances_long-format.csv - Andropogon gerardii abundance data
#   - 0_Mg_3_log-normalized-abundances_long-format.csv - Miscanthus x giganteus abundance data
#   - 0_Pv-DD_log-normalized-abundances_long-format.csv - Panicum virgatum abundance data
#   - 0_Sn_log-normalized-abundances_long-format.csv - Sorghastrum nutans abundance data
#   - 0_Trip2019_log-normalized-abundances_long-format.csv - Tripsacum dactyloides 2019 abundance data
#   - 0_Trip2022_log-normalized-abundances_long-format.csv - Tripsacum dactyloides 2022 abundance data
#   - 3a_summary-by-OG_long_log2FC_DAP.csv - Log2 fold change and DAP status by orthogroup
#   - 0_Orthogroups_Ag_Mg_Pv_Sc_Sn_Td-v1.0.txt - Orthogroup-to-protein mapping
#
# Output:
#   - Evolution_Final_Orthogroups_relaxed.pdf - Boxplot comparing correlations
#   - Evolution_Statistics_Table_relaxed.csv - Summary statistics with medians and p-values
#
# Dependencies: tidyverse, ggpubr
#
# Notes: 
#   - Tripsacum 2019 and 2022 samples are merged into a single "Td" group.
#   - Relaxed filter: Orthogroups must be present in at least 4 of 5 lineages.
#   - Conserved Cold: Orthogroups with DAP status in ≥2 species.
#   - Statistical test: Wilcoxon rank-sum (Mann-Whitney U) test.

# ==============================================================================
# 1. SETUP
# ==============================================================================

# Clear environment
rm(list = ls())

# Load required libraries
library(tidyverse)    # Data manipulation and visualization
library(ggpubr)       # Publication-ready plots

# Set working directory and output path
# setwd("YOUR_WORKING_DIRECTORY")
output_dir <- "/Users/eo235/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/PACMAD-Rhizome-Proteomics/manuscript_revised/anaylsis"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Define input files
files_abundance <- list(
  Ag = "0_Ag_log-normalized-abundances_long-format.csv",
  Mg = "0_Mg_3_log-normalized-abundances_long-format.csv",
  Pv = "0_Pv-DD_log-normalized-abundances_long-format.csv",
  Sn = "0_Sn_log-normalized-abundances_long-format.csv",
  Td19 = "0_Trip2019_log-normalized-abundances_long-format.csv",
  Td22 = "0_Trip2022_log-normalized-abundances_long-format.csv"
)
file_logfc <- "3a_summary-by-OG_long_log2FC_DAP.csv"
file_og_map <- "0_Orthogroups_Ag_Mg_Pv_Sc_Sn_Td-v1.0.txt"

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

#' Filter orthogroups present in at least 4 of 5 lineages
#' @param mat Data frame with orthogroup column and species columns
#' @return Filtered data frame
keep_if_4_of_5 <- function(mat) {
  mat_vals <- mat %>% select(-orthogroup)
  n_present <- rowSums(!is.na(mat_vals))
  mat %>% filter(n_present >= 4)
}

#' Extract pairwise Spearman correlations from a matrix
#' @param data_matrix Numeric matrix (rows = orthogroups, cols = species)
#' @return Vector of upper-triangle correlation values
get_pairwise_corrs <- function(data_matrix) {
  if (ncol(data_matrix) < 2) return(numeric(0))
  corr_mat <- cor(data_matrix, method = "spearman", use = "pairwise.complete.obs")
  return(corr_mat[upper.tri(corr_mat)])
}

#' Perform Wilcoxon test and extract p-value
#' @param group_back List with 'corr' vector for background group
#' @param group_cold List with 'corr' vector for conserved cold group
#' @return P-value from Mann-Whitney U test
get_wilcox_p <- function(group_back, group_cold) {
  test_result <- wilcox.test(group_back$corr, group_cold$corr)
  return(test_result$p.value)
}

# ==============================================================================
# 3. LOAD DATA
# ==============================================================================

# Load orthogroup mapping
og_map <- read_tsv(file_og_map, show_col_types = FALSE)[, c(1, 4)]
colnames(og_map) <- c("orthogroup", "protein_id")
og_map <- og_map %>% distinct(protein_id, .keep_all = TRUE)

# ==============================================================================
# 4. PREPARE LOG2FC MATRIX (Relaxed 4/5 Filter)
# ==============================================================================

df_logfc <- read_csv(file_logfc, show_col_types = FALSE)

# Pivot to wide format
if ("species" %in% names(df_logfc)) {
  mat_logfc <- df_logfc %>%
    select(orthogroup, species, log2FC) %>%
    pivot_wider(names_from = species, values_from = log2FC)
} else {
  mat_logfc <- df_logfc %>% 
    select(orthogroup, any_of(c("Ag", "Mg", "Pv", "Sn", "Trip19", "Trip22")))
}

# Merge Tripsacum replicates
if ("Trip19" %in% names(mat_logfc) & "Trip22" %in% names(mat_logfc)) {
  mat_logfc$Td <- rowMeans(mat_logfc[, c("Trip19", "Trip22")], na.rm = TRUE)
  mat_logfc <- mat_logfc %>% select(-Trip19, -Trip22)
}

# Ensure all species columns exist and apply filter
for (col in c("Ag", "Mg", "Pv", "Sn", "Td")) {
  if (!col %in% names(mat_logfc)) mat_logfc[[col]] <- NA
}
mat_logfc <- mat_logfc %>% 
  select(orthogroup, Ag, Mg, Pv, Sn, Td) %>% 
  keep_if_4_of_5()

# ==============================================================================
# 5. PREPARE ABUNDANCE MATRICES (Relaxed 4/5 Filter)
# ==============================================================================

#' Build abundance matrix for a given season
#' @param season_name Character string: "Winter" or "Summer"
#' @return Wide-format data frame with orthogroup abundances per species
get_abundance_matrix <- function(season_name) {
  df_list <- list()
  
  for (sp in names(files_abundance)) {
    fpath <- files_abundance[[sp]]
    if (file.exists(fpath)) {
      df <- read_csv(fpath, show_col_types = FALSE)
      if ("season" %in% names(df)) {
        df_filt <- df %>% filter(tolower(season) == tolower(season_name))
        if (nrow(df_filt) > 0) {
          df_agg <- df_filt %>%
            group_by(protein_id) %>%
            summarise(abundance = mean(normalized_abundance, na.rm = TRUE), .groups = 'drop') %>%
            mutate(Species = sp)
          df_list[[sp]] <- df_agg
        }
      }
    }
  }
  
  if (length(df_list) == 0) return(NULL)
  
  mat <- bind_rows(df_list) %>%
    inner_join(og_map, by = "protein_id") %>%
    group_by(orthogroup, Species) %>%
    summarise(val = sum(abundance, na.rm = TRUE), .groups = 'drop') %>%
    pivot_wider(names_from = Species, values_from = val)
  
  # Merge Tripsacum replicates
  if ("Td19" %in% names(mat) & "Td22" %in% names(mat)) {
    mat$Td <- rowMeans(mat[, c("Td19", "Td22")], na.rm = TRUE)
    mat <- mat %>% select(-Td19, -Td22)
  }
  
  # Ensure all species columns exist and apply filter
  for (col in c("Ag", "Mg", "Pv", "Sn", "Td")) {
    if (!col %in% names(mat)) mat[[col]] <- NA
  }
  mat <- mat %>% 
    select(orthogroup, Ag, Mg, Pv, Sn, Td) %>% 
    keep_if_4_of_5()
  
  return(mat)
}

# Build seasonal abundance matrices
mat_winter <- get_abundance_matrix("Winter")
mat_summer <- get_abundance_matrix("Summer")

# ==============================================================================
# 6. HARMONIZE DATASETS (Strict Intersection)
# ==============================================================================

# Find orthogroups present in all three matrices
common_ogs <- intersect(mat_logfc$orthogroup, mat_winter$orthogroup)
common_ogs <- intersect(common_ogs, mat_summer$orthogroup)
cat("Final Harmonized Orthogroup Count:", length(common_ogs), "\n")

# Filter matrices to common orthogroups
mat_logfc_final <- mat_logfc %>% filter(orthogroup %in% common_ogs)
mat_winter_final <- mat_winter %>% filter(orthogroup %in% common_ogs)
mat_summer_final <- mat_summer %>% filter(orthogroup %in% common_ogs)

# ==============================================================================
# 7. DEFINE ORTHOGROUP GROUPS
# ==============================================================================

# Count DAP occurrences per orthogroup
dap_counts <- df_logfc %>%
  filter(orthogroup %in% common_ogs) %>%
  filter(species %in% c("Ag", "Mg", "Pv", "Sn", "Trip19", "Trip22", "Td")) %>%
  mutate(is_dap = if_else(DAP != "NS" & !is.na(DAP), 1, 0)) %>%
  group_by(orthogroup) %>%
  summarise(n_dap = sum(is_dap))

# Define groups: Conserved Cold (≥2 DAP) vs Background (≤1 DAP)
ogs_conserved <- dap_counts %>% filter(n_dap >= 2) %>% pull(orthogroup)
ogs_background <- dap_counts %>% filter(n_dap <= 1) %>% pull(orthogroup)

# ==============================================================================
# 8. CALCULATE CORRELATIONS
# ==============================================================================

#' Filter matrix by orthogroup subset and compute correlations
#' @param mat Data frame with orthogroup and species columns
#' @param og_subset Character vector of orthogroups to include
#' @return List with 'corr' (correlation values) and 'n' (sample size)
filter_and_corr <- function(mat, og_subset) {
  m <- mat %>% filter(orthogroup %in% og_subset) %>% select(-orthogroup)
  list(corr = get_pairwise_corrs(m), n = nrow(m))
}

# Calculate correlations for each metric and group
res_ab_win_c <- filter_and_corr(mat_winter_final, ogs_conserved)
res_ab_win_b <- filter_and_corr(mat_winter_final, ogs_background)
res_ab_sum_c <- filter_and_corr(mat_summer_final, ogs_conserved)
res_ab_sum_b <- filter_and_corr(mat_summer_final, ogs_background)
res_fc_c <- filter_and_corr(mat_logfc_final, ogs_conserved)
res_fc_b <- filter_and_corr(mat_logfc_final, ogs_background)

# ==============================================================================
# 9. PREPARE PLOT DATA
# ==============================================================================

plot_data <- bind_rows(
  data.frame(Metric = "Abundance (Summer)", Group = "Background", 
             Corr = res_ab_sum_b$corr, N = res_ab_sum_b$n),
  data.frame(Metric = "Abundance (Summer)", Group = "Conserved Cold", 
             Corr = res_ab_sum_c$corr, N = res_ab_sum_c$n),
  data.frame(Metric = "Abundance (Winter)", Group = "Background", 
             Corr = res_ab_win_b$corr, N = res_ab_win_b$n),
  data.frame(Metric = "Abundance (Winter)", Group = "Conserved Cold", 
             Corr = res_ab_win_c$corr, N = res_ab_win_c$n),
  data.frame(Metric = "Response (Log2FC)", Group = "Background", 
             Corr = res_fc_b$corr, N = res_fc_b$n),
  data.frame(Metric = "Response (Log2FC)", Group = "Conserved Cold", 
             Corr = res_fc_c$corr, N = res_fc_c$n)
)

# Set factor levels for plotting order
plot_data$Metric <- factor(plot_data$Metric, 
                           levels = c("Abundance (Summer)", "Abundance (Winter)", "Response (Log2FC)"))

# Create legend labels with sample sizes
n_back <- unique(plot_data$N[plot_data$Group == "Background"])
n_cons <- unique(plot_data$N[plot_data$Group == "Conserved Cold"])

label_back <- paste0("Background (n=", n_back, ")")
label_cons <- paste0("Conserved Cold (n=", n_cons, ")")

plot_data$Group <- factor(plot_data$Group, 
                          levels = c("Background", "Conserved Cold"),
                          labels = c(label_back, label_cons))

# ==============================================================================
# 10. GENERATE PLOT
# ==============================================================================

p <- ggboxplot(plot_data, x = "Metric", y = "Corr", fill = "Group", 
               palette = c("#95a5a6", "#e74c3c")) +
  stat_compare_means(aes(group = Group), label = "p.signif") +
  labs(
    y = "Spearman Correlation", 
    fill = "Group", 
    title = "Evolutionary Conservation: Seasonal Abundance vs Response", 
    subtitle = "Orthogroups conserved across 4 of 5 species. n=1683"
  ) +
  theme(
    axis.text.x = element_text(size = 11, face = "bold"), 
    legend.position = "bottom",
    legend.title = element_text(face = "bold")
  )

print(p)
ggsave(file.path(output_dir, "Evolution_Final_Orthogroups_relaxed.pdf"), 
       width = 8, height = 6)

# ==============================================================================
# 11. GENERATE STATISTICS TABLE
# ==============================================================================

final_stats_table <- data.frame(
  Metric = c("Abundance (Summer)", "Abundance (Winter)", "Response (Log2FC)"),
  
  # Sample sizes
  N_Background = c(res_ab_sum_b$n, res_ab_win_b$n, res_fc_b$n),
  N_Conserved  = c(res_ab_sum_c$n, res_ab_win_c$n, res_fc_c$n),
  
  # Median correlations
  Median_Background_rho = c(
    median(res_ab_sum_b$corr, na.rm = TRUE), 
    median(res_ab_win_b$corr, na.rm = TRUE), 
    median(res_fc_b$corr, na.rm = TRUE)
  ),
  Median_Conserved_rho = c(
    median(res_ab_sum_c$corr, na.rm = TRUE), 
    median(res_ab_win_c$corr, na.rm = TRUE), 
    median(res_fc_c$corr, na.rm = TRUE)
  ),
  
  # P-values from Wilcoxon test
  P_Value = c(
    get_wilcox_p(res_ab_sum_b, res_ab_sum_c),
    get_wilcox_p(res_ab_win_b, res_ab_win_c),
    get_wilcox_p(res_fc_b, res_fc_c)
  )
)

# Add significance stars
final_stats_table$Significance <- case_when(
  final_stats_table$P_Value < 0.001 ~ "***",
  final_stats_table$P_Value < 0.01  ~ "**",
  final_stats_table$P_Value < 0.05  ~ "*",
  TRUE ~ "ns"
)

# Print and save
print(final_stats_table)
write_csv(final_stats_table, file.path(output_dir, "Evolution_Statistics_Table_relaxed.csv"))
