# Script Title: Generic Differentially Abundant Protein (DAP) Analysis Function
# 
# Description: Creates a reusable function to analyze differentially abundant
#              proteins in proteomics data. The function identifies DAPs based
#              on log2 fold change and adjusted p-value thresholds, creates
#              volcano plot variables, and outputs a simplified dataset with
#              DAP classifications. Includes extensive error handling and
#              debugging information.
#
# Author: Elad Oren
# Created: 2025-01-01
# Last Modified: 2025-04-01
#
# Input: 
#   - Any CSV with protein abundance data that includes:
#     * Protein ID column
#     * Abundance ratio column (e.g., "Abundance_Ratio_Winter_Summer")
#     * Adjusted p-value column
#
# Output:
#   - CSV with simplified data containing:
#     * Protein IDs
#     * log2FC values
#     * Adjusted p-values
#     * -log10(padj) for volcano plots
#     * DAP classification ("Up", "Down", or "NS")
#
# Dependencies: dplyr, readr
#
# Usage Example:
#   result <- analyze_dap_debug(
#     input_file = "/path/to/input.csv",
#     output_file = "/path/to/output.csv",
#     protein_col = "protein_id",
#     ratio_col = "Abundance_Ratio_Winter_Summer",
#     padj_col = "adj_pval_Winter_Summer",
#     log2fc_threshold = 1,
#     padj_threshold = 0.05
#   )


library(dplyr)
library(readr)

#' Analyze Differentially Abundant Proteins with enhanced debugging
#' 
#' @param input_file Path to input CSV file
#' @param output_file Path to output CSV file
#' @param protein_col Name of the protein ID column
#' @param ratio_col Name of the ratio column (mean ratio)
#' @param padj_col Name of the adjusted p-value column
#' @param log2fc_threshold Threshold for significant log2FC
#' @param padj_threshold Threshold for significant padj
#' @param print_summary Whether to print summary statistics
#' @return Processed dataframe with DAP annotations
analyze_dap_debug <- function(
    input_file,
    output_file,
    protein_col = "protein_id", 
    ratio_col = "Abundance_Ratio_Winter_Summer",
    padj_col = "adj_pval_Winter_Summer",
    log2fc_threshold = 1,
    padj_threshold = 0.05,
    print_summary = TRUE
) {
  # Print all parameters for debugging
  cat("Function parameters:\n")
  cat("  input_file:", input_file, "\n")
  cat("  output_file:", output_file, "\n")
  cat("  protein_col:", protein_col, "\n")
  cat("  ratio_col:", ratio_col, "\n")
  cat("  padj_col:", padj_col, "\n")
  cat("  log2fc_threshold:", log2fc_threshold, "\n")
  cat("  padj_threshold:", padj_threshold, "\n")
  
  # Check if input file exists
  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }
  
  # Try reading the file with error handling
  tryCatch({
    # Read input data
    cat("Attempting to read data from:", input_file, "\n")
    summary_df <- read_csv(input_file)
    cat("Successfully read the file. Dimensions:", nrow(summary_df), "rows x", ncol(summary_df), "columns\n")
    
    # Print column names for debugging
    cat("Column names in the file:\n")
    print(colnames(summary_df))
    
    # Check if required columns exist
    required_cols <- c(protein_col, ratio_col, padj_col)
    missing_cols <- required_cols[!required_cols %in% colnames(summary_df)]
    if (length(missing_cols) > 0) {
      stop("Missing columns in data: ", paste(missing_cols, collapse = ", "))
    }
    
    # Check for NAs in ratio column
    na_ratio_count <- sum(is.na(summary_df[[ratio_col]]))
    if (na_ratio_count > 0) {
      cat("Warning:", na_ratio_count, "NA values found in", ratio_col, "column\n")
    }
    
    # Check for zeros in ratio column (which would cause problems with log2)
    zero_ratio_count <- sum(summary_df[[ratio_col]] == 0, na.rm = TRUE)
    if (zero_ratio_count > 0) {
      cat("Warning:", zero_ratio_count, "zero values found in", ratio_col, "column. This will cause -Inf in log2FC\n")
    }
    
    # Compute log2FC, -log10(padj), and DAP status
    cat("Computing log2FC, -log10(padj), and DAP status\n")
    summary_df <- summary_df %>%
      mutate(
        log2FC = log2(!!sym(ratio_col)),
        neg_log10padj = -log10(!!sym(padj_col)),
        DAP = case_when(
          log2FC >= log2fc_threshold & !!sym(padj_col) <= padj_threshold ~ "Up",
          log2FC <= -log2fc_threshold & !!sym(padj_col) <= padj_threshold ~ "Down",
          TRUE ~ "NS"
        )
      )
    
    # Count DAP categories and check for NAs
    up_regulated <- sum(summary_df$DAP == "Up", na.rm = TRUE)
    down_regulated <- sum(summary_df$DAP == "Down", na.rm = TRUE)
    not_significant <- sum(summary_df$DAP == "NS", na.rm = TRUE)
    na_dap <- sum(is.na(summary_df$DAP))
    
    # Print summary
    cat("\nDAP Analysis Results:\n")
    cat("Up-regulated proteins (log2FC ≥", log2fc_threshold, ", padj ≤", padj_threshold, "):", up_regulated, "\n")
    cat("Down-regulated proteins (log2FC ≤", -log2fc_threshold, ", padj ≤", padj_threshold, "):", down_regulated, "\n")
    cat("Non-significant proteins:", not_significant, "\n")
    if (na_dap > 0) {
      cat("NA DAP values (potential issues):", na_dap, "\n")
    }
    
    # Check for Inf/-Inf in log2FC
    inf_count <- sum(is.infinite(summary_df$log2FC), na.rm = TRUE)
    if (inf_count > 0) {
      cat("Warning:", inf_count, "Inf/-Inf values in log2FC\n")
    }
    
    # Prepare simplified output
    cat("Preparing simplified output\n")
    summary_df_short <- summary_df %>%
      select(
        !!sym(protein_col), 
        log2FC, 
        padj = !!sym(padj_col), 
        neg_log10padj, 
        DAP
      ) %>%
      filter(!is.na(padj))
    
    cat("Final data dimensions:", nrow(summary_df_short), "rows x", ncol(summary_df_short), "columns\n")
    
    # Check output directory
    output_dir <- dirname(output_file)
    if (!dir.exists(output_dir)) {
      cat("Warning: Output directory does not exist:", output_dir, "\n")
      cat("Attempting to create directory...\n")
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Save simplified output with error handling
    cat("Attempting to write results to:", output_file, "\n")
    tryCatch({
      write_csv(summary_df_short, output_file)
      cat("Successfully wrote output file\n")
    }, error = function(e) {
      cat("Error writing output file:", conditionMessage(e), "\n")
    })
    
    # Show first few rows of output for verification
    cat("\nFirst few rows of output:\n")
    print(head(summary_df_short))
    
    # Assign result to global environment for inspection
    assign("dap_result", summary_df_short, envir = .GlobalEnv)
    cat("\nOutput data has been assigned to 'dap_result' in the global environment\n")
    
    # Return the processed dataframe
    return(summary_df_short)
  }, error = function(e) {
    cat("Error in function:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# Script example with full file paths
# Copy and paste this entire script into R, then run the example below

# Specify your input and output files with full paths
input_file <- "/Users/eo235/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/0_raw_data/Sn_2022/0_Sn_clean-data_removed-NA.csv"
output_file <- "/Users/eo235/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/1_tables/1_Sn_summary_log2FC_DAP.csv"

# For your Mg_Ag_peptides.csv file, adjust these column names as needed
result <- analyze_dap_debug(
  input_file = input_file,
  output_file = output_file,
  protein_col = "protein_id",  # Adjust this if your column name is different
  ratio_col = "Abundance_Ratio_Sorg_Winter_Summer",  # Adjust this if your column name is different
  padj_col = "adj_pval_Sorg_Winter_Summer",  # Adjust this if your column name is different
  log2fc_threshold = 1,
  padj_threshold = 0.05
)

# Check if result was successful
if (!is.null(result)) {
  cat("\nFunction completed successfully. Check 'dap_result' variable.\n")
} else {
  cat("\nFunction failed. See error messages above.\n")
}
