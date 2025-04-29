# Script Title: Analyze Tripsacum 2019 Differentially Abundant Proteins with Tukey's Test
# 
# Description: Performs differential abundance analysis on Tripsacum 2019 
#              proteomics data using a more robust statistical approach.
#              For each protein, performs ANOVA with Tukey's test to compare
#              winter vs. summer samples, adjusts p-values using Benjamini-Hochberg
#              method, calculates log2FC, and classifies proteins as differentially
#              abundant based on threshold criteria.
#
# Author: Elad Oren
# Created: 2025-01-01
# Last Modified: 2025-04-01
#
# Input: 
#   - 0_Trip2019_clean-data_removed-NA.csv - Filtered Tripsacum 2019 dataset
#
# Output:
#   - 1_Trip2019_summary_log2FC_DAP.csv - Final DAP analysis results with
#     robust statistical testing
#
# Dependencies: tidyr, multcomp, car, dplyr
#
# Notes: Needed to recalculate ratio-based pvalues for Trip2019 data set. 
#        Uses Tukey's HSD test Requires at least 3 data points per protein to perform the test. 
#        Classifies proteins as differentially abundant using thresholds of |log2FC| ≥ 1 and padj ≤ 0.05.


# Load necessary libraries
library(tidyr)
library(multcomp)  # For Tukey's test
library(car)       # For Anova function
library(dplyr)

# Read the processed data file
protein_data <- read.csv("0_Trip2019_clean-data_removed-NA.csv")

# First identify which columns are normalized abundances for winter and summer
winter_cols <- grep("^Abundances_Normalized_.*_winter$", names(protein_data), value = TRUE)
summer_cols <- grep("^Abundances_Normalized_.*_summer$", names(protein_data), value = TRUE)

# Create an empty dataframe to store results
tukey_results <- data.frame(
  Protein_id = character(),
  Abundance_Ratio_P_Value_Winter_Summer = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each protein
for (i in 1:nrow(protein_data)) {
  # Get current protein ID
  current_protein <- protein_data$Protein_id[i]
  
  # Create a data frame for this protein's values
  protein_values <- data.frame(
    value = c(
      as.numeric(protein_data[i, winter_cols]),
      as.numeric(protein_data[i, summer_cols])
    ),
    season = c(
      rep("winter", length(winter_cols)),
      rep("summer", length(summer_cols))
    )
  )
  
  # Remove any NA values
  protein_values <- protein_values[!is.na(protein_values$value), ]
  
  # Only run the test if we have enough data points
  if (nrow(protein_values) >= 3) {
    # Fit ANOVA model
    aov_model <- aov(value ~ season, data = protein_values)
    
    # Run Tukey's test
    tukey_result <- TukeyHSD(aov_model)
    
    # Extract p-value
    tukey_pval <- tukey_result$season["winter-summer", "p adj"]
    
    # Add to results
    tukey_results <- rbind(tukey_results, 
                           data.frame(Protein_id = current_protein, 
                                      pvalue = tukey_pval))
  } else {
    # Not enough data points for this protein
    tukey_results <- rbind(tukey_results, 
                           data.frame(Protein_id = current_protein, 
                                      pvalue = NA))
  }
}

# Join the Tukey p-values back to the original data
protein_data_with_tukey <- protein_data %>%
  left_join(tukey_results, by = "Protein_id")

# Check the column names to see what we actually have
print("Column names after join:")
print(names(protein_data_with_tukey))


# Join the Tukey p-values back to the original data
protein_data_with_tukey <- protein_data_with_tukey %>%
  rename(
    original_pvalue = Abundance_Ratio_P_Value_Winter_Summer,
    tukey_pvalue = pvalue
  )

# Calculate FDR-adjusted p-values using Benjamini-Hochberg method
protein_data_with_tukey <- protein_data_with_tukey %>%
  mutate(
    padj = p.adjust(tukey_pvalue, method = "BH")
  )

# Calculate log2FC and -log10(padj)
protein_data_with_tukey <- protein_data_with_tukey %>%
  mutate(
    log2FC = log2(Abundance_Ratio_Winter_Summer),
    neg_log10padj = -log10(padj)
  )

# Add DAP (Differentially Abundant Protein) column
protein_data_with_tukey <- protein_data_with_tukey %>%
  mutate(
    DAP = ifelse(log2FC >= 1 & padj <= 0.05, "Up", 
                 ifelse(log2FC <= -1 & padj <= 0.05, "Down", "NS"))
  )

# Count proteins in each category
up_regulated <- sum(protein_data_with_tukey$DAP == "Up", na.rm = TRUE)
down_regulated <- sum(protein_data_with_tukey$DAP == "Down", na.rm = TRUE)
not_significant <- sum(protein_data_with_tukey$DAP == "NS", na.rm = TRUE)

# Print summary
print(paste("Up-regulated proteins (log2FC >= 1, padj <= 0.05):", up_regulated))
print(paste("Down-regulated proteins (log2FC <= -1, padj <= 0.05):", down_regulated))
print(paste("Non-significant proteins:", not_significant))


# Extract essential columns and save the smaller dataset
protein_data_with_tukey <- protein_data_with_tukey %>%
  dplyr::select(
    protein_id = Protein_id,
    log2FC,
    padj,
    neg_log10padj,
    DAP
    )

# Save the updated data
write.csv(protein_data_with_tukey, "1_Trip2019_summary_log2FC_DAP.csv", row.names = FALSE)





