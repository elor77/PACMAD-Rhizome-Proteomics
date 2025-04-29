# Script Title: Cross-Species Correlation Matrix of Differentially Abundant Proteins
# 
# Description: Generates a correlation matrix of log2 fold changes across species 
#              based on orthogroup-level differential abundance analysis. Computes 
#              Spearman correlations, highlights key orthogroups, and visualizes 
#              pairwise species comparisons in a multi-panel scatterplot.
#
# Author: Elad Oren
# Created: 2025-01-01
# Last Modified: 2025-04-27
#
# Input: 
#   - 3a_summary-by-OG_long_log2FC_DAP.csv - Long-format orthogroup summary file
#
# Output:
#   - ../Figures/correlation_matrix.pdf - Multi-panel figure showing species pairwise correlations
#
# Dependencies: dplyr, readr, tidyr, ggplot2, ggrepel, gridExtra, grid, ggpubr
#
# Notes: Correlation matrix is generated using lower triangular plots only. Points are 
#        colored by consistent DAP status ("Both Up", "Both Down", "Other"). Specific 
#        orthogroups (LEA3, BAM, HSP1, GST, PEPB) are labeled manually. Correlation 
#        values are calculated using Spearman's rank method and summarized separately.

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(gridExtra)
library(grid)

# Read the CSV data
data <- read.csv("3a_summary-by-OG_long_log2FC_DAP.csv", stringsAsFactors = FALSE)

# Create a wide-format data frame for pairwise comparisons
wide_data <- data %>%
  # Create separate columns for log2FC and DAP for each species
  pivot_wider(
    id_cols = orthogroup,
    names_from = species,
    values_from = c(log2FC, DAP),
    names_sep = "_"
  )

# Extract species names
species <- unique(data$species)

# Create a function to generate scatter plots for each pair of species
generate_comparison_plots <- function(data, species1, species2, x_range = c(-4.25, 5.25), y_range = c(-4.25, 5.25), highlight_orthogroup = NULL) {
  # Extract log2FC and DAP columns for the two species
  cols_needed <- c(
    "orthogroup",
    paste0("log2FC_", species1),
    paste0("log2FC_", species2),
    paste0("DAP_", species1),
    paste0("DAP_", species2)
  )
  
  # Filter for rows that have data for both species
  plot_data <- data %>%
    select(all_of(cols_needed)) %>%
    filter(complete.cases(.))
  
  # Rename columns for easier reference
  colnames(plot_data) <- c("orthogroup", "x", "y", "DAP_x", "DAP_y")
  
  # Calculate Spearman correlation
  corr <- cor(plot_data$x, plot_data$y, method = "spearman", use = "pairwise.complete.obs")
  
  # Add color based on DAP status - highlighting only if both are up or both are down
  plot_data <- plot_data %>%
    mutate(color_group = case_when(
      DAP_x == "UP" & DAP_y == "UP" ~ "Both Up",
      DAP_x == "DOWN" & DAP_y == "DOWN" ~ "Both Down",
      TRUE ~ "Other"
    ))
  
  # Filter for specified orthogroups
  highlight_data <- NULL
  if (!is.null(highlight_orthogroup)) {
    highlight_data <- plot_data %>% filter(orthogroup %in% highlight_orthogroup)
    
    # If no matches, prevent error
    if (nrow(highlight_data) == 0) {
      highlight_data <- NULL
    } else {
      # Assign manual labels for specific orthogroups
      highlight_data <- highlight_data %>%
        mutate(custom_label = case_when(
          orthogroup == "OG0022470" ~ "LEA3",
          orthogroup == "OG0015799" ~ "BAM",
          orthogroup == "OG0000288" ~ "HSP1",
          orthogroup == "OG0000642" ~ "GST",
          orthogroup == "OG0013692" ~ "PEPB",
          TRUE ~ orthogroup  # Default: keep original OG name
        ))
    }
  }
  
  # Create the scatter plot
  p <- ggplot(plot_data, aes(x = x, y = y, color = color_group)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    scale_color_manual(values = c("Both Up" = "#658FFF", "Both Down" = "#DC247F", "Other" = "#C0C0C0")) +
    labs(
      subtitle = paste("ρ =", round(corr, 2), "\nn =", nrow(plot_data)),
      x = species1,
      y = species2
    ) +
    coord_cartesian(xlim = x_range, ylim = y_range) +
    theme_classic() +
    theme(
      plot.subtitle = element_text(size = 8),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 7),
      legend.position = "none",
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")
    )
  
  # If highlight_data is not empty, add labels
  if (!is.null(highlight_data)) {
    p <- p + 
      geom_text_repel(
        data = highlight_data,
        aes(x = x, y = y, label = custom_label),
        nudge_y = 0.5,
        size = 3,
        color = "black",
        segment.color = "black",
        segment.size = 0.25,
        min.segment.length = 0
      )
  }
  
  return(list(plot = p, correlation = corr, n = nrow(plot_data)))
}

# Create a matrix of plots
# We'll use a list to store them first
matrix_plots <- list()
correlation_df <- data.frame(
  Species1 = character(),
  Species2 = character(),
  Spearman_Correlation = numeric(),
  Sample_Size = numeric(),
  stringsAsFactors = FALSE
)

# Define common x and y ranges
axis_range <- c(-4.25, 5.25)

# Generate plots for lower triangular matrix
for (i in 1:length(species)) {
  for (j in 1:length(species)) {
    if (i > j) {  # Only lower triangular part
      sp1 <- species[j]  # column
      sp2 <- species[i]  # row
      
      result <- generate_comparison_plots(wide_data, sp1, sp2, axis_range, axis_range, 
                                          highlight_orthogroup = c("OG0022470","OG0015799","OG0000288","OG0000642","OG0013692"))
      
      # Store plot in matrix position
      matrix_plots[[paste(i, j, sep = "_")]] <- result$plot
      
      # Add correlation to data frame
      correlation_df <- rbind(correlation_df, data.frame(
        Species1 = sp1,
        Species2 = sp2,
        Spearman_Correlation = result$correlation,
        Sample_Size = result$n
      ))
    }
  }
}

# Print correlation summary
print("Spearman Correlation Summary:")
print(correlation_df)

# Create an empty matrix to define the layout
n <- length(species)
layout_matrix <- matrix(NA, nrow = n, ncol = n)

# Fill diagonal with text grobs
text_grobs <- list()
for (i in 1:n) {
  text_grobs[[i]] <- textGrob(species[i], gp = gpar(fontface = "bold"))
  layout_matrix[i, i] <- i
}

# Fill lower triangular with plot indices
plot_count <- length(text_grobs)
for (i in 1:n) {
  for (j in 1:n) {
    if (i > j) {
      plot_count <- plot_count + 1
      layout_matrix[i, j] <- plot_count
    }
  }
}

# Create a list with all grobs (text and plots)
all_grobs <- c(text_grobs, matrix_plots)

# Fix the indices in the layout matrix to match the actual plot positions
for (i in 1:n) {
  for (j in 1:n) {
    if (i > j) {
      # Find the corresponding plot in matrix_plots
      plot_key <- paste(i, j, sep = "_")
      plot_index <- which(names(matrix_plots) == plot_key)
      
      if (length(plot_index) > 0) {
        # Adjust the index to account for the text grobs
        layout_matrix[i, j] <- n + plot_index
      }
    }
  }
}

# Create the multiplot using grid.arrange
create_multiplot <- function() {
  # Create a list of all grobs to display
  display_grobs <- list()
  
  # Add diagonal text elements
  for (i in 1:n) {
    display_grobs[[i]] <- text_grobs[[i]]
  }
  
  # Add plots
  for (i in 1:length(matrix_plots)) {
    display_grobs[[n + i]] <- matrix_plots[[i]]
  }
  
  # Create empty plots for upper triangular
  for (i in 1:n) {
    for (j in 1:n) {
      if (i < j) {
        # This is upper triangular - should be empty
        index <- length(display_grobs) + 1
        display_grobs[[index]] <- nullGrob()
        layout_matrix[i, j] <- index
      }
    }
  }
  
  # Create the grid arrangement
  g <- grid.arrange(
    grobs = display_grobs,
    layout_matrix = layout_matrix,
    widths = unit(rep(1, n), "null"),
    heights = unit(rep(1, n), "null")
  )
  return(g)
}

# Create a combined legend
legend_data <- data.frame(
  x = 1:3,
  y = 1:3,
  color_group = c("Both Up", "Both Down", "Other")
)

legend_plot <- ggplot(legend_data, aes(x, y, color = color_group)) +
  geom_point() +
  scale_color_manual(values = c("Both Up" = "#658FFF", "Both Down" = "#DC247F", "Other" = "#C0C0C0")) +
  labs(color = "DAP Status") +
  theme_void()

# Extract the legend
legend_grob <- ggpubr::get_legend(legend_plot)

# Create the plot matrix
multi_plot <- create_multiplot()

# Final plot with legend
final_plot <- grid.arrange(multi_plot, legend_grob, heights = c(10, 1))

# Display correlations in sorted order
correlation_df <- correlation_df %>%
  arrange(desc(Spearman_Correlation))

# Print sorted correlations
print("Sorted Spearman Correlations:")
print(correlation_df)

# Also display the highlighted orthogroups data
highlight_ogs <- c("OG0022470", "OG0015799", "OG0000288", "OG0000642", "OG0013692")
og_data <- data %>%
  filter(orthogroup %in% highlight_ogs) %>%
  select(orthogroup, species, log2FC, DAP) %>%
  arrange(orthogroup, species)

print("Highlighted Orthogroups Data:")
print(og_data)


ggsave("../Figures/correlation_matrix.pdf", plot = final_plot, width = 12, height = 14, units = "in", dpi = 300  )
