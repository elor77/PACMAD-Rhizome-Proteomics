# Script Title: Correlate Rhizome Proteomic Changes with Tissue-specific mRNA Expression
#
# Description: Compares two-year mean proteomic changes (Winter/Summer log2 fold change) 
#              in Tripsacum rhizomes with cold-responsive mRNA expression across tissues 
#              in Tripsacum and maize. Calculates Spearman correlations, classifies 
#              protein and mRNA differential expression (DAPs and DEGs), highlights key 
#              targets (e.g., LEA3), and summarizes overlap statistics.
#
# Author: Elad Oren
# Created: 2025-03-20
# Last Modified: 2025-04-27
#
# Input:
#   - ../1_tables/2_Merged_Trip2019_Trip2022_log2FC_fisher.csv - Two-year mean protein abundance data
#   - RNAseq_DEG_Td.RData - Tripsacum tissue-specific mRNA differential expression data
#   - RNAseq_DEG_maize_tissues.RData - Maize tissue-specific mRNA differential expression data
#   - mRNA_DEG_tol_7dvs22C.txt - F₂ leaf mRNA differential expression table
#   - Td_Zm_conversion.txt - Ortholog mapping table between Tripsacum and maize
#
# Output:
#   - Fig6_tissue_correlations.pdf - Multi-panel scatterplot figure summarizing protein-mRNA correlations
#   - td_*_correlation_data.csv - Per-tissue correlation tables for Tripsacum
#   - zm_*_correlation_data.csv - Per-tissue correlation tables for maize
#
# Dependencies: ggplot2, ggpubr, cowplot, dplyr, ggrepel, grid
#
# Notes:
#   - Tripsacum ‘Pete’ and F₂ leaf datasets are analyzed separately.
#   - Tripsacum-to-maize mRNA matching is based on ortholog mapping.
#   - Only proteins classified as DAPs (Up/Down) are included in RNA-protein overlap statistics.
#   - The LEA3 protein (Td00001aa022594) is manually labeled in each plot where available.

library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)
library(ggrepel) 
library(grid)

rm(list=ls())

# Set your working directory
# setwd('/your/working/directory')
setwd('/Users/eo235/Library/CloudStorage/OneDrive-Personal/Projects/BucklerLab/Tripsacum_freezing_tolerance/manuscript_V2/temp/')


# Load the TdDEG and ZmDEG data
temp_env <- new.env()
load('RNAseq_DEG_Td.RData', envir = temp_env)
TdDEG <- get(ls(temp_env)[1], envir = temp_env)
rm(temp_env)

temp_env <- new.env()
load('RNAseq_DEG_maize_tissues.RData', envir = temp_env)
ZmDEG <- get(ls(temp_env)[1], envir = temp_env)
rm(temp_env)

# Load the specific leaf data file for Td
leaf_DEG <- read.table(file = 'mRNA_DEG_tol_7dvs22C.txt',
                       sep = '\t', header = TRUE, quote = '', stringsAsFactors = FALSE)


# 1. rename the existing leaf DEGs to make them explicit
TdDEG$leaf_Pete <- TdDEG$leaf   # keep Pete seedling (5 °C) version

# 2. add the F2 leaf data as a new element
TdDEG$leaf_F2   <- leaf_DEG     # 7 °C data you read from file

# you now have both:
names(TdDEG)
# [1] "leaf_Pete" "leaf_F2" "top-crown" "bottom-crown" "root"


# # Replace leaf data in TdDEG with the specific file data
# # First, make sure the column names match
# if(!identical(colnames(TdDEG$leaf), colnames(leaf_DEG))) {
#   # Check if leaf_DEG has the required columns
#   required_cols <- c("log2FoldChange", "padj")
#   if(all(required_cols %in% colnames(leaf_DEG))) {
#     # We can replace the leaf data and maintain required structure
#     TdDEG$leaf <- leaf_DEG
#   } else {
#     # Manually ensure required columns exist
#     warning("Column names in leaf_DEG don't match TdDEG$leaf structure. Attempting to adapt.")
#     # We'd need to adapt/rename columns here if needed
#   }
# } else {
#   # Direct replacement if column structures match
#   TdDEG$leaf <- leaf_DEG
# }

# # Add this after loading the custom leaf data
# cat("Checking custom leaf data for LEA3:\n")
# if ("Td00001aa022594" %in% rownames(leaf_DEG)) {
#   cat("  Found LEA3 in rownames\n")
#   cat("  log2FC:", leaf_DEG["Td00001aa022594", "log2FoldChange"], "\n")
# } else {
#   cat("  LEA3 not found in rownames\n")
#   # Check if it might exist with a different format
#   possible_matches <- grep("022594", rownames(leaf_DEG), value = TRUE)
#   if (length(possible_matches) > 0) {
#     cat("  Possible matches in rownames:", paste(possible_matches, collapse=", "), "\n")
#   }
#   
#   # Try looking in a protein_ID column if it exists
#   if ("protein_ID" %in% colnames(leaf_DEG)) {
#     if ("Td00001aa022594" %in% leaf_DEG$protein_ID) {
#       cat("  Found LEA3 in protein_ID column\n")
#       row_index <- which(leaf_DEG$protein_ID == "Td00001aa022594")
#       cat("  log2FC:", leaf_DEG$log2FoldChange[row_index], "\n")
#     }
#   }
# }
# 
# # After replacing the leaf data in TdDEG
# cat("\nVerifying LEA3 in updated TdDEG$leaf:\n")
# if ("Td00001aa022594" %in% rownames(TdDEG$leaf)) {
#   cat("  Found LEA3 in TdDEG$leaf rownames\n")
# } else {
#   cat("  LEA3 not found in TdDEG$leaf rownames\n")
# }

# Load the gene ID conversion dictionary
Td_Zm_conversion <- read.delim("Td_Zm_conversion.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Create a mapping from Td IDs to Zm IDs
Td_to_Zm <- setNames(Td_Zm_conversion$Zm_ID, Td_Zm_conversion$Td_ID)
Zm_to_Td <- setNames(Td_Zm_conversion$Td_ID, Td_Zm_conversion$Zm_ID)

# Remove duplicates (keep the first occurrence of each ID)
Td_to_Zm <- Td_to_Zm[!duplicated(names(Td_to_Zm))]
Zm_to_Td <- Zm_to_Td[!duplicated(names(Zm_to_Td))]





proteins_mean <- read.csv(file = '../1_tables/2_Merged_Trip2019_Trip2022_log2FC_fisher.csv',
                          header = TRUE, quote = '', stringsAsFactors = FALSE)

# 1. Remove the _Txxx part
proteins_mean$short_protein_id <- sub("_T\\d+$", "", proteins_mean$protein_id)

# 2. Check for duplicates
any_duplicated <- any(duplicated(proteins_mean$short_protein_id))

# 3. (Optional) List duplicated IDs if there are any
duplicated_ids <- unique(proteins_mean$short_protein_id[duplicated(proteins_mean$short_protein_id)])

proteins_mean <- proteins_mean %>%
  distinct(short_protein_id, .keep_all = TRUE)

# Function to create correlation plots for each tissue
create_correlation_plot <- function(tissue_name, gene_dataset, dataset_name, proteins_mean, is_zm = FALSE,  y_label = tissue_name) {
  if (is_zm) {
    # For ZmDEG, map protein IDs to Zm IDs first
    zm_ids <- Td_to_Zm[proteins_mean$short_protein_id]
    zm_ids <- zm_ids[!is.na(zm_ids)]  # Remove NA mappings
    
    # Match RNA data with Zm IDs
    rna_data <- gene_dataset[[tissue_name]][match(zm_ids, 
                                                  rownames(gene_dataset[[tissue_name]])), ]
    
    # Create a combined dataframe with mapped IDs
    combined_data <- data.frame(
      gene_td = names(zm_ids),
      gene_zm = zm_ids,
      log2FC_protein = proteins_mean$mean_log2FC[match(names(zm_ids), proteins_mean$short_protein_id)],
      log2FC_RNA = rna_data$log2FoldChange,
      padj_RNA = rna_data$padj,
      DAP = proteins_mean$DAP[match(names(zm_ids), proteins_mean$short_protein_id)]
    )
  } else {
    # For TdDEG, directly match protein IDs
    rna_data <- gene_dataset[[tissue_name]][match(proteins_mean$short_protein_id, 
                                                  rownames(gene_dataset[[tissue_name]])), ]
    
    # Create a combined dataframe
    combined_data <- data.frame(
      gene = proteins_mean$short_protein_id,
      log2FC_protein = proteins_mean$mean_log2FC,
      log2FC_RNA = rna_data$log2FoldChange,
      padj_RNA = rna_data$padj,
      DAP = proteins_mean$DAP
    )
  }
  
  # Remove rows with NA values
  combined_data <- combined_data[!is.na(combined_data$log2FC_RNA), ]
  
  # Classify RNA as differentially expressed genes (DEGs)
  combined_data$DEG <- ifelse(combined_data$log2FC_RNA >= 1 & 
                                combined_data$padj_RNA <= 0.05, "Up", 
                              ifelse(combined_data$log2FC_RNA <= -1 & 
                                       combined_data$padj_RNA <= 0.05, 'Down', 'NS'))
  
  # Highlight genes that are consistent between protein and RNA
  combined_data$highlight <- ifelse(combined_data$DEG == 'Up' & combined_data$DAP == 'Up', "Up",
                                    ifelse(combined_data$DEG == 'Down' & combined_data$DAP == 'Down', "Down", "NS"))
  
  # Create the correlation plot
  p <- ggplot(combined_data, aes(x = log2FC_protein, y = log2FC_RNA, color = highlight)) +
    # Add horizontal dashed lines at y = +1 and y = -1
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", size = 0.25) +
    geom_hline(yintercept = -1, linetype = "dashed", color = "black", size = 0.25) +
    geom_point(alpha = 0.5, size = 3) +
    scale_color_manual(values = c("Up" = "#648FFF", "Down" = "#DC267F", "NS" = "#CCCCCC")) +
    xlab('Tripsacum two-year mean log2FC(Winter/Summer)') +
    ylab(paste(tissue_name, dataset_name, 'log2FC (5C/22C)')) +
    theme_classic() +
    theme(legend.position = 'none',
          text = element_text(size = 10, face = 'bold'),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste(dataset_name, tissue_name)) +
    stat_cor(aes(x = log2FC_protein, y = log2FC_RNA),
             method = "spearman",
             inherit.aes = FALSE,
             label.x = c(-3), label.y = c(8),
             cor.coef.name = expression(rho)) +
    annotate("text", x = -3, y = 7, 
             label = paste0("n=", format(nrow(combined_data), big.mark=",")), 
             hjust = 0, size = 3.5) +
    # Set fixed x-axis limits
    xlim(-3.5, 5.5) + 
    # Set y-axis limits v=based on dataset
    ylim(if(dataset_name == "Td") c(-6.5, 8.5) else c(-6.5, 10.5))
 
  # Prepare labels using geom_text_repel instead of annotate
  
  # Create a data frame for the points to label
  label_data <- data.frame(
    x = numeric(0),
    y = numeric(0),
    label = character(0)
  )
  
  # Always try to add LEA3 label (present in both Td and Zm datasets)
  lea3_col <- if(dataset_name == "Td") "gene" else "gene_td"
  lea3_row <- which(combined_data[[lea3_col]] == "Td00001aa022594")
  if (length(lea3_row) > 0) {
    label_data <- rbind(label_data, data.frame(
      x = combined_data$log2FC_protein[lea3_row],
      y = combined_data$log2FC_RNA[lea3_row],
      label = "LEA3"
    ))
  }
  
  # Add the repelled labels if any exist
  if (nrow(label_data) > 0) {
    p <- p + 
      ggrepel::geom_text_repel(
        data = label_data,
        aes(x = x, y = y, label = label),
        size = 3.5,               
        fontface = "bold",        
        color = "black",          
        segment.color = "black",  
        segment.size = 0.25,      
        min.segment.length = 0,   
        #box.padding = 0.5,       
        #point.padding = 0.3,     
        #force = 2,               )
        nudge_y = -0.3,           
        nudge_x = 0.3,            
        direction = "both"        
      )
  } 

  return(list(plot = p, data = combined_data))
}

# List of tissues to analyze
tissues <- c("leaf", "top-crown", "bottom-crown", "root")

# element names in TdDEG
leaf_sources_td <- c(F2 = "leaf_F2", Pete = "leaf_Pete")

# Create plots for TdDEG data
td_plots <- list()
for (tissue in tissues) {
  
  if (tissue == "leaf") {                          # two data sets for leaf
    for (src in names(leaf_sources_td)) {
      tissue_key <- leaf_sources_td[src]           # "leaf_F2" or "leaf_Pete"
      
      res <- create_correlation_plot(
        tissue_name   = tissue_key,         # points to correct TdDEG slot
        gene_dataset  = TdDEG,
        dataset_name  = "Td",
        proteins_mean = proteins_mean,
        is_zm         = FALSE,
        y_label       = paste0("leaf (", src, ")") # label becomes “leaf (F2)” etc.
      )
      td_plots[[paste(tissue, src, sep = "_")]] <- res$plot
      assign(paste0("td_", tissue, "_", src, "_data"), res$data)
    }
    
  } else {                                         # one data set for non-leaf tissues
    res <- create_correlation_plot(
      tissue_name   = tissue,
      gene_dataset  = TdDEG,
      dataset_name  = "Td",
      proteins_mean = proteins_mean,
      is_zm         = FALSE
    )
    td_plots[[tissue]] <- res$plot
    assign(paste0("td_", gsub("-", "_", tissue), "_data"), res$data)
  }
}

# Create plots for ZmDEG data
zm_plots <- list()
for (tissue in tissues) {
  result <- create_correlation_plot(tissue, ZmDEG, "Zm", proteins_mean, is_zm = TRUE)
  zm_plots[[tissue]] <- result$plot
  # Save the data for further analysis if needed
  assign(paste0("zm_", gsub("-", "_", tissue), "_data"), result$data)
}

# ── Row 1 ────────────────────────────────────────────────
row1 <- plot_grid(
  td_plots$leaf_F2,        # A
  nullGrob(), nullGrob(), nullGrob(),
  labels = c("A", "", "", ""),
  ncol  = 4,
  align = "hv"
)

# ── Row 2 (Tripsacum ‘Pete’ tissues) ─────────────────────
row2 <- plot_grid(
  td_plots$leaf_Pete,
  td_plots$`top-crown`,
  td_plots$`bottom-crown`,
  td_plots$root,
  labels = c("B","C","D","E"),
  ncol  = 4,
  align = "hv"
)

# ── Row 3 (maize tissues) ────────────────────────────────
row3 <- plot_grid(
  zm_plots$leaf,
  zm_plots$`top-crown`,
  zm_plots$`bottom-crown`,
  zm_plots$root,
  labels = c("F","G","H","I"),
  ncol  = 4,
  align = "hv"
)

# ── Assemble the final figure ────────────────────────────
fig6 <- plot_grid(
  row1, row2, row3,
  ncol        = 1,        # stack rows vertically
  rel_heights = c(1, 1, 1)
)

ggsave("Fig6_tissue_correlations.pdf",
       fig6, width = 20, height = 12)

print(fig6)


# Save all Td correlation data
write.csv(td_leaf_F2_data, "td_leaf_F2_data_correlation_data.csv", row.names = FALSE)
write.csv(td_leaf_Pete_data, "td_leaf_Pete_data_correlation_data.csv", row.names = FALSE)
write.csv(td_top_crown_data, "td_top_crown_Pete_correlation_data.csv", row.names = FALSE)
write.csv(td_bottom_crown_data, "td_bottom_crown_Pete_correlation_data.csv", row.names = FALSE)
write.csv(td_root_data, "td_root_Pete_correlation_data.csv", row.names = FALSE)

# Save all Zm correlation data
write.csv(zm_leaf_data, "zm_leaf_correlation_data.csv", row.names = FALSE)
write.csv(zm_top_crown_data, "zm_top_crown_correlation_data.csv", row.names = FALSE)
write.csv(zm_bottom_crown_data, "zm_bottom_crown_correlation_data.csv", row.names = FALSE)
write.csv(zm_root_data, "zm_root_correlation_data.csv", row.names = FALSE)

# ───────────────────────────────────────────────────────────────
#  Helper: how many DAPs (protein Up/Down) are also DEGs at RNA?
# ───────────────────────────────────────────────────────────────
calc_overlap <- function(df) {
  # rows classified by create_correlation_plot()
  # DAP  column: "Up", "Down", "NS"
  # DEG  column: "Up", "Down", "NS"
  
  dap_rows   <- df %>% dplyr::filter(DAP != "NS")
  n_dap      <- nrow(dap_rows)
  
  overlap    <- dap_rows %>% dplyr::filter(DEG != "NS")
  n_overlap  <- nrow(overlap)
  
  pct_overlap <- if (n_dap == 0) NA_real_ else round(n_overlap / n_dap * 100, 1)
  
  tibble::tibble(n_dap, n_overlap, pct_overlap)
}

td_summary <- bind_rows(
  leaf_F2        = calc_overlap(td_leaf_F2_data),
  leaf_Pete      = calc_overlap(td_leaf_Pete_data),
  top_crown      = calc_overlap(td_top_crown_data),
  bottom_crown   = calc_overlap(td_bottom_crown_data),
  root           = calc_overlap(td_root_data),
  .id = "tissue"
)

td_summary


zm_summary <- bind_rows(
  leaf_5C      = calc_overlap(zm_leaf_data),
  top_crown    = calc_overlap(zm_top_crown_data),
  bottom_crown = calc_overlap(zm_bottom_crown_data),
  root         = calc_overlap(zm_root_data),
  .id = "tissue"
)

zm_summary
