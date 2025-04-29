# PACMAD-Rhizome-Proteomics

This repository contains data and scripts for the manuscript: "Cross-Species Rhizome Proteomics Uncovers Conserved Freezing-Tolerance Strategies In Temperate Prairie Grasses: Insights For Improving Maize Cold Tolerance"

## Overview

This study analyzes proteomics data from rhizomes of multiple grass species (Tripsacum, Andropogon, Miscanthus, Panicum, and Sorghastrum) collected in summer and winter to identify conserved and species-specific cold adaptation mechanisms. The findings highlight the convergent evolution of cold tolerance strategies across independently evolved PACMAD grasses, with particular focus on the role of LEA3 proteins and other cryoprotectants.

## Important Note About File Paths

The scripts in this repository contain hardcoded paths that reflect my original development environment. You should adjust code where necessary to match your local file structure. 

## Directory Structure

- **data/**: Raw proteomics data files and RNA-seq data from different grass species
  - Includes essential RNA-seq data files `RNAseq_DEG_maize_tissues.RData` and `RNAseq_DEG_Td.RData`
  
- **tables/**: Result tables from analyses (differential abundance, GO enrichment, etc.)
  - Contains orthogroup assignments, log2FC values, and functional annotations
  
- **scripts/**: R scripts organized by analysis stage:
  - **1_data_import/**: Scripts for importing and cleaning TMT proteomics data
  - **2_data_transformation/**: Scripts for transforming protein abundance data to long format
  - **3_differential_analysis/**: Scripts for differential abundance analysis (DAPs)
  - **4_multispecies_comparison/**: Scripts for cross-species comparisons and PCA analysis
  - **5_functional_analysis/**: Scripts for GO term enrichment and cold-category classification
  - **6_cold_categories_analysis/**: Scripts for cold response categorization of top proteins
  - **7_phylogenetic_analysis/**: Scripts for protein family tree visualization (LEA, HSP, EF1)
  - **8_rna_integration/**: Scripts for correlating protein and RNA abundance across tissues
  - **utils/**: Utility scripts for sequence extraction and annotation
  
- **phylogeny/**: RAxML tree files for HSP, LEA, and EF1 protein families
  - Contains best trees for phylogenetic analysis of key protein families
  
- **annotations/**: Functional annotations from PANNZER
  - Includes PANNZER GO term assignments for differentially accumulated proteins
  
- **sequences/**: Aligned and unaligned protein sequences
  - Contains protein sequences for LEA, HSP, and EF1 families used in phylogenetic analyses
  
- **uniprot/**: UniProt annotations for top cold-responsive proteins
  - Includes results from DIAMOND BLASTp against SwissProt

  ### Supplemental Code S1–S3 (Protein Reassignment Pipeline for *M. giganteus*)

| ID | Script | Description |
|:--:|:---|:---|
| S1 | [`1_import_Ag_Mg_data.R`](scripts/1_data_import/1_import_Ag_Mg_data.R) | Imports original Ag/Mg multiplexed peptide intensities and metadata |
| S2 | [`2_map_Mg_peptides_to_proteins.R`](scripts/1_data_import/2_map_Mg_peptides_to_proteins.R) | Re-aligns *M. giganteus* peptides to the *M. giganteus* proteome using DIAMOND |
| S3 | [`3_summarize_Mg_proteins.R`](scripts/1_data_import/3_summarize_Mg_proteins.R) | Re-aggregates peptides into protein-level summaries using a multi-tiered tie-breaker |


## Species Analyzed

- **Tripsacum dactyloides** (Eastern gamagrass)
- **Andropogon gerardii** (Big bluestem)
- **Miscanthus × giganteus** (Giant miscanthus)
- **Panicum virgatum** (Switchgrass)
- **Sorghastrum nutans** (Indian grass)

## Dependencies

### R Packages
- readxl, dplyr, tidyr, readr, stringr (data manipulation)
- ggplot2, FactoMineR, factoextra (visualization and multivariate analysis)
- ggtree, phytools, ape (phylogenetic analysis)
- topGO (Gene Ontology enrichment)

### External Tools
- seqkit (for sequence extraction)
- RAxML (for phylogenetic tree construction)
- PANNZER (for functional annotation)
- DIAMOND (for protein sequence comparison)

## Citation

[Will update soon]

## License

This project is licensed under the MIT License - see the LICENSE file for details.