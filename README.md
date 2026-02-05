# PACMAD-Rhizome-Proteomics

Code and data repository for the manuscript:

**"Constrained evolution of a core winter proteome across independently cold-adapted PACMAD grasses"**

Elad Oren, Jingjing Zhai, Travis E. Rooney, Ruthie Angelovici, Charles O. Hale, Lara J. Brindisi, Sheng-Kai Hsu, Christine M. Gault, Jian Hua, Thuy La, Nicholas Lepak, Qin Fu, Edward S. Buckler, M. Cinta Romay

*Corresponding authors:* Elad Oren (eladoren@volcani.agri.gov.il), M. Cinta Romay (mcr72@cornell.edu)

## Overview

This study compares seasonal rhizome proteomes (winter vs. summer) from five independently cold-adapted PACMAD grass species grown in a common garden in Ithaca, NY (USDA Hardiness Zone 6a). Using TMTpro-labeled quantitative proteomics, we show that shared cold-responsive proteins exhibit substantially higher cross-species fold-change correlation (ρ = 0.80) than background proteins (ρ = 0.45), revealing evolutionary constraint on protein-level response magnitude. LEA3 was the only ortholog elevated across all five species. Structural comparison of the Tripsacum and maize LEA3 orthologs revealed a 15-amino-acid insertion in maize that disrupts the conserved 11-mer motif architecture, suggesting that transcriptional induction alone does not ensure freezing tolerance.

## Species Analyzed

| Species | Common Name | Cultivar/Entry |
|---|---|---|
| *Tripsacum dactyloides* × *T. floridanum* | Eastern gamagrass × Florida gamagrass | Multiple hybrid genotypes |
| *Andropogon gerardii* | Big bluestem | Sentinel |
| *Miscanthus* × *giganteus* | Giant miscanthus | — |
| *Panicum virgatum* | Switchgrass | Dust Devil |
| *Sorghastrum nutans* | Indiangrass | Golden Sunset |

## Important Note About File Paths

The scripts contain hardcoded paths from the original development environment. Adjust paths as needed to match your local file structure.

## Directory Structure

```
├── data/                    # Raw proteomics and RNA-seq data files
│   ├── RNAseq_DEG_maize_tissues.RData
│   └── RNAseq_DEG_Td.RData
│
├── scripts/                 # R scripts organized by analysis stage
│   ├── 1_data_import/       # TMT proteomics data import and cleaning
│   ├── 2_data_transformation/  # Protein abundance → long format
│   ├── 3_differential_analysis/  # Differential abundance (DAPs)
│   ├── 4_multispecies_comparison/  # Cross-species PCA and correlations
│   ├── 5_functional_analysis/  # GO enrichment and cold-category classification
│   ├── 6_cold_categories_analysis/  # Cold response categorization of top proteins
│   ├── 7_phylogenetic_analysis/  # Protein family tree visualization (LEA, HSP, EF1)
│   ├── 8_rna_integration/   # Protein–RNA abundance correlation across tissues
│   └── utils/               # Sequence extraction and annotation utilities
│
├── tables/                  # Result tables (DAPs, GO enrichment, log2FC, orthogroups)
├── annotations/             # PANNZER GO term assignments for DAPs
├── sequences/               # Aligned and unaligned protein sequences (LEA, HSP, EF1)
├── uniprot/                 # DIAMOND BLASTp results against SwissProt
├── phylogeny/               # RAxML best trees for EF1, LEA, and HSP families
│
├── LEA3_structure/          # Structural modeling and analysis of LEA3 proteins
│   ├── AlphaFold_output/    # AlphaFold3 3D structure predictions (Tripsacum & maize)
│   ├── AmphipaSeek/         # Amphipathic helix identification results
│   └── ProtScale/           # Kyte-Doolittle hydropathy analysis outputs
│
└── manuscript/              # Manuscript files and figures
```

### Supplemental Code S1–S3 (Protein Reassignment Pipeline for *M. giganteus*)

| ID | Script | Description |
|:--:|:---|:---|
| S1 | [`1_import_Ag_Mg_data.R`](scripts/1_data_import/1_import_Ag_Mg_data.R) | Imports original Ag/Mg multiplexed peptide intensities and metadata |
| S2 | [`2_map_Mg_peptides_to_proteins.R`](scripts/1_data_import/2_map_Mg_peptides_to_proteins.R) | Re-aligns *M. giganteus* peptides to the *M. giganteus* proteome using DIAMOND |
| S3 | [`3_summarize_Mg_proteins.R`](scripts/1_data_import/3_summarize_Mg_proteins.R) | Re-aggregates peptides into protein-level summaries using a multi-tiered tie-breaker |

## Data Availability

- **Proteomics data**: ProteomeXchange / PRIDE, accession [PXD063668](https://www.ebi.ac.uk/pride/archive/projects/PXD063668) (public access).
  To match PRIDE file names to species, seasons, batches, and TMT channels, see the sample key file: [`Supplemental_Table_S12_PRIDE_File_Key_PXD063668.xlsx`](manuscript/Supplemental_Tables/Supplemental_Table_S12_PRIDE_File_Key_PXD063668.xlsx).
- **RNA-seq data (maize)**: NCBI BioProject [PRJNA705456](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA705456)
- **RNA-seq data (Tripsacum)**: NCBI BioProject [PRJNA1260937](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1260937)
- **Genome assemblies**:
  - *T. dactyloides*, *A. gerardii*, *S. nutans* — [MaizeGDB PanAnd](https://maizegdb.org/PanAnd_project)
  - *M. giganteus* (Msinensis_497_v7.0), *P. virgatum* (Pvirgatumvar_WBCHAP1_778_v1.0) — [JGI Data Portal](https://data.jgi.doe.gov/)

## Dependencies

### R Packages
- readxl, dplyr, tidyr, readr, stringr (data manipulation)
- ggplot2, FactoMineR, factoextra (visualization and multivariate analysis)
- ggtree, phytools, ape (phylogenetic analysis)
- topGO (Gene Ontology enrichment)
- DESeq2 (differential gene expression for RNA-seq)

### External Tools
- [seqkit](https://bioinf.shenwei.me/seqkit/) v2.0 — sequence extraction
- [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/) v8.2.13 — phylogenetic tree construction (PROTGAMMAJTT, 100 bootstraps)
- [PANNZER2](http://ekhidna2.biocenter.helsinki.fi/sanspanz/) — functional annotation
- [DIAMOND](https://github.com/bbuchfink/diamond) v2.1.9 — protein sequence comparison
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) v7.475 — multiple sequence alignment
- [AlphaFold3](https://alphafoldserver.com/) — protein structure prediction
- [STAR](https://github.com/alexdobin/STAR) v2.7.10b — RNA-seq read alignment
- [Proteome Discoverer](https://www.thermofisher.com/order/catalog/product/OPTON-31795) 2.5 — MS data processing

## Citation

[Manuscript under review — citation will be updated upon publication]

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
