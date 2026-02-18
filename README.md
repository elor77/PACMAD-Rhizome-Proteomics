# PACMAD-Rhizome-Proteomics

Code and data repository for:

**"Constrained evolution of a core winter proteome across independently cold-adapted PACMAD grasses"**

Oren E, Zhai J, Rooney TE, Angelovici R, Hale CO, Brindisi LJ, Hsu S-K, Gault CM, Hua J, La T, Lepak N, Fu Q, Buckler ES, Romay MC

*Corresponding authors:* Elad Oren (eladoren@volcani.agri.gov.il), M. Cinta Romay (mcr72@cornell.edu)

## Directory Structure

Scripts contain hardcoded paths from the original development environment; adjust as needed.

```
├── data/                    # Raw proteomics and RNA-seq data files
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
├── tables/                  # Result tables (DAPs, GO enrichment, log2FC, orthogroups)
├── annotations/             # PANNZER GO term assignments for DAPs
├── sequences/               # Aligned and unaligned protein sequences (LEA, HSP, EF1)
├── uniprot/                 # DIAMOND BLASTp results against SwissProt
├── phylogeny/               # RAxML best trees for EF1, LEA, and HSP families
├── LEA3_structure/          # AlphaFold3 predictions, AmphipaSeek, and ProtScale outputs
└── manuscript/              # Manuscript files and figures
```

### Supplemental Code S1–S3 (Protein Reassignment Pipeline for *M. giganteus*)

| ID | Script | Description |
|:--:|:---|:---|
| S1 | [`1_import_Ag_Mg_data.R`](scripts/1_data_import/1_import_Ag_Mg_data.R) | Imports original Ag/Mg multiplexed peptide intensities and metadata |
| S2 | [`2_map_Mg_peptides_to_proteins.R`](scripts/1_data_import/2_map_Mg_peptides_to_proteins.R) | Re-aligns *M. giganteus* peptides to the *M. giganteus* proteome using DIAMOND |
| S3 | [`3_summarize_Mg_proteins.R`](scripts/1_data_import/3_summarize_Mg_proteins.R) | Re-aggregates peptides into protein-level summaries using a multi-tiered tie-breaker |

## Data Availability

- **Proteomics**: ProteomeXchange / PRIDE, accession [PXD063668](https://www.ebi.ac.uk/pride/archive/projects/PXD063668). See [`Supplemental_Table_S12_PRIDE_File_Key_PXD063668.xlsx`](manuscript/Supplemental_Tables/Supplemental_Table_S13_PRIDE_File_Key_PXD063668.xlsx) for sample-to-file mapping.
- **RNA-seq (maize)**: NCBI BioProject [PRJNA705456](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA705456)
- **RNA-seq (Tripsacum)**: NCBI BioProject [PRJNA1260937](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1260937)
- **Genomes**: *T. dactyloides*, *A. gerardii*, *S. nutans* — [MaizeGDB PanAnd](https://maizegdb.org/PanAnd_project); *M. giganteus*, *P. virgatum* — [JGI Data Portal](https://data.jgi.doe.gov/)

## License

MIT License — see [LICENSE](LICENSE) for details.
