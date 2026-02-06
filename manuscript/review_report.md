# Peer Review Report

**Manuscript:** Constrained evolution of a core winter proteome across independently cold-adapted PACMAD grasses

**Journal:** *New Phytologist*

**Reviewer Expertise:** Plant evolutionary physiology, comparative genomics, quantitative proteomics (Poaceae)

---

## Recommendation: Minor Revision

---

## Summary of Research

This study uses comparative TMT-based proteomics of rhizomes from five independently cold-adapted PACMAD grass species sampled under natural winter and summer conditions to test whether protein-level cold responses are more evolutionarily constrained than previously observed at the transcript level. The authors report that shared cold-responsive orthogroups show substantially higher cross-species fold-change correlation (ρ = 0.80) than background proteins (ρ = 0.45), and identify LEA3 as the only ortholog elevated in all five species. A structural analysis of the maize LEA3 ortholog reveals a 15-amino-acid insertion that disrupts the conserved 11-mer motif architecture, offering a candidate mechanism for why transcriptional induction alone does not confer freezing tolerance in maize. The work extends the Cope et al. (2025) framework—that proteins face stronger evolutionary constraints than transcripts—into an adaptive stress context, representing a meaningful conceptual advance.

---

## Major Concerns

### 1. Unbalanced experimental design and statistical power asymmetry

The Tripsacum dataset dominates the study (n ≈ 32 samples across two years) while the other four species have n = 2–3 biological replicates per season. The authors acknowledge this (p. 10–11) and argue that the intersection-based comparative analysis is constrained by the weakest dataset. However, this creates two problems that are not adequately addressed:

**(a)** The DAP lists for small-n species are almost certainly incomplete, meaning the "core winter proteome" is biased toward very large-effect proteins. The authors should explicitly model the statistical power of their small-n species (e.g., what fold-change magnitude is detectable at n = 3 with their variance structure?) and discuss how many conserved DAPs they might be missing. The current statement that this represents "a high-confidence set rather than an exhaustive catalog" is correct but insufficient—a power analysis or simulation would substantially strengthen the claim.

**(b)** *Sorghastrum nutans* has only n = 2 biological replicates per season due to winter mortality. At n = 2, formal statistical testing is severely underpowered. The authors should clarify whether S. nutans was included in the Wilcoxon rank-sum test for shared DAPs vs. background (Fig. 2B), and if so, whether its inclusion or exclusion changes the result. A sensitivity analysis dropping S. nutans would be informative.

### 2. Common-garden design confounds species with genotype

All five species were grown at two sites in Ithaca, NY (USDA Zone 6a), ensuring identical winter conditions—a clear strength. However, these species naturally occur across vastly different climatic ranges. The observed proteomic responses reflect each species' reaction to *this specific* winter, not necessarily the response at the edge of each species' cold tolerance. For example, *Miscanthus × giganteus* may be well within its comfort zone at −16°C, while another species may be near its lethal threshold. The authors should discuss whether the high fold-change correlation might partly reflect the fact that all species experienced the same absolute stressor rather than being pushed to equivalent physiological limits. This does not invalidate the finding but nuances its interpretation.

### 3. The parallel vs. convergent evolution argument needs sharper framing

The Discussion distinguishes "parallel evolution" (same orthologs, e.g., LEA3) from "convergent evolution" (different HSP family members achieving similar function). This is a valuable distinction, but the manuscript does not rigorously test it. Specifically:

- The claim of "parallel recruitment" of LEA3 assumes that LEA3 cold-responsiveness arose independently in each lineage rather than being ancestral. If the common ancestor of these PACMAD species already cold-induced LEA3 (even at lower latitudes, e.g., in tropical montane environments), this would represent *retention of an ancestral response* rather than parallel evolution. The authors touch on this (p. 34: "the ancestral protein was already suited to cryoprotection"), but this contradicts the "parallel evolution" framing used throughout. The distinction between "parallel independent recruitment" and "shared ancestral retention under continued constraint" is critical for the evolutionary narrative and should be resolved explicitly.

- The convergent HSP argument would be strengthened by ancestral state reconstruction or at minimum a more formal mapping of which HSP subfamilies are cold-responsive in each lineage onto the species phylogeny.

### 4. Cross-clade comparison with Pooideae transcript data is problematic

The comparison between PACMAD protein-level correlations (ρ = 0.80) and Pooideae transcript-level correlations (r = 0.09–0.40; Schubert et al., 2019) is central to the paper's narrative but compares across multiple confounded variables simultaneously: clade (PACMAD vs. Pooideae), molecular level (protein vs. transcript), tissue (rhizome vs. seedling), stress regime (natural seasonal vs. controlled cold treatment), and developmental stage. The authors acknowledge this (p. 32: "although this comparison spans different clades, tissues, and experimental designs") but then proceed to draw strong conclusions from it regardless.

To strengthen this comparison, the authors should either:
- Generate transcript data from the same rhizome samples (even for a subset of species) to enable a within-study protein-vs-transcript comparison, or
- More prominently caveat that the protein > transcript constraint conclusion relies on the Cope et al. framework as indirect evidence, and that their cross-clade comparison is illustrative rather than definitive.

The existing Tripsacum seedling RNA-seq comparison (ρ = 0.29 for leaf mRNA vs. rhizome protein; Fig. S5A) actually provides *some* within-system evidence, but the tissue difference makes it hard to interpret. This limitation should be discussed more transparently.

### 5. LEA3 structural analysis: compelling but over-interpreted

The LEA3 story is the manuscript's strongest element, but several claims extend beyond the supporting evidence:

**(a)** The hydropathy comparison between Tripsacum and maize LEA3 is descriptive. The claim that the 15-aa insertion "disrupts" function is a hypothesis, not a demonstrated result. The language should consistently reflect this (e.g., "may impair" rather than "disrupts").

**(b)** The statement that "mean hydropathy at Dure's hydrophobic positions averaged approximately −1.1 in Tripsacum, whereas the maize ortholog showed values 15–20% more hydrophilic" (p. 30) requires a formal statistical test or at least confidence intervals, not just a percentage difference.

**(c)** The AlphaFold3 predictions are mentioned but the multimer prediction results (Fig. S8) are referred to only in passing. If the authors are going to invoke structural predictions, the limitations of AlphaFold3 for intrinsically disordered proteins (IDPs) should be explicitly stated. LEA3 proteins are IDPs that adopt structure only upon binding; AlphaFold3's confidence for such proteins is generally low. Were pLDDT scores reported?

**(d)** The comparison to wheat TdLEA3 and Arabidopsis AtLEA7 (Fig. S6) is valuable but these are from different LEA subfamilies (LEA_4 for TdLEA3, LEA_4 for AtLEA7). The authors should clarify whether these outgroup comparisons are within the same orthogroup or represent functional analogy.

### 6. TMT-based quantification: batch effects and normalization transparency

The study uses TMTpro labeling across multiple batches (TMT16 for 2019, TMT18 for 2022, seven separate batches for 2022). While the authors describe batch correction for Tripsacum PCA (linear model residuals), several points require clarification:

- How were cross-batch comparisons handled for the four non-Tripsacum species analyzed in three separate TMT18 batches? Were species confounded with batches? If *Andropogon* and *Miscanthus* were in the same batch (as implied by the shared search against the *A. gerardii* proteome, p. 11), were winter and summer samples of each species within the same TMT plex, or split across plexes? This is critical for accurate fold-change estimation.

- The "total peptide amount" normalization is a basic approach. Given the multi-species, multi-batch design, was any internal reference channel or bridge sample used to enable cross-batch normalization? If not, the authors should justify why simple normalization is adequate.

- The DIAMOND BLASTp re-mapping of *Miscanthus* peptides from the *A. gerardii* proteome to *M. sinensis* (p. 11) is a non-standard procedure. How many proteins were gained or lost in this re-mapping? What was the false reassignment rate? This should be validated more thoroughly (e.g., by showing that fold-changes are robust to the choice of reference proteome).

---

## Minor Concerns

### Nomenclature and formatting

1. The species name "*Miscanthus × giganteus*" is used throughout, but the proteome searched is *M. sinensis* (v7.0). This discrepancy should be noted more prominently, as *M. × giganteus* is a sterile triploid hybrid between *M. sinensis* and *M. sacchariflorus*, and peptide assignment to a diploid reference may systematically miss alleles from the *M. sacchariflorus* subgenome.

2. Table 1: "*Angropogon gerardii*" — typo, should be "*Andropogon gerardii*."

3. The title uses "PACMAD grasses" but only five of the six PACMAD subfamilies are represented (no Arundinoideae). The title is technically accurate (these are PACMAD grasses) but the reader might expect broader subfamily coverage. Consider "Andropogoneae and allied PACMAD grasses" or similar, or clarify in the Introduction that the sampling covers four subfamilies.

4. Figure 1E: The phylogeny includes *Z. mays* for reference but does not include branch lengths or a scale bar beyond "5 Mya." Clarify the source of this tree.

### Statistical and analytical points

5. The threshold of |log₂FC| ≥ 1 and adjusted P < 0.05 is described as "stringent" (p. 10), but a 2-fold change cutoff is actually quite standard. The rationale for this threshold is well-argued (addressing power imbalance), but calling it "stringent" is misleading—"conservative" or "robust" would be more accurate.

6. The chi-square test on functional categories (χ²[40] = 58.1, P = 0.032) has a large number of cells. Were expected frequencies ≥ 5 in all cells? If not, a Fisher's exact test or permutation test would be more appropriate.

7. For the cross-species correlation analysis (Fig. 2B), the comparison uses the Wilcoxon rank-sum test on pairwise Spearman correlations. However, pairwise correlations from the same set of species are not independent. A permutation-based test or bootstrap procedure would better account for this non-independence.

8. The OrthoFinder analysis identified 51,442 orthogroups across five proteomes, of which 4,888 contained quantified proteins. How many orthogroups were species-specific singletons, and how does this affect the comparative analysis? The filtering to orthogroups detected in ≥4 of 5 species for the correlation analysis should be justified more explicitly.

### Figures

9. Figure 2A: The scatterplots are dense and the labeled proteins (LEA3, GST, HSP1, PEBP) are difficult to distinguish. Consider using larger point sizes for labeled DAPs or inset zoom panels for the most conserved orthogroups.

10. Figure 3B: The table of top 50 proteins per species is informative but dense. The color-coding by functional category is helpful, but the small font makes it difficult to read. Consider showing only the top 20 in the main figure and moving the full table to supplementary material.

11. Figure 4C (HSP tree): The tree is very large and the heat map bars are small. The main point—that different HSP clades are recruited in different species—could be made more effectively with a simplified summary figure showing which HSP subfamilies are cold-responsive in each species.

12. Figure 5: The hydropathy plots (panel C) are clear, but the structural models (panel B) would benefit from a side-by-side comparison at the same scale and orientation, with the insertion region in maize explicitly highlighted.

### Writing and presentation

13. The Introduction is well-structured but slightly long (621 words). The paragraph on CBF/DREB signaling and species-specific adaptations (p. 4–5) could be condensed, as the paper does not deeply engage with signaling pathways.

14. The Discussion's final paragraph on "future work" (p. 36) mentioning transgenic complementation or CRISPR editing of ZmLEA3 is appropriate but could note the practical challenges (e.g., that freezing tolerance is polygenic, and single-gene edits are unlikely to produce dramatic phenotypic shifts in the field).

15. P. 6, "interspecific and intraspecific hybrids of *Tripsacum dactyloides*": Clarify earlier that the Tripsacum material includes both inter- and intraspecific crosses, as this is potentially confusing when the species is later referred to simply as "Tripsacum."

16. The Data Availability section is commendably thorough.

17. P. 29: "maximum‐likelihood and maximum-likelihood trees" appears to be a duplication error. Was one intended to be "neighbor-joining" or "Bayesian"?

18. Supporting Information listing mentions "Fig. S8 Phylogenetic analysis of 14-3-3 proteins" but 14-3-3 proteins are not discussed in the main text. Either briefly mention this analysis or remove the reference.

---

## Summary Assessment

This manuscript presents a well-conceived comparative proteomics study that addresses a genuine gap in our understanding of how freezing tolerance evolves in warm-season grasses. The central finding—that protein-level response magnitudes are more conserved across independently cold-adapted species than expected from transcript-level studies—is novel and aligns with recent macroevolutionary theory (Cope et al., 2025). The LEA3 structural story provides a compelling mechanistic anchor. These strengths place the work above a purely descriptive proteomics survey.

However, the study's impact is limited by: (1) the power asymmetry across species that makes the "core proteome" likely incomplete, (2) the confounded comparison with Pooideae transcript data, and (3) some over-interpretation of the LEA3 structural analysis. Addressing these concerns—primarily through additional statistical analyses, sensitivity tests, and more careful language—should be achievable without new experiments. The work is suitable for *New Phytologist* after revision.
