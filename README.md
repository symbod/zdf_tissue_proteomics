# Proteomics of diabetic bone healing in ZDF rats

Bone tissue from diabetic and non-diabetic rats with an empty bone defect or a polycaprolactone (PCL) scaffold was split into an acid-insoluble (AI) and an acid-soluble (AS) fraction and measured by mass spectrometry. The pipeline harmonizes and normalizes each fraction, runs differential expression per fraction for the two comparisons (diabetic vs. non-diabetic, with empty defect and with PCL scaffold), consolidates the fraction-specific DEPs into tissue-level DEPs, adds connector proteins from the human PPI network, and runs pathway enrichment, comparison with bone-healing reference sets, and drug repurposing.

## Input (`input/`)

The unprocessed proteomics data is not provided in this repository.

## Setup

- R 4.6 and the R packages at the versions listed in `renv.lock` (Bioconductor 3.23, `proharmed` from GitHub), for example via `renv::restore()`.
- Python 3.11: `pip install pandas==2.3.1 numpy==2.3.1 scipy==1.17.0 matplotlib==3.10.3 matplotlib-venn==1.1.2 openpyxl==3.1.5 drugstone==1.0.2`.
- Steps 1, 5, 6 and 8 query online databases (UniProt and gProfiler via proharmed, Drugst.One, KEGG) and step 7 the GO ontology. Results of these steps can drift as the databases are updated; the numbers reported in the paper are preserved in its supplementary tables.

## Pipeline

1. `Proteomics_Harmonization.Rmd`: filters the protein IDs to rat and maps them to gene names and human orthologs (proharmed).
2. `Proteomics_Normalization.Rmd`: per fraction, RobNorm normalization, scaling by the sham medians and outlier detection (PRONE, POMA).
3. `Proteomics_ExpressionAnalysis.Rmd`: per fraction, differential expression with limma for both comparisons.
4. `tissue_level_DEP_collection.R`: consolidates the fraction-specific DEPs into first-level and second-level tissue-level DEPs for all combinations of the AR and ΔAR thresholds (Supplementary Figure S1).
5. `ppi_enrichment_must.py`: connector proteins that link the tissue-level DEPs in the human PPI network (Drugst.One, Multi-Steiner tree).
6. `pathway_enrichment.R`: GO and KEGG enrichment of the first-level DEPs, the tissue-level DEPs and the tissue-level DEPs with connector proteins (clusterProfiler).
7. Figures and overlap tests:
   - `AI_AS_proteome_investigation.py`: proteins detected and their abundances in AI vs. AS (Figure 4)
   - `AI_AS_dep_per_fraction_investigation.py`: DEPs per fraction and their classification (Figure 5)
   - `protein_sets_overlap_bone_CAPs.py`: overlap of the protein sets with the bone-healing reference proteins (Figure 6, Supplementary Table S4)
   - `scaffold_label_permutation_test.R`: permutation test of the differences between the two comparisons in the number of tissue-level DEPs and in their overlap significance (Section 3.3)
   - `enriched_pathway_overlap_with_bone_pathways.py`: overlap of the enriched terms with the bone-healing reference pathways (Figure 7, Supplementary Table S6)
   - `figure_pathway_hierarchy.py`: top enriched terms per protein set and comparison (Figure 8)
8. `drug_repurposing.py`: drug candidates by TrustRank on the NeDRex network (Drugst.One, Supplementary Table S7).
9. `supplements/export_*.py`: Supplementary Tables S1, S2, S3 and S5.

The main results of the paper use the AR threshold 5.59 (80th percentile of the background abundance-ratio distribution) and ΔAR 1.2.
