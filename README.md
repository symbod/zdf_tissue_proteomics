# Proteomics of diabetic bone healing in ZDF rats

Bone tissue from diabetic and non-diabetic rats with an empty bone defect or a polycaprolactone (PCL) scaffold was split into an acid-insoluble (AI) and an acid-soluble (AS) fraction and measured by mass spectrometry. The pipeline harmonizes and normalizes each fraction, runs differential abundance per fraction for the two comparisons (diabetic vs. non-diabetic, with empty defect and with PCL scaffold), consolidates the fraction-specific DAPs into tissue-level DAPs, adds connector proteins from the human PPI network, and runs pathway enrichment, comparison with bone-regeneration reference sets, and drug repurposing.

## Input (`input/`)

The unprocessed proteomics data is not provided in this repository.

## Setup

- R 4.6 and the R packages at the versions listed in `renv.lock` (Bioconductor 3.23, `proharmed` from GitHub), for example via `renv::restore()`.
- Python 3.11: `pip install pandas==2.3.1 numpy==2.3.1 scipy==1.17.0 matplotlib==3.10.3 matplotlib-venn==1.1.2 openpyxl==3.1.5 drugstone==1.0.2`.
- Steps 1, 5, 6 and 8 query online databases (UniProt and gProfiler via proharmed, Drugst.One, KEGG) and step 7 the GO ontology. Results of these steps can drift as the databases are updated; the numbers reported in the paper are preserved in its supplementary tables.

## Pipeline

1. `Proteomics_Harmonization.Rmd`: filters the protein IDs to rat and maps them to gene names and human orthologs (proharmed).
2. `Proteomics_Normalization.Rmd`: per fraction, RobNorm normalization and outlier detection (PRONE, POMA).
3. `sham_model_differential_abundance.R`: per fraction, differential abundance with limma using all six groups (genotype × sham / empty / PCL). The diabetic effect is the sham-baseline-corrected contrast (diabetic defect − diabetic sham) − (non-diabetic defect − non-diabetic sham); significance at the fraction level is nominal (p < 0.05, |log2FC| ≥ 1). As a robustness check, the script repeats the test on the LoessF-normalized data and reports the agreement of the DAP sets (Jaccard index); only the RobNorm results are used downstream.
4. `tissue_level_DEP_collection.R`: consolidates the fraction-specific DAPs into first-level and second-level tissue-level DAPs for all combinations of the AR and ΔAR thresholds (Supplementary Figure S1). The validation rules are implemented once, in `tissue_level_rules.R`, which this script and both permutation tests use. Since the sham model, the rules look at the samples of the four groups of the sham-corrected contrast (defect and sham animals of both genotypes): exclusivity and dominance are judged over all of them, and the stability criterion of the second level is the sham-corrected contrast applied to the log2 share of the protein in the fraction where it is a DAP, i.e. the part of the fraction-level effect attributable to a change in fraction distribution, which must not exceed the ΔAR threshold (1.2-fold). `tissue_level_permutation_check.R` estimates how many tissue-level DAPs are expected by chance (genotype labels shuffled per animal, same shuffle for both fractions).
5. `ppi_enrichment_must.py`: connector proteins that link the tissue-level DAPs in the human PPI network (Drugst.One, Multi-Steiner tree).
6. `pathway_enrichment.R`: GO and KEGG enrichment of the first-level DAPs, the tissue-level DAPs and the tissue-level DAPs with connector proteins (clusterProfiler, on the human orthologs).
7. Figures and overlap tests:
   - `AI_AS_proteome_investigation.py`: proteins detected and their abundances in AI vs. AS (Figure 4)
   - `AI_AS_dep_per_fraction_investigation.py`: DAPs per fraction and their classification (Figure 5)
   - `protein_sets_overlap_bone_CAPs.py`: overlap of the protein sets with the bone-regeneration reference proteins (Figure 6, Supplementary Table S4)
   - `scaffold_label_permutation_test.R`: permutation test of the difference between the two comparisons in the number of tissue-level DAPs (Section 3.3)
   - `enriched_pathway_overlap_with_bone_pathways.py`: overlap of the enriched terms with the bone-regeneration reference pathways (Figure 7, Supplementary Table S6)
   - `figure_pathway_hierarchy.py`: top enriched terms per protein set and comparison (Figure 8)
8. `drug_repurposing.py`: drug candidates by TrustRank on the NeDRex network (Drugst.One, Supplementary Table S7).
   - `figure_drug_repurposing.py`: scatter plot of all drug candidates (Supplementary Table S7), mean TrustRank score in the empty defect vs. with the PCL scaffold, with the top three drugs of each comparison labelled by name (drugs tied with the third included); drugs with identical scores stack on one marker (`valid_DEPs/FC*/drug_repurposing/drug_candidates_scatter.png`)
9. `supplements/export_*.py`: Supplementary Tables S1, S2, S3 and S5. All supplementary tables are written to `supplements/` under their paper names (`Supplementary_Table_S1.xlsx` to `Supplementary_Table_S7.xlsx`); `supplements/Supplementary_Figure_S1.png` is the captioned version of `valid_DEPs/FC_distribution.png`.

The main results of the paper use the AR threshold 5.59 (80th percentile of the background abundance-ratio distribution) and ΔAR 1.2.
