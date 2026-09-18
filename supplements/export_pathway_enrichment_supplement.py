"""
Exports the pathway enrichment results into a single Excel file.
Both comparisons are combined in each sheet — the 'Comparison' column distinguishes them.

set1 (tis level) and set2 (tis+netwk): 9 sheets each (3 AR × 3 stability), named "AR{fc}-ΔAR{stab} {label}"
set3 (1st level): 1 sheet only — threshold-independent.

Output: pathway_enrichment_supplement.xlsx  — 19 sheets total
"""

import os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

FC_THRESHOLDS   = [4.04, 5.59, 9.32]
STAB_THRESHOLDS = [1.1, 1.2, 1.3]

COMPARISONS = [
    "diabetic_empty_42-nondiabetic_empty_42",
    "diabetic_PCL_42-nondiabetic_PCL_42",
]

OUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pathway_enrichment_supplement.xlsx")

README = pd.DataFrame([
    ("Sheet name",         "Content"),
    ("AR<value>-ΔAR<value> tis level",  "Tissue-level DEPs; one sheet per AR × stability threshold (9 total)."),
    ("AR<value>-ΔAR<value> tis+netwk",  "Tissue-level DEPs + network proteins; one sheet per AR × stability threshold (9 total)."),
    ("1st level",                "First-level DEPs; single sheet, threshold-independent."),
    ("",             ""),
    ("Threshold",           "Description"),
    ("AR",                 "Abundance ratio threshold: minimum comparison-level abundance ratio between fractions required to qualify a protein as a second-level DEP."),
    ("ΔAR",                "Stability threshold: maximum ratio of group-level geometric mean abundance ratios allowed for a second-level DEP to be considered stable across groups."),
    ("",             ""),
    ("Protein set",  "Description"),
    ("tis level",    "Tissue-level DEPs (first- and second-level combined)."),
    ("tis+netwk",    "Tissue-level DEPs + network proteins."),
    ("1st level",    "First-level DEPs only."),
    ("",             ""),
    ("Column name",  "Description"),
    ("ID",           "GO term or KEGG pathway ID."),
    ("Description",  "Pathway or term name."),
    ("GeneRatio",    "Hits / total query genes."),
    ("BgRatio",      "Background ratio."),
    ("RichFactor",   "Rich factor."),
    ("FoldEnrichment","Observed / expected overlap."),
    ("zScore",       "Enrichment z-score."),
    ("pvalue",       "Unadjusted p-value."),
    ("p.adjust",     "BH-adjusted p-value."),
    ("qvalue",       "q-value (FDR)."),
    ("geneID",       "Contributing genes."),
    ("Count",        "Hit count."),
    ("Category",     "BP / MF / CC / KEGG."),
    ("Comparison",   "Comparison identifier."),
    ("Significant",  "TRUE if p.adjust < 0.05."),
], columns=["_", "__"])

THRESHOLD_DEPENDENT_SETS = {
    "set1_tissue_level":        "tis level",
    "set2_tissue_plus_network": "tis+netwk",
}
THRESHOLD_INDEPENDENT_SETS = {
    "set3_first_level": "1st level",
}

with pd.ExcelWriter(OUT_FILE, engine="openpyxl") as writer:
    README.to_excel(writer, sheet_name="README", index=False, header=False, startrow=2)
    writer.sheets["README"]["A1"] = "Pathway Enrichment Analysis Results"

    for fc in FC_THRESHOLDS:
        for stab in STAB_THRESHOLDS:
            for set_folder, set_label in THRESHOLD_DEPENDENT_SETS.items():
                frames = []
                for comp_key in COMPARISONS:
                    path = os.path.join(
                        ROOT, "valid_DEPs", f"FC{fc}_Stab{stab}",
                        "enrichment", set_folder,
                        f"{comp_key}_{set_folder}_enrichment_complete.csv",
                    )
                    if not os.path.exists(path):
                        print(f"  MISSING: {path}")
                        continue
                    frames.append(pd.read_csv(path, index_col=0))
                if not frames:
                    continue
                df = pd.concat(frames, ignore_index=True)
                sheet_name = f"AR{fc}-ΔAR{stab} {set_label}"
                df.to_excel(writer, sheet_name=sheet_name, index=False)
                print(f"  {sheet_name}: {len(df)} rows")

    for set_folder, set_label in THRESHOLD_INDEPENDENT_SETS.items():
        frames = []
        for comp_key in COMPARISONS:
            path = os.path.join(
                ROOT, "valid_DEPs", f"FC{FC_THRESHOLDS[0]}_Stab{STAB_THRESHOLDS[0]}",
                "enrichment", set_folder,
                f"{comp_key}_{set_folder}_enrichment_complete.csv",
            )
            if not os.path.exists(path):
                print(f"  MISSING: {path}")
                continue
            frames.append(pd.read_csv(path, index_col=0))
        if not frames:
            continue
        df = pd.concat(frames, ignore_index=True)
        df.to_excel(writer, sheet_name=set_label, index=False)
        print(f"  {set_label}: {len(df)} rows")

print(f"\nSaved → {OUT_FILE}")