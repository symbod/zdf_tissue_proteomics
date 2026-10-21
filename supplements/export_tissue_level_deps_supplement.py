"""
Exports the validated DAPs for every threshold combination and comparison into
a single Excel file, one sheet per (threshold combination × comparison).

Output: tissue_level_deps_supplement.xlsx  — 18 sheets (9 thresholds × 2 comparisons)
Sheet naming: "AR{fc} ΔAR{stab} {comparison label}"
"""

import os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

FC_THRESHOLDS   = [4.04, 5.59, 9.32]
STAB_THRESHOLDS = [1.1, 1.2, 1.3]

COMPARISONS = {
    "diabetic_empty_42-nondiabetic_empty_42": "Empty defect",
    "diabetic_PCL_42-nondiabetic_PCL_42":     "PCL scaffold",
}

OUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "tissue_level_deps_supplement.xlsx")

README = pd.DataFrame([
    ("Sheet name",           "Content"),
    ("AR<value> ΔAR<value> <comparison>", "e.g. AR5.59 ΔAR1.2 Empty defect. One sheet per threshold × comparison (18 total)."),
    ("",                     ""),
    ("Threshold",            "Description"),
    ("AR",                   "Abundance ratio threshold: minimum comparison-level abundance ratio between fractions required to qualify a protein as a second-level DAP."),
    ("ΔAR",                  "Stability threshold: the largest fold by which a change in fraction distribution may contribute to the sham-corrected fraction-level effect of a second-level DAP (see ShareContrast_Fold)."),
    ("",                     ""),
    ("Column name",          "Description"),
    ("Protein.IDs",          "UniProt accessions."),
    ("Gene.Names",           "Gene symbols."),
    ("Orthologs",            "Human ortholog gene symbols."),
    ("Direction",            "Up or Down."),
    ("DAP_in_AI",            "Boolean whether protein is DAP in AI."),
    ("DAP_in_AS",            "Boolean whether protein is DAP in AS."),
    ("AR_between_fractions", "Abundance ratio between fractions."),
    ("Dominant_Fraction",    "Dominant fraction (AI or AS) according to comparison-level abundance ratio."),
    ("Validation",           "Validation reason for the protein being a tissue-level DAP."),
    ("AIShare_Group1_Defect", "Geometric-mean share of the protein in the AI fraction, AI/(AI+AS), across the diabetic defect samples of the comparison."),
    ("AIShare_Group1_Sham",   "Same across the diabetic sham samples."),
    ("AIShare_Group2_Defect", "Same across the non-diabetic defect samples."),
    ("AIShare_Group2_Sham",   "Same across the non-diabetic sham samples."),
    ("ShareContrast_log2",    "Sham-corrected contrast (diabetic defect − diabetic sham) − (non-diabetic defect − non-diabetic sham) of the log2 share in the fraction in which the protein is a DAP: the part of the fraction-level effect attributable to a change in fraction distribution. Empty for DAPs in both fractions."),
    ("ShareContrast_Fold",    "2^|ShareContrast_log2|; a second-level DAP requires this to be at most ΔAR."),
], columns=["_", "__"])

COLUMN_RENAME = {
    "DEP_in_AI":            "DAP_in_AI",
    "DEP_in_AS":            "DAP_in_AS",
    "FC_between_fractions": "AR_between_fractions",
}

with pd.ExcelWriter(OUT_FILE, engine="openpyxl") as writer:
    README.to_excel(writer, sheet_name="README", index=False, header=False, startrow=2)
    writer.sheets["README"]["A1"] = "Validated Differentially Abundant Proteins (DAPs)"
    for fc in FC_THRESHOLDS:
        for stab in STAB_THRESHOLDS:
            path = os.path.join(ROOT, "valid_DEPs", f"FC{fc}_Stab{stab}", "validated_DEPs.csv")
            df   = pd.read_csv(path)

            for comp_key, comp_label in COMPARISONS.items():
                subset     = df[df["Comparison"] == comp_key].drop(columns=["Comparison", "Comparison_Label"]) \
                                                              .rename(columns=COLUMN_RENAME)
                sheet_name = f"AR{fc} ΔAR{stab} {comp_label}"
                subset.to_excel(writer, sheet_name=sheet_name, index=False)
                print(f"  {sheet_name}: {len(subset)} validated DAPs")

print(f"\nSaved → {OUT_FILE}")