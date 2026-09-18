"""
Exports the DE analysis results (all tested proteins, including non-significant)
with ortholog annotations for AI and AS fractions into a single Excel file,
for use as a supplement.

Output: de_results_supplement.xlsx
  - Sheet "AI DE results" : full DE results for the AI fraction
  - Sheet "AS DE results" : full DE results for the AS fraction
"""

import os
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

AI_DE = os.path.join(ROOT, "de_analysis_results_AI_RobNorm_scaled/de_results_raw_with_orthologs.csv")
AS_DE = os.path.join(ROOT, "de_analysis_results_AS_RobNorm_scaled/de_results_raw_with_orthologs.csv")

OUT_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "de_results_supplement.xlsx")

ai_de = pd.read_csv(AI_DE)
as_de = pd.read_csv(AS_DE)

README = pd.DataFrame([
    ("Sheet name",  "Content"),
    ("AI DE results", "Limma DE results, AI fraction."),
    ("AS DE results", "Limma DE results, AS fraction."),
    ("",            ""),
    ("Column name", "Description"),
    ("Protein.IDs", "UniProt accessions."),
    ("Gene.Names",  "Gene symbols."),
    ("logFC",       "Log2 fold change (diabetic / non-diabetic)."),
    ("P.Value",     "Unadjusted p-value."),
    ("adj.P.Val",   "BH-adjusted p-value."),
    ("Change",      "Up / Down / No Change."),
    ("Comparison",  "Comparison identifier."),
    ("Assay",       "Normalized data the analysis was run on (RobNorm_scaled in all rows); the fraction is given by the sheet name."),
    ("IDs",         "Running row number of the protein in the dataset, assigned by PRONE; no biological meaning."),
    ("Orthologs",   "Human ortholog gene symbols."),
], columns=["_", "__"])

with pd.ExcelWriter(OUT_FILE, engine="openpyxl") as writer:
    README.to_excel(writer, sheet_name="README", index=False, header=False, startrow=2)
    writer.sheets["README"]["A1"] = "Differential Expression Analysis Results"
    ai_de.to_excel(writer, sheet_name="AI DE results", index=False)
    as_de.to_excel(writer, sheet_name="AS DE results", index=False)

print(f"Saved → {OUT_FILE}")
for label, df in [("AI", ai_de), ("AS", as_de)]:
    for comp, group in df.groupby("Comparison"):
        n_sig = (group["Change"] != "No Change").sum()
        print(f"  {label} | {comp}: {len(group)} proteins tested, {n_sig} significant")