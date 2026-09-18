import drugstone
import pandas as pd
import json
import time
import os

drugstone.set_api('https://api.stable.drugst.one/')
drugstone.print_license()
drugstone.accept_license()

# ============================================================================
# CONFIGURATION
# ============================================================================
FC_THRESHOLDS   = [4.04, 5.59, 9.32]
STAB_THRESHOLDS = [1.1, 1.2, 1.3]
BASE_DIR = 'valid_DEPs'

COMPARISONS = {
    'diabetic_empty_42-nondiabetic_empty_42': 'Empty defect',
    'diabetic_PCL_42-nondiabetic_PCL_42':     'PCL scaffold',
}

PROTEIN_SETS = [
    'set1_tissue_level',
    'set2_tissue_plus_network',
    'set3_first_level',
]

SET_SCORE_LABELS = {
    'set1_tissue_level':        'Score (tissue-level DEPs)',
    'set2_tissue_plus_network': 'Score (tissue-level DEPs + connectors)',
    'set3_first_level':         'Score (first-level DEPs)',
}

# Drug search parameters
PARAMETERS = {
    "target": "drug",
    "algorithm": "trustrank",
    "ppi_dataset": "nedrex",
    "damping_factor": 0.85,
    "includeIndirectDrugs": False,
    "includeNonApprovedDrugs": True,
    "filterPaths": True, # include only shortest connections in the result
    "resultSize": 100
}


# ============================================================================

def get_genes(comparison, protein_set, base_dir, fc_thresh, stab_thresh):
    """Get genes for a specific protein set"""

    # Load validated DEPs
    validated_file = os.path.join(base_dir, f'FC{fc_thresh}_Stab{stab_thresh}', 'validated_DEPs.csv')
    df = pd.read_csv(validated_file)
    df = df[df['Comparison'] == comparison]

    # Rat-to-human orthology can be ambiguous (a protein may list several
    # candidate orthologs); pick one deterministically per protein so seed
    # genes are reproducible rather than depending on upstream mapping order.
    if protein_set == 'set1_tissue_level':
        # All tissue-level DEPs (first-level + second-level)
        genes = [sorted(gene.split(';'))[0] for gene in df['Orthologs'].dropna().tolist()]

    elif protein_set == 'set2_tissue_plus_network':
        # All tissue-level DEPs + connector proteins (from both_levels network run)
        genes = [sorted(gene.split(';'))[0] for gene in df['Orthologs'].dropna().tolist()]

        exception_file = f'network_enrichment_results/{comparison}/both_levels/exception_proteins.csv'
        if os.path.exists(exception_file):
            exception_df = pd.read_csv(exception_file)
            exception_genes = exception_df['Gene'].dropna().tolist()
            genes.extend(exception_genes)
        else:
            print(f"WARNING: Exception file not found: {exception_file}")

    elif protein_set == 'set3_first_level':
        # First-level DEPs only
        df = df[df['Validation'].str.contains('Exclusive to |(same direction)', regex=True)]
        genes = [sorted(gene.split(';'))[0] for gene in df['Orthologs'].dropna().tolist()]

    else:
        raise ValueError(f"Unknown protein set: {protein_set}")

    genes = list(set(genes))  # Remove duplicates
    return genes


def run_drug_search(genes, comparison, set_name, output_dir):
    """Run drug repurposing search"""

    print(f"\n{'=' * 70}")
    print(f"Running drug search: {comparison} - {set_name}")
    print(f"{'=' * 70}")
    print(f"Number of seed genes: {len(genes)}")

    # Run DrugStone
    task = drugstone.new_task(genes, PARAMETERS)
    result = task.get_result()

    # Save files with unique names
    json_file = os.path.join(output_dir, f'{comparison}_{set_name}_drugs.json')

    # Save directly
    result.download_json()

    # Wait and rename
    time.sleep(2)
    os.rename('result.json', json_file)

    print(f"Saved: {json_file}")

    # Parse results
    with open(json_file, 'r') as f:
        data = json.load(f)

    # Extract drug information
    drug_data = []
    for drug_info in data['drugs'].values():
        drug_record = {
            'Comparison': comparison,
            'Set': set_name,
            'Drug': drug_info['label'],
            'Score': drug_info['score'],
            'DrugBank_ID': drug_info['drugId'],
            'Status': drug_info['status'],
            'Num_Targets': len(drug_info['hasEdgesTo']),
            'Targets': ', '.join(drug_info['hasEdgesTo']),
            'Score_Raw': drug_info['score_raw']
        }
        drug_data.append(drug_record)

    drug_df = pd.DataFrame(drug_data)
    drug_df = drug_df.sort_values(['Score', 'Num_Targets'], ascending=[False, False])

    # Save CSV
    csv_file = os.path.join(output_dir, f'{comparison}_{set_name}_drugs.csv')
    drug_df.to_csv(csv_file, index=False)
    print(f"Saved: {csv_file}")

    # Print summary
    print(f"\nSummary:")
    print(f"  Total drugs: {len(drug_df)}")
    print(f"  Multi-target drugs (3+): {len(drug_df[drug_df['Num_Targets'] >= 3])}")
    print(
        f"  Top drug: {drug_df.iloc[0]['Drug']} (score={drug_df.iloc[0]['Score']:.4f}, targets={drug_df.iloc[0]['Num_Targets']})")

    return drug_df


# ============================================================================
# MAIN ANALYSIS
# ============================================================================

# Collects intersection drugs across all thresholds for the supplement export
_supplement_rows = []

for fc in FC_THRESHOLDS:
    for stab in STAB_THRESHOLDS:
        output_dir = os.path.join(BASE_DIR, f'FC{fc}_Stab{stab}', 'drug_repurposing')
        os.makedirs(output_dir, exist_ok=True)

        for comparison in COMPARISONS:
            set_dfs = {}  # set_name → drug_df

            for set_name in PROTEIN_SETS:
                try:
                    genes = get_genes(comparison, set_name, BASE_DIR, fc, stab)
                    if not genes:
                        print(f"WARNING: No genes found for {comparison} - {set_name}")
                        continue
                    drug_df = run_drug_search(genes, comparison, set_name, output_dir)
                    set_dfs[set_name] = drug_df
                except Exception as e:
                    print(f"ERROR in {comparison} - {set_name}: {e}")
                    continue

            # ── Intersection: drugs present in all 3 protein sets ─────────────
            if len(set_dfs) == 3:
                common_drugs = set(set_dfs[PROTEIN_SETS[0]]['Drug'])
                for s in PROTEIN_SETS[1:]:
                    common_drugs &= set(set_dfs[s]['Drug'])

                rows = []
                for drug in common_drugs:
                    scores = {s: float(set_dfs[s].loc[set_dfs[s]['Drug'] == drug, 'Score'].iloc[0])
                              for s in PROTEIN_SETS}
                    ref_row = set_dfs[PROTEIN_SETS[0]][set_dfs[PROTEIN_SETS[0]]['Drug'] == drug].iloc[0]
                    rows.append({
                        'AR threshold':       fc,
                        'ΔAR':                stab,
                        'Comparison':         comparison,
                        'Drug':               drug,
                        'DrugBank_ID':        ref_row['DrugBank_ID'],
                        'Status':             ref_row['Status'],
                        SET_SCORE_LABELS['set1_tissue_level']:        scores['set1_tissue_level'],
                        SET_SCORE_LABELS['set2_tissue_plus_network']:  scores['set2_tissue_plus_network'],
                        SET_SCORE_LABELS['set3_first_level']:         scores['set3_first_level'],
                        'Mean score':         sum(scores.values()) / 3,
                    })

                intersection_df = (
                    pd.DataFrame(rows)
                    .sort_values('Mean score', ascending=False)
                    .reset_index(drop=True)
                )

                top_file = os.path.join(output_dir, f'{comparison}_top_drugs_intersection.csv')
                intersection_df.to_csv(top_file, index=False)
                print(f"\nIntersection ({comparison}): {len(intersection_df)} drugs in all 3 sets → {top_file}")

                _supplement_rows.append(intersection_df)
            else:
                print(f"WARNING: Only {len(set_dfs)}/3 sets completed for {comparison} FC{fc}_Stab{stab} — skipping intersection")


# ============================================================================
# SUPPLEMENT EXPORT
# ============================================================================

if _supplement_rows:
    all_top = pd.concat(_supplement_rows, ignore_index=True)

    SUPP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'supplements')
    os.makedirs(SUPP_DIR, exist_ok=True)
    OUT_FILE = os.path.join(SUPP_DIR, 'top_drugs_supplement.xlsx')
    README = pd.DataFrame([
        ("Sheet name",  "Content"),
        ("Empty defect", "Diabetic vs. non-diabetic, empty bone defect."),
        ("PCL scaffold", "Diabetic vs. non-diabetic, PCL scaffold."),
        ("",            ""),
        ("Column name", "Description"),
        ("AR threshold",                       "Abundance ratio threshold: minimum comparison-level abundance ratio between fractions required to qualify a protein as a second-level DEP."),
        ("ΔAR",                                "Stability threshold: maximum ratio of group-level geometric mean abundance ratios allowed for a second-level DEP to be considered stable across groups."),
        ("Drug",                               "Drug name."),
        ("DrugBank_ID",                        "DrugBank identifier."),
        ("Status",                             "Approval status (e.g. approved, investigational)."),
        ("Score (tissue-level DEPs)",          "TrustRank score using tissue-level DEPs as seeds."),
        ("Score (tissue-level DEPs + connectors)", "TrustRank score using tissue-level DEPs + network connectors as seeds."),
        ("Score (first-level DEPs)",           "TrustRank score using first-level DEPs only as seeds."),
        ("Mean score",                         "Mean TrustRank score across all three protein sets. Used for ranking."),
        ("",                                   ""),
        ("Note",                               "Only drugs appearing in the top 100 of all three protein sets are included."),
    ], columns=["_", "__"])

    with pd.ExcelWriter(OUT_FILE, engine='openpyxl') as writer:
        README.to_excel(writer, sheet_name='README', index=False, header=False, startrow=2)
        writer.sheets['README']['A1'] = "Drug Repurposing Candidates"
        for comp_key, comp_label in COMPARISONS.items():
            subset = (
                all_top[all_top['Comparison'] == comp_key]
                .drop(columns='Comparison')
                .sort_values(['AR threshold', 'ΔAR', 'Mean score'],
                             ascending=[True, True, False])
                .reset_index(drop=True)
            )
            subset.to_excel(writer, sheet_name=comp_label, index=False)
            print(f"  Sheet '{comp_label}': {len(subset)} rows")

    print(f"\nSaved → {OUT_FILE}")

print("\n✓ Drug repurposing analysis complete!")