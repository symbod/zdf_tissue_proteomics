import drugstone
import pandas as pd
import json
import time
import os
import glob

drugstone.set_api('https://api.stable.drugst.one/')
drugstone.print_license()
drugstone.accept_license()

VALIDATED_DEPS = 'valid_DEPs/FC5.59_Stab1.2/validated_DEPs.csv'

COMPARISONS = [
    "diabetic_empty_42-nondiabetic_empty_42",
    "diabetic_PCL_42-nondiabetic_PCL_42",
]

# First-level DAPs: validated by fraction agreement or exclusivity
FIRST_LEVEL = [
    "Significant in both fractions (same direction)",
    "Exclusive to AI fraction",
    "Exclusive to AS fraction",
]

# Two filter sets: first-level only, or all validated (first + second level)
FILTER_SETS = {
    "first_level":   FIRST_LEVEL,
    "both_levels": None,   # None → keep all rows where Validation is not NaN
}

ENRICHMENT_PARAMETERS = {
    "target":      "drug-target",
    "algorithm":   "multisteiner",
    "ppiDataset":  "NeDRex",
    "num_trees":   5,
    "tolerance":   5,
    "hubPenalty":  0.5,
    "resultSize":  500,
}

df       = pd.read_csv(VALIDATED_DEPS)
root_dir = os.getcwd()

for comparison in COMPARISONS:
    for filter_name, validation_filters in FILTER_SETS.items():

        print(f"\n{'='*60}")
        print(f"Comparison : {comparison}")
        print(f"Filter set : {filter_name}")
        print(f"{'='*60}")

        comp_df     = df[df["Comparison"] == comparison]
        filtered_df = (
            comp_df[comp_df["Validation"].isin(validation_filters)]
            if validation_filters is not None
            else comp_df[comp_df["Validation"].notna()]
        )

        # Rat-to-human orthology can be ambiguous (a protein may list several
        # candidate orthologs); pick one deterministically per protein so seed
        # genes are reproducible rather than depending on upstream mapping order.
        my_genes = list(set(
            sorted(gene.split(";"))[0]
            for gene in filtered_df["Orthologs"].dropna()
        ))

        print(f"Proteins after filtering : {len(filtered_df)}")
        print(f"Seed proteins            : {len(my_genes)}")

        output_dir = os.path.join("network_enrichment_results", comparison, filter_name)
        os.makedirs(output_dir, exist_ok=True)

        print("\n=== Performing Network Enrichment ===")
        enrichment_task   = drugstone.new_task(my_genes, ENRICHMENT_PARAMETERS)
        enrichment_result = enrichment_task.get_result()

        # download_json/download_graph write into the current working directory and never overwrite:
        # with a result.json from a previous run present they would write result(1).json, which is
        # not read below. Remove the files of the previous run first.
        for stale in glob.glob(os.path.join(output_dir, "result*.json")) + glob.glob(os.path.join(output_dir, "graph*.html")):
            os.remove(stale)
        os.chdir(output_dir)
        enrichment_result.download_json()
        enrichment_result.download_graph()
        os.chdir(root_dir)

        print("Waiting for enrichment files to be written...")
        time.sleep(3)

        with open(os.path.join(output_dir, "result.json"), "r") as f:
            enriched_data = json.load(f)

        print("\n=== Network Enrichment Results ===")

        if "genes" not in enriched_data:
            print("No genes found in result.")
            continue

        all_nodes           = enriched_data["genes"]
        seed_lower          = {g.lower() for g in my_genes}
        exception_proteins  = [
            {"Gene": info.get("label", nid), "Node_ID": nid, "Type": "Exception"}
            for nid, info in all_nodes.items()
            if info.get("label", nid).lower() not in seed_lower
        ]

        print(f"Total proteins in enriched network : {len(all_nodes)}")
        print(f"Seed proteins (input)              : {len(my_genes)}")
        print(f"Exception proteins (added)         : {len(exception_proteins)}")

        pd.DataFrame(exception_proteins).to_csv(
            os.path.join(output_dir, "exception_proteins.csv"), index=False
        )

        summary = {
            "Total_Nodes":        len(all_nodes),
            "Seed_Proteins":      len(my_genes),
            "Exception_Proteins": len(exception_proteins),
            "Total_Edges":        len(enriched_data.get("edges", {})),
        }
        with open(os.path.join(output_dir, "network_summary.txt"), "w") as f:
            f.write("=== Network Summary ===\n")
            for key, value in summary.items():
                f.write(f"{key.replace('_', ' ')}: {value}\n")

        print(f"Results saved to {output_dir}/")