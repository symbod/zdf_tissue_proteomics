import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy.stats import hypergeom

# ── Config ────────────────────────────────────────────────────────────────────

FC_THRESHOLDS   = [4.04, 5.59, 9.32]
STAB_THRESHOLDS = [1.1, 1.2, 1.3]
NET_DIR         = "network_enrichment_results"

COMPARISONS = {
    "diabetic_empty_42-nondiabetic_empty_42": "Empty defect",
    "diabetic_PCL_42-nondiabetic_PCL_42":     "PCL scaffold",
}

labels = [
    "First-level\nDEPs",
    "Tissue-level\nDEPs",
    "Tissue-level DEPs +\nconnector proteins",
]

COLORS = {
    "Empty defect": "#4A90C4",
    "PCL scaffold":  "#E07B54",
    "Shared":        "#73B87C",
}
LINE_COLOR = "#1A1A1A"
BAR_WIDTH  = 0.22
GROUP_GAP  = 0.08

# Load data that is shared across all threshold combinations
ai_de     = pd.read_csv("de_analysis_results_AI_RobNorm_scaled/de_results_raw_with_orthologs.csv")
as_de     = pd.read_csv("de_analysis_results_AS_RobNorm_scaled/de_results_raw_with_orthologs.csv")
bone_caps = set(pd.read_csv("input/bone_caps_meta_analysis.csv")["Gene"].dropna().str.upper())


def human_orthologs(ortholog_strings):
    """Set of human gene symbols named by a column of ';'-separated ortholog lists."""
    out = set()
    for s in pd.Series(ortholog_strings).dropna().astype(str):
        out |= {g.strip().upper() for g in s.split(";") if g.strip()}
    return out

# Collect results for supplement export
_results = []

# ── Loop over all threshold combinations ──────────────────────────────────────

for fc in FC_THRESHOLDS:
    for stab in STAB_THRESHOLDS:
        THRESH_DIR = f"valid_DEPs/FC{fc}_Stab{stab}"

        print(f"\n{'='*60}")
        print(f"FC={fc}  Stab={stab}  ({THRESH_DIR})")
        print(f"{'='*60}")

        # ── Compute protein counts and p-values from source files ─────────────

        # Reads:
        #   - summary_statistics.csv  → tissue-level / first-level counts and overlaps
        #   - validated_DEPs.csv      → per-protein gene sets for the "shared" bars
        #   - DE results              → total tested proteome size per comparison
        #   - exception_proteins.csv  → connector proteins from both-levels network enrichment
        #   - bone_caps_meta_analysis → reference gene list for hypergeometric test

        summary = pd.read_csv(os.path.join(THRESH_DIR, "summary_statistics.csv")).set_index("Comparison")

        PROTEIN_N     = {}
        protein_pvals = {}
        valid_deps    = pd.read_csv(os.path.join(THRESH_DIR, "validated_DEPs.csv"))
        _gene_sets    = {}

        for comp_key, comp_label in COMPARISONS.items():
            row = summary.loc[comp_key]

            # Total proteome size = all proteins tested in either fraction for this comparison
            ai_comp  = ai_de[ai_de["Comparison"] == comp_key]
            as_comp  = as_de[as_de["Comparison"] == comp_key]
            ai_prots = set(ai_comp["Protein.IDs"].dropna())
            as_prots = set(as_comp["Protein.IDs"].dropna())
            M = len(ai_prots | as_prots)       # population size
            n = int(row["Bone_Caps_Tested"])   # successes in population

            # Human genes represented by the tested proteome, used below to decide
            # which connector proteins are genuinely new members of the population.
            tested_genes = human_orthologs(
                pd.concat([ai_comp[["Protein.IDs", "Orthologs"]],
                           as_comp[["Protein.IDs", "Orthologs"]]])
                  .drop_duplicates("Protein.IDs")["Orthologs"]
            )

            # Set 1: tissue-level DEPs (first-level + second-level combined)
            s1 = int(row["Validated_Both_Fractions"] + row["Validated_Exclusive"] + row["Validated_Second_Level"])
            o1 = int(row["Overlap_First_Level"] + row["Overlap_Second_Level"])
            # hypergeom.sf(k, M, n, N) with n the number of successes in M and N the number of draws
            p1 = hypergeom.sf(o1 - 1, M, n, s1)

            # Set 2: tissue-level DEPs + connector proteins (from both-levels network).
            # Connectors come from the human interactome, so most of them are not
            # members of the tested proteome. Adding them to the draw without also
            # adding them to the population would make the draw exceed what the
            # population can supply, so the population grows by the connectors that
            # are genuinely new (and its successes by any bone-healing reference
            # protein among them).
            exc_path  = os.path.join(NET_DIR, comp_key, "both_levels", "exception_proteins.csv")
            exc_prots = set(pd.read_csv(exc_path)["Gene"].dropna().str.upper())
            new_conn  = exc_prots - tested_genes
            M2 = M + len(new_conn)
            n2 = n + len((exc_prots & bone_caps) - tested_genes)
            s2 = s1 + len(exc_prots)
            o2 = o1 + len(exc_prots & bone_caps)
            p2 = hypergeom.sf(o2 - 1, M2, n2, s2)

            # Set 3: first-level DEPs only
            s3 = int(row["Validated_Both_Fractions"] + row["Validated_Exclusive"])
            o3 = int(row["Overlap_First_Level"])
            p3 = hypergeom.sf(o3 - 1, M, n, s3)

            # Verify file-derived counts match summary statistics, then store gene sets
            cd         = valid_deps[valid_deps["Comparison"] == comp_key]
            first_mask = ~cd["Validation"].str.startswith("Dominant")
            assert first_mask.sum() == s3, (
                f"{comp_key}: first-level DEP count mismatch — file={first_mask.sum()}, summary={s3}"
            )
            assert len(cd) == s1, (
                f"{comp_key}: tissue-level DEP count mismatch — file={len(cd)}, summary={s1}"
            )
            # Cross-comparison "shared" bars are computed on human gene symbols so
            # that DEPs (measured as rat proteins) and connectors (human interactome
            # genes) live in one namespace; upper-casing the rat symbols instead
            # would match the two only by spelling coincidence.
            first_lvl  = human_orthologs(cd.loc[first_mask, "Orthologs"])
            tissue_lvl = human_orthologs(cd["Orthologs"])
            _gene_sets[comp_key] = [first_lvl, tissue_lvl, tissue_lvl | exc_prots]

            PROTEIN_N[comp_label]     = [s3, s1, s2]
            protein_pvals[comp_label] = [p3, p1, p2]

            for set_label, pop, succ, s, o, p in [
                ("Tissue-level DEPs",                     M,  n,  s1, o1, p1),
                ("Tissue-level DEPs + connector proteins", M2, n2, s2, o2, p2),
                ("First-level DEPs",                      M,  n,  s3, o3, p3),
            ]:
                expected = round((s * succ) / pop, 2) if pop > 0 else 0
                _results.append({
                    "AR threshold":             fc,
                    "ΔAR":                      stab,
                    "Comparison":               comp_label,
                    "Protein set":              set_label,
                    "Proteome size (M)":        pop,
                    "Bone-healing reference proteins in proteome (n)": succ,
                    "Protein set size (N)":     s,
                    "Bone-healing reference protein overlap (k)":    o,
                    "Expected overlap":         expected,
                    "p-value":                  p,
                })

            print(f"\n{comp_label}:  population M={M} (successes n={n})")
            print(f"  Tissue-level DEPs                  N={s1:>4}  overlap={o1:>3}  p={p1:.3e}")
            print(f"  Tissue-level DEPs + connectors     N={s2:>4}  overlap={o2:>3}  p={p2:.3e}"
                  f"   [population expanded to M={M2} by {len(new_conn)} new connectors]")
            print(f"  First-level DEPs only              N={s3:>4}  overlap={o3:>3}  p={p3:.3e}")

        # ── Shared proteins (intersection between comparisons) ────────────────

        comp_keys = list(COMPARISONS.keys())
        PROTEIN_N["Shared"] = [
            len(_gene_sets[comp_keys[0]][i] & _gene_sets[comp_keys[1]][i])
            for i in range(3)
        ]

        # ── Figure ────────────────────────────────────────────────────────────

        x          = np.arange(len(labels))
        bar_comps  = list(PROTEIN_N.keys())       # Empty defect, PCL scaffold, Shared
        pval_comps = list(protein_pvals.keys())   # Empty defect, PCL scaffold
        n_bar      = len(bar_comps)
        n_pval     = len(pval_comps)

        fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(10.5, 4.5))
        fig.subplots_adjust(left=0.08, right=0.98, top=0.95, bottom=0.30, wspace=0.32)

        # ── Panel (a): protein set sizes ──────────────────────────────────────

        for ci, comp in enumerate(bar_comps):
            offset    = (ci - (n_bar - 1) / 2) * (BAR_WIDTH + GROUP_GAP / 2)
            positions = x + offset
            ax_a.bar(positions, PROTEIN_N[comp], width=BAR_WIDTH,
                     color=COLORS[comp], alpha=0.90,
                     edgecolor="white", linewidth=0.6, zorder=2)
            for pos, n in zip(positions, PROTEIN_N[comp]):
                ax_a.text(pos, n + 3, f"n={n}", ha="center", va="bottom",
                          fontsize=7.5, color=COLORS[comp], fontweight="bold")

        ymax_a = max(v for vals in PROTEIN_N.values() for v in vals) * 1.28
        ax_a.set_ylabel("Number of proteins", fontsize=10, labelpad=6)
        ax_a.set_ylim(0, ymax_a)
        ax_a.set_xlim(-0.6, len(labels) - 0.4)
        ax_a.tick_params(axis="y", labelsize=9)
        ax_a.yaxis.grid(True, linestyle=":", linewidth=0.7, alpha=0.5, zorder=0)
        ax_a.set_axisbelow(True)
        ax_a.spines[["top", "right"]].set_visible(False)
        ax_a.set_xticks(x)
        ax_a.set_xticklabels(labels, fontsize=9)
        ax_a.set_xlabel("Protein set", fontsize=10, labelpad=14)
        ax_a.text(-0.01, 1.02, "(a)", transform=ax_a.transAxes,
                  fontsize=11, fontweight="bold", va="bottom", ha="left")

        legend_handles = [
            mpatches.Patch(color=COLORS["Empty defect"], label="Empty defect"),
            mpatches.Patch(color=COLORS["PCL scaffold"],  label="PCL scaffold"),
            mpatches.Patch(color=COLORS["Shared"],        label="Shared (both comparisons)"),
        ]

        # ── Panel (b): overlap significance — lollipop chart ─────────────────

        all_neg_log = [-np.log10(p) for pvals in protein_pvals.values() for p in pvals]
        ymax_b      = max(all_neg_log) * 1.22

        for ci, comp in enumerate(pval_comps):
            neg_log   = [-np.log10(p) for p in protein_pvals[comp]]
            offset    = (ci - (n_pval - 1) / 2) * (BAR_WIDTH + GROUP_GAP / 2)
            positions = x + offset
            for pos, val in zip(positions, neg_log):
                val_clipped = min(val, ymax_b)
                ax_b.scatter(pos, val_clipped, color=COLORS[comp],
                             s=100, zorder=3, edgecolors="white", linewidths=0.6)
                if val > ymax_b:
                    ax_b.annotate(
                        "", xy=(pos, ymax_b - 0.02), xytext=(pos, ymax_b - 0.4),
                        arrowprops=dict(arrowstyle="-|>", color=COLORS[comp], lw=2.0),
                        zorder=6,
                    )

        sig_line = -np.log10(0.05)
        ax_b.axhline(sig_line, color=LINE_COLOR, linewidth=1.8, linestyle="--", zorder=5)
        ax_b.text(-0.55, sig_line + ymax_b * 0.02, "p = 0.05",
                  ha="left", va="bottom", fontsize=8, color=LINE_COLOR)

        ax_b.set_ylabel("−log₁₀(p-value)\n(overlap with bone-healing reference proteins)", fontsize=10, labelpad=6)
        ax_b.set_ylim(0, ymax_b)
        ax_b.set_xlim(-0.6, len(labels) - 0.4)
        ax_b.tick_params(axis="y", labelsize=9)
        ax_b.yaxis.grid(True, linestyle=":", linewidth=0.7, alpha=0.5, zorder=0)
        ax_b.set_axisbelow(True)
        ax_b.spines[["top", "right"]].set_visible(False)
        ax_b.set_xticks(x)
        ax_b.set_xticklabels(labels, fontsize=9)
        ax_b.set_xlabel("Protein set", fontsize=10, labelpad=14)
        ax_b.text(-0.01, 1.02, "(b)", transform=ax_b.transAxes,
                  fontsize=11, fontweight="bold", va="bottom", ha="left")

        # ── Shared legend centered below both panels ──────────────────────────

        fig.legend(handles=legend_handles, fontsize=9, frameon=False,
                   loc="lower center", bbox_to_anchor=(0.53, 0.01), ncol=3,
                   handlelength=1.2, handletextpad=0.5, columnspacing=1.2,
                   title="Comparison",
                   title_fontproperties={"weight": "bold", "size": 10})

        # ── Save ──────────────────────────────────────────────────────────────

        out_path = os.path.join(THRESH_DIR, "validation_protein_overlap.png")
        plt.savefig(out_path, dpi=180, bbox_inches="tight")
        plt.close()
        print(f"Saved {out_path}")

# ── Supplement export ─────────────────────────────────────────────────────────

PROTEIN_SET_ORDER = [
    "Tissue-level DEPs",
    "Tissue-level DEPs + connector proteins",
    "First-level DEPs",
]

results_df = pd.DataFrame(_results)

SUPP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "supplements")
os.makedirs(SUPP_DIR, exist_ok=True)
OUT_FILE = os.path.join(SUPP_DIR, "protein_bone_overlap_supplement.xlsx")
README = pd.DataFrame([
    ("Sheet name",   "Content"),
    ("Empty defect", "Diabetic vs. non-diabetic, empty bone defect."),
    ("PCL scaffold", "Diabetic vs. non-diabetic, PCL scaffold."),
    ("",             ""),
    ("Column name",  "Description"),
    ("AR threshold",              "Abundance ratio threshold: minimum comparison-level abundance ratio between fractions required to qualify a protein as a second-level DEP."),
    ("ΔAR",                       "Stability threshold: minimum ratio of group-level geometric mean abundance ratios required for a second-level DEP to be considered stable across groups."),
    ("Protein set size (N) | Tissue-level DEPs",                    "Number of tissue-level DEPs."),
    ("Protein set size (N) | Tissue-level DEPs + connector proteins","Number of tissue-level DEPs + network connector proteins."),
    ("Protein set size (N) | First-level DEPs",                     "Number of first-level DEPs."),
    ("Bone-healing reference protein overlap (k) | Tissue-level DEPs",                    "Tissue-level DEPs also in the bone-healing reference proteins."),
    ("Bone-healing reference protein overlap (k) | Tissue-level DEPs + connector proteins","Tissue-level DEPs + connectors also in the bone-healing reference proteins."),
    ("Bone-healing reference protein overlap (k) | First-level DEPs",                     "First-level DEPs also in the bone-healing reference proteins."),
    ("Expected overlap | Tissue-level DEPs",                    "Expected overlap under hypergeometric null, tissue-level DEPs."),
    ("Expected overlap | Tissue-level DEPs + connector proteins","Expected overlap under hypergeometric null, tissue-level DEPs + connectors."),
    ("Expected overlap | First-level DEPs",                     "Expected overlap under hypergeometric null, first-level DEPs."),
    ("p-value | Tissue-level DEPs",                    "Hypergeometric p-value, tissue-level DEPs."),
    ("p-value | Tissue-level DEPs + connector proteins","Hypergeometric p-value, tissue-level DEPs + connectors."),
    ("p-value | First-level DEPs",                     "Hypergeometric p-value, first-level DEPs."),
], columns=["_", "__"])

with pd.ExcelWriter(OUT_FILE, engine="openpyxl") as writer:
    README.to_excel(writer, sheet_name="README", index=False, header=False, startrow=2)
    writer.sheets["README"]["A1"] = "Protein Set Overlap with Bone-Healing Reference Proteins"
    for comp_label in COMPARISONS.values():
        subset = results_df[results_df["Comparison"] == comp_label].drop(columns="Comparison")
        pivot  = subset.pivot(
            index=["AR threshold", "ΔAR"],
            columns="Protein set",
            values=["Protein set size (N)", "Bone-healing reference protein overlap (k)", "Expected overlap", "p-value"],
        )
        pivot = pivot.reindex(columns=pd.MultiIndex.from_product([
            ["Protein set size (N)", "Bone-healing reference protein overlap (k)", "Expected overlap", "p-value"],
            PROTEIN_SET_ORDER,
        ]))
        pivot.columns = [f"{metric} | {pset}" for metric, pset in pivot.columns]
        pivot = pivot.reset_index()

        pivot.to_excel(writer, sheet_name=comp_label, index=False)

print(f"Saved → {OUT_FILE}")
print(f"  Sheets: {list(COMPARISONS.values())}")