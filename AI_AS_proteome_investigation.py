"""
Single combined figure comparing proteins detected in AI vs AS fractions
across all non-sham samples. Left: 2×2 grid of Venn diagrams, one per
experimental group. Right: abundance scatter plot with one point per protein
per group, coloured by group. Sham samples (scaffold == "sham") are excluded.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib_venn import venn2
from matplotlib.transforms import offset_copy
from scipy.stats import pearsonr

OUT_DIR = "ai_as_proteome_investigation_results"
os.makedirs(OUT_DIR, exist_ok=True)

AI_NORM = "normalization_results_AI/RobNorm_normalized_data.csv"
AS_NORM = "normalization_results_AS/RobNorm_normalized_data.csv"
META    = "input/SampleDescription.csv"

N_META = 9

# Fraction colours: blue / vermillion
COLOR_AI = "#0072B2"
COLOR_AS = "#D55E00"
ALPHA = 0.6

GROUPS = [
    ("diabetic",    "empty", "Diabetic empty bone defect"),
    ("nondiabetic", "empty", "Non-diabetic empty bone defect"),
    ("diabetic",    "PCL",   "Diabetic PCL scaffold"),
    ("nondiabetic", "PCL",   "Non-diabetic PCL scaffold"),
]
# Group colours: brown, magenta, green, purple — none overlap with blue/vermillion
GROUP_COLORS = ["#8c564b", "#e377c2", "#2ca02c", "#9467bd"]

ai_df = pd.read_csv(AI_NORM)
as_df = pd.read_csv(AS_NORM)
meta  = pd.read_csv(META, sep=";")

id_to_gene = pd.concat([
    ai_df[["Protein.IDs", "Gene.Names"]],
    as_df[["Protein.IDs", "Gene.Names"]],
]).drop_duplicates("Protein.IDs").set_index("Protein.IDs")["Gene.Names"]


def samples_for(fraction, genotype, scaffold):
    mask = (
        (meta["fraction"] == fraction) &
        (meta["genotype"] == genotype) &
        (meta["scaffold"] == scaffold)
    )
    return set(meta.loc[mask, "sample_name"])


def detected_proteins(df, sample_names):
    cols = [c for c in df.columns[N_META:] if c in sample_names]
    return set(df.loc[df[cols].notna().any(axis=1), "Protein.IDs"])


def label_venn_with_percentages(v, counts, ax, fontsize=13, pct_fontsize=11.5):
    """Label each venn2 subset with its size plus, on a smaller second line, its
    share of the union of both sets (i.e. of all three subset counts)."""
    total = sum(counts)
    if total == 0:
        return
    for subset_id, n in zip(("10", "01", "11"), counts):
        lbl = v.get_label_by_id(subset_id)
        if lbl is None:
            continue
        lbl.set_fontsize(fontsize)
        base = lbl.get_transform()
        lbl.set_transform(offset_copy(base, fig=ax.figure, y=5, units="points"))
        ax.text(*lbl.get_position(), f"({n / total:.1%})",
                transform=offset_copy(base, fig=ax.figure, y=-5, units="points"),
                ha="center", va="top", fontsize=pct_fontsize)


# ── Figure layout: outer 1×2, inner 2×2 for Venns ────────────────────────────
fig = plt.figure(figsize=(12, 6))
outer_gs = fig.add_gridspec(
    1, 2,
    width_ratios=[2, 1.2],
    left=0.04, right=0.97, top=0.87, bottom=0.15,
    wspace=0.25,
)
inner_gs = outer_gs[0, 0].subgridspec(2, 2, hspace=0.22, wspace=0.05)

venn_axes = [
    fig.add_subplot(inner_gs[0, 0]),
    fig.add_subplot(inner_gs[0, 1]),
    fig.add_subplot(inner_gs[1, 0]),
    fig.add_subplot(inner_gs[1, 1]),
]
ax_scatter = fig.add_subplot(outer_gs[0, 1])

# ── Venn diagrams ─────────────────────────────────────────────────────────────
venn_legend = [
    mpatches.Patch(facecolor=COLOR_AI, alpha=ALPHA, label="AI fraction"),
    mpatches.Patch(facecolor=COLOR_AS, alpha=ALPHA, label="AS fraction"),
]

both_sets    = []
group_tables = {}   # {title: {ai_only, as_only, both}}

for ax_v, (genotype, scaffold, title) in zip(venn_axes, GROUPS):
    ai_cols = samples_for("AI", genotype, scaffold)
    as_cols = samples_for("AS", genotype, scaffold)
    det_ai  = detected_proteins(ai_df, ai_cols)
    det_as  = detected_proteins(as_df, as_cols)

    ai_only = det_ai - det_as
    as_only = det_as - det_ai
    both    = det_ai & det_as

    v = venn2([det_ai, det_as], set_labels=None,
              set_colors=(COLOR_AI, COLOR_AS), alpha=ALPHA, ax=ax_v)
    label_venn_with_percentages(v, (len(ai_only), len(as_only), len(both)), ax_v)
    ax_v.set_title(title, fontsize=14, color="black", pad=2)

    both_sets.append(both)
    group_tables[title] = {"ai_only": ai_only, "as_only": as_only, "both": both}
    print(f"{title}: AI-only={len(ai_only)}, AS-only={len(as_only)}, both={len(both)}")

both_all_groups = both_sets[0].intersection(*both_sets[1:])
print(f"Detected in both fractions across all four groups: {len(both_all_groups)}")

as_only_all_groups = set.intersection(*[s["as_only"] for s in group_tables.values()])
print(f"AS-exclusive across all four groups: {len(as_only_all_groups)}")
pd.DataFrame([
    {"Protein.IDs": p, "Gene.Names": id_to_gene.get(p, "")}
    for p in sorted(as_only_all_groups)
]).to_csv(os.path.join(OUT_DIR, "AS_exclusive_all_groups.csv"), index=False)
print(f"Saved AS-exclusive list → {os.path.join(OUT_DIR, 'AS_exclusive_all_groups.csv')}")

# ── Scatter ───────────────────────────────────────────────────────────────────
all_x, all_y = [], []
scatter_legend = []

for (genotype, scaffold, title), color in zip(GROUPS, GROUP_COLORS):
    ai_cols = samples_for("AI", genotype, scaffold)
    as_cols = samples_for("AS", genotype, scaffold)
    det_ai  = detected_proteins(ai_df, ai_cols)
    det_as  = detected_proteins(as_df, as_cols)
    both    = det_ai & det_as

    ai_med = (ai_df[ai_df["Protein.IDs"].isin(both)]
              .set_index("Protein.IDs")[[c for c in ai_df.columns[N_META:] if c in ai_cols]]
              .median(axis=1, skipna=True))
    as_med = (as_df[as_df["Protein.IDs"].isin(both)]
              .set_index("Protein.IDs")[[c for c in as_df.columns[N_META:] if c in as_cols]]
              .median(axis=1, skipna=True))

    x = ai_med.reindex(as_med.index).values
    y = as_med.values
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]

    all_x.extend(x)
    all_y.extend(y)

    r = pearsonr(x, y)[0] if len(x) > 1 else float("nan")
    ax_scatter.scatter(x, y, color=color, s=5, alpha=0.65,
                       linewidths=0, rasterized=True)
    scatter_legend.append(
        mpatches.Patch(facecolor=color, alpha=0.8,
                       label=f"{title}  (r = {r:.2f})")
    )
    print(f"{title}: n={len(x)}, r={r:.3f}")

all_x = np.array(all_x)
all_y = np.array(all_y)
lo = min(all_x.min(), all_y.min()) - 0.3
hi = max(all_x.max(), all_y.max()) + 0.3

ax_scatter.plot([lo, hi], [lo, hi], color="#AAAAAA", lw=0.8, ls="--", zorder=0)
ax_scatter.set_xlim(lo, hi)
ax_scatter.set_ylim(lo, hi)
ax_scatter.set_xlabel("Median normalized log₂ abundance\nin AI fraction", fontsize=13)
ax_scatter.set_ylabel("Median normalized log₂ abundance\nin AS fraction", fontsize=13)
ax_scatter.tick_params(labelsize=12)
ax_scatter.spines[["top", "right"]].set_visible(False)
ax_scatter.set_aspect("equal")

# ── Legends and panel labels (after canvas.draw() for real positions) ─────────
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
inv = fig.transFigure.inverted()

# AI/AS legend: centred below the Venn block
vx0 = inv.transform((venn_axes[0].get_window_extent(renderer).x0, 0))[0]
vx1 = inv.transform((venn_axes[1].get_window_extent(renderer).x1, 0))[0]
venn_cx = (vx0 + vx1) / 2
venn_by = min(
    inv.transform((0, ax.get_window_extent(renderer).y0))[1]
    for ax in venn_axes[2:]
)
fig.legend(handles=venn_legend, loc="upper center",
           bbox_to_anchor=(venn_cx, venn_by - 0.01),
           ncol=2, fontsize=12, frameon=True)

ax_scatter.legend(handles=scatter_legend, loc="upper center",
                  bbox_to_anchor=(0.5, -0.22), ncol=1, fontsize=11, frameon=True)

# Panel labels
for ax, label in zip([venn_axes[0], ax_scatter], ["(a)", "(b)"]):
    x = inv.transform((ax.get_window_extent(renderer).x0, 0))[0]
    fig.text(x, 0.92, label, fontsize=16, fontweight="bold",
             va="bottom", ha="left")

out_path = os.path.join(OUT_DIR, "venn_AI_vs_AS_by_group.png")
plt.savefig(out_path, dpi=150, bbox_inches="tight")
plt.close()
print(f"Saved {out_path}")