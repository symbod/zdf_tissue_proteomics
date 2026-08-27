"""
DEP fraction breakdown — single combined figure.
One row per comparison: empty defect on top (panels a, b), PCL scaffold below
(panels c, d).
Left column (a, c):  Venn diagram of the DEP overlap between AI and AS fractions.
Right column (b, d): stacked bars broken down by category
                     Both | AI-only | AS-only, with concordant/discordant
                     and detected/exclusive sub-segments.

Reads validation_table CSVs exported by tissue_level_DEP_collection.R.
dep_type / ai_exclusive / as_exclusive are threshold-independent.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.font_manager import FontProperties
from matplotlib_venn import venn2
from matplotlib.transforms import offset_copy

OUT_DIR        = "ai_as_proteome_investigation_results"
VALIDATION_DIR = "valid_DEPs/FC5.59_Stab1.2"

os.makedirs(OUT_DIR, exist_ok=True)

COMPARISONS = [
    "diabetic_empty_42-nondiabetic_empty_42",
    "diabetic_PCL_42-nondiabetic_PCL_42",
]
PANEL_LABELS = {
    "diabetic_empty_42-nondiabetic_empty_42": "Empty defect",
    "diabetic_PCL_42-nondiabetic_PCL_42":     "PCL scaffold",
}

COLORS = {
    "concordant":    "#009E73",  # bluish green   — both fractions, same direction
    "discordant":    "#D55E00",  # vermillion     — both fractions, opposite direction
    "ai_only_det":   "#56B4E9",  # sky blue       — AI DEP, detected in AS
    "ai_only_excl":  "#0072B2",  # blue           — AI DEP, never detected in AS
    "as_only_det":   "#E69F00",  # orange         — AS DEP, detected in AI
    "as_only_excl":  "#CC79A7",  # reddish purple — AS DEP, never detected in AI
}
LEGEND_LABELS = {
    "concordant":   "same direction",
    "discordant":   "opposite direction",
    "ai_only_det":  "detected in AS",
    "ai_only_excl": "not detected in AS",
    "as_only_det":  "detected in AI",
    "as_only_excl": "not detected in AI",
}

COLOR_AI = "#0072B2"
COLOR_AS = "#D55E00"

BAR_W     = 0.70
MIN_LABEL = 15
CAT_X     = [0, 1.0, 2.0]
CAT_GROUPS = [
    ("Both\nfractions",   ["concordant",  "discordant"]),
    ("AI fraction\nonly", ["ai_only_det", "ai_only_excl"]),
    ("AS fraction\nonly", ["as_only_det", "as_only_excl"]),
]


def load_validation_table(comp):
    df = pd.read_csv(os.path.join(VALIDATION_DIR, f"validation_table_{comp}.csv"))
    for col in ("ai_exclusive", "as_exclusive"):
        if df[col].dtype == object:
            df[col] = df[col].str.upper() == "TRUE"
    return df


def compute_breakdown(comp):
    vt = load_validation_table(comp)
    concordant   = (vt["Validation"] == "Significant in both fractions (same direction)").sum()
    ai_only_excl = (vt["Validation"] == "Exclusive to AI fraction").sum()
    as_only_excl = (vt["Validation"] == "Exclusive to AS fraction").sum()
    discordant   = (vt["dep_type"] == "Both_opposite").sum()
    ai_only_det  = ((vt["dep_type"] == "AI_only") & (vt["Validation"] != "Exclusive to AI fraction")).sum()
    as_only_det  = ((vt["dep_type"] == "AS_only") & (vt["Validation"] != "Exclusive to AS fraction")).sum()
    d = {
        "concordant":   int(concordant),
        "discordant":   int(discordant),
        "ai_only_det":  int(ai_only_det),
        "ai_only_excl": int(ai_only_excl),
        "as_only_det":  int(as_only_det),
        "as_only_excl": int(as_only_excl),
    }
    d["shared"]  = d["concordant"] + d["discordant"]
    d["ai_only"] = d["ai_only_det"] + d["ai_only_excl"]
    d["as_only"] = d["as_only_det"] + d["as_only_excl"]
    return d


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


all_data = {comp: compute_breakdown(comp) for comp in COMPARISONS}

ylim_cat = max(
    max(d["shared"], d["ai_only"], d["as_only"])
    for d in all_data.values()
) * 1.18


# ── Figure & gridspec ─────────────────────────────────────────────────────────
# Col 0 = Venns stacked (panel a), Col 1 = bars stacked (panel b)
fig = plt.figure(figsize=(12, 7))
gs = fig.add_gridspec(
    2, 2,
    width_ratios=[1, 1.55],
    left=0.05, right=0.70,
    top=0.87, bottom=0.08,
    hspace=0.30, wspace=0.28,
)
venn_axes = [fig.add_subplot(gs[i, 0]) for i in range(2)]
bar_axes  = [fig.add_subplot(gs[0, 1])]
bar_axes.append(fig.add_subplot(gs[1, 1], sharey=bar_axes[0]))


# ── Panel (a): Venn diagrams ──────────────────────────────────────────────────
for ax, comp in zip(venn_axes, COMPARISONS):
    d = all_data[comp]
    v = venn2(
        subsets=(d["ai_only"], d["as_only"], d["shared"]),
        set_labels=("AI fraction", "AS fraction"),
        set_colors=(COLOR_AI, COLOR_AS),
        alpha=0.55, ax=ax,
    )
    label_venn_with_percentages(v, (d["ai_only"], d["as_only"], d["shared"]), ax)
    for lbl in v.set_labels:
        if lbl:
            lbl.set_fontsize(11)


# ── Panel (b): Category stacked bars ─────────────────────────────────────────
for ax, comp in zip(bar_axes, COMPARISONS):
    d = all_data[comp]

    for x_pos, (_, keys) in zip(CAT_X, CAT_GROUPS):
        bottom = 0
        for key in keys:
            count = d[key]
            if count == 0:
                continue
            ax.bar(x_pos, count, bottom=bottom, width=BAR_W,
                   color=COLORS[key], edgecolor="white", linewidth=1.0, zorder=2)
            if count >= MIN_LABEL:
                ax.text(x_pos, bottom + count / 2, str(count),
                        ha="center", va="center", fontsize=10,
                        color="white", fontweight="bold", zorder=3)
            else:
                ax.annotate(
                    str(count),
                    xy=(x_pos + BAR_W / 2, bottom + count / 2),
                    xytext=(x_pos + BAR_W / 2 + 0.10, bottom + count / 2),
                    ha="left", va="center", fontsize=9, fontweight="bold",
                    color=COLORS[key], zorder=3,
                    arrowprops=dict(arrowstyle="-", color=COLORS[key], lw=0.8),
                )
            bottom += count
        ax.text(x_pos, bottom + ylim_cat * 0.015, f"n={bottom}",
                ha="center", va="bottom", fontsize=10.5, fontweight="bold")

    ax.set_xticks(CAT_X)
    ax.set_xticklabels([g[0] for g in CAT_GROUPS], fontsize=11)
    ax.set_ylim(0, ylim_cat)
    ax.set_xlim(-0.5, 2.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.yaxis.grid(True, linestyle=":", linewidth=0.6, alpha=0.5, zorder=0)
    ax.set_axisbelow(True)
    ax.set_ylabel("Number of DEPs", fontsize=12)


# ── Panel labels (a)(b)(c)(d) + row titles ────────────────────────────────────
# Anchored to the gridspec cell geometry (not the rendered Axes bbox): venn2
# applies aspect="equal", which shrinks the Venn axes box around its center,
# so a title/label hung off that axes' own top edge drifts out of line with
# the bar-chart title in the same row. The gridspec cell top is unaffected by
# that shrink, so both columns line up exactly.
fig.canvas.draw()
renderer = fig.canvas.get_renderer()
inv = fig.transFigure.inverted()

row_labels = [("a", "b"), ("c", "d")]

for i, comp in enumerate(COMPARISONS):
    pos_venn = gs[i, 0].get_position(fig)
    pos_bar  = gs[i, 1].get_position(fig)
    row_top  = pos_venn.y1  # == pos_bar.y1, same gridspec row

    for pos, label in zip((pos_venn, pos_bar), row_labels[i]):
        x_center = (pos.x0 + pos.x1) / 2
        fig.text(x_center, row_top + 0.006, PANEL_LABELS[comp],
                  fontsize=13, fontweight="bold", ha="center", va="bottom")
        fig.text(pos.x0 - 0.012, row_top + 0.032, f"({label})",
                  fontsize=16, fontweight="bold", ha="right", va="bottom")


# ── Legend for panel (b) ──────────────────────────────────────────────────────
title_fp = FontProperties(weight="bold", size=10.5)
legend_groups = [
    ("DEPs in both fractions:",   "concordant",  "discordant"),
    ("DEPs in AI fraction only:", "ai_only_det", "ai_only_excl"),
    ("DEPs in AS fraction only:", "as_only_det", "as_only_excl"),
]

# Anchor legend groups within the vertical span of the bar column
top_bbox    = bar_axes[0].get_window_extent(renderer)
bottom_bbox = bar_axes[1].get_window_extent(renderer)
col_top    = inv.transform((0, top_bbox.y1))[1]
col_bottom = inv.transform((0, bottom_bbox.y0))[1]
col_span   = col_top - col_bottom
leg_y = [col_top - 0.01,
         col_top - col_span * 0.37,
         col_top - col_span * 0.72]

for (title_text, k1, k2), y_anc in zip(legend_groups, leg_y):
    fig.legend(
        handles=[
            mpatches.Patch(facecolor=COLORS[k1], edgecolor="white", label=LEGEND_LABELS[k1]),
            mpatches.Patch(facecolor=COLORS[k2], edgecolor="white", label=LEGEND_LABELS[k2]),
        ],
        title=title_text,
        title_fontproperties=title_fp,
        loc="upper left",
        bbox_to_anchor=(0.715, y_anc),
        ncol=1, fontsize=10.5, frameon=False,
        handlelength=1.2, handletextpad=0.5, borderpad=0.3, labelspacing=0.3,
    )


out = os.path.join(OUT_DIR, "dep_fraction_breakdown.png")
plt.savefig(out, dpi=180, bbox_inches="tight")
plt.close()
print(f"Saved {out}")