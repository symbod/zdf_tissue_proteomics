"""
Dotplot of top enriched pathways per protein set and comparison.

Rows = union of top-10 pathways (significant in ALL 3 protein sets) per comparison,
sorted by mean −log10(p.adj) across the six cells. Columns = 3 protein sets × 2 comparisons.
Dot size = −log10(p.adj); grey dot = not significant.
GO depths downloaded once from go-basic.obo and cached at ~/.cache/go-basic.obo.
"""

import os
import re
import textwrap
import urllib.request
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict, deque


def normalize_pathway_id(pid):
    """Convert organism-specific KEGG IDs (e.g. rno00010) to KEGG: format (KEGG:00010).
    GO IDs pass through unchanged."""
    if re.match(r'^[a-z]{3}\d+$', str(pid)):
        return "KEGG:" + pid[3:]
    return pid

FC_THRESHOLDS   = [4.04, 5.59, 9.32]
STAB_THRESHOLDS = [1.1, 1.2, 1.3]

SETS      = ["set3_first_level", "set1_tissue_level", "set2_tissue_plus_network"]
SET_LABELS = [
    "First-level\nDEPs",
    "Tissue-level\nDEPs",
    "Tissue-level DEPs +\nconnector proteins",
]
COMP_KEYS = [
    "diabetic_empty_42-nondiabetic_empty_42",
    "diabetic_PCL_42-nondiabetic_PCL_42",
]
COMP_LBLS = ["Empty defect", "PCL scaffold"]
TOP_N     = 10

OBO_CACHE = os.path.expanduser("~/.cache/go-basic.obo")
OBO_URL   = "https://current.geneontology.org/ontology/go-basic.obo"
BONE_FILE = "input/bone_enrichments_meta_analysis.csv"

CAT_COLORS  = {"BP": "#0072B2", "MF": "#E69F00", "CC": "#56B4E9", "KEGG": "#D55E00"}
CAT_ORDER   = ["BP", "MF", "CC", "KEGG"]
CAT_LABELS  = {"BP": "GO:BP", "MF": "GO:MF", "CC": "GO:CC", "KEGG": "KEGG"}  # legend display names
BONE_STAR  = " ★"  # appended to pathway labels that are in the bone-healing reference

# 6 columns: 3 sets × 2 comparisons, with a visual gap between the two comparison groups
COL_KEYS      = [(ck, s) for ck in COMP_KEYS for s in SETS]
COL_POS       = [0.0, 0.42, 0.84, 1.42, 1.84, 2.26]
COL_TICK_LBLS = SET_LABELS * 2
COMP_CTR      = [0.42, 1.84]   # x-centre of each comparison group, for the header text
COMP_SEP      = 1.13           # x-position of the dashed divider between comparison groups


# ── GO depth ──────────────────────────────────────────────────────────────────

def ensure_obo() -> str:
    """
    Downloads the GO ontology file (go-basic.obo) if it is not already cached.
    Returns: path to the local .obo file (str).
    """
    if not os.path.exists(OBO_CACHE):
        os.makedirs(os.path.dirname(OBO_CACHE), exist_ok=True)
        print(f"Downloading go-basic.obo → {OBO_CACHE}  (~27 MB) …")
        urllib.request.urlretrieve(OBO_URL, OBO_CACHE)
    return OBO_CACHE


def parse_go_depths(obo_path: str) -> dict:
    """
    Receives: path to a go-basic.obo file (str).
    Returns:  dict mapping each GO term ID (e.g. "GO:0006096") to its depth (int),
              where depth = minimum number of is_a steps from one of the three
              root terms (biological_process / molecular_function / cellular_component).

    How it works:
      1. Reads all "is_a" parent–child edges from the OBO file.
      2. Runs BFS outward from the three root nodes, assigning depth level by level.
         BFS guarantees the shortest path, so each term gets its minimum depth.
    """
    # Pass 1: build parent → {children} map from all "is_a" lines
    children: dict = defaultdict(set)
    cur = None
    with open(obo_path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line == "[Term]":
                cur = None
            elif line.startswith("id: GO:"):
                cur = line[4:]
            elif line.startswith("is_a:") and cur:
                children[line.split()[1]].add(cur)

    # Pass 2: BFS from the three roots, depth 0 → 1 → 2 → …
    roots = {"GO:0008150", "GO:0003674", "GO:0005575"}
    depth: dict = {r: 0 for r in roots}
    q: deque = deque((r, 0) for r in roots)
    while q:
        node, d = q.popleft()
        for child in children.get(node, ()):
            if child not in depth:
                depth[child] = d + 1
                q.append((child, d + 1))
    return depth


# ── data loading ──────────────────────────────────────────────────────────────

def load_all(go_depths: dict, base: str) -> pd.DataFrame:
    """
    Receives: go_depths — dict {GO_ID: depth} from parse_go_depths().
              base      — path to the enrichment directory for one threshold combination.
    Returns:  a single DataFrame with one row per pathway per (comparison, set),
              i.e. 6 × (number of tested pathways) rows in total.
              Columns from the CSV are kept as-is; three columns are added:
                - "Set"         : which protein set this row belongs to
                - "Comparison"  : which comparison (e.g. diabetic_empty_42-...)
                - "go_depth"    : GO depth from go_depths (NaN for KEGG terms)
                - "neg_log_padj": −log10(p.adjust), used as dot size in the figure
    """
    frames = []
    for ck in COMP_KEYS:
        for s in SETS:
            fp = os.path.join(base, s, f"{ck}_{s}_enrichment_complete.csv")
            df = pd.read_csv(fp, index_col=0)
            df["Set"] = s
            df["Comparison"] = ck
            frames.append(df)
    data = pd.concat(frames, ignore_index=True)

    # KEGG pathways have no GO:ID, so their depth is undefined → NaN
    data["go_depth"] = data.apply(
        lambda r: float(go_depths.get(r["ID"], np.nan))
        if r["Category"] != "KEGG" else np.nan,
        axis=1,
    )
    # clip prevents log(0) for extremely small p-values from clusterProfiler
    data["neg_log_padj"] = -np.log10(data["p.adjust"].clip(lower=1e-300))
    return data


# ── top-pathway selection ─────────────────────────────────────────────────────

def top_ids_for_comp(data: pd.DataFrame, comp: str) -> list:
    """
    Receives: data — the full stacked DataFrame from load_all().
              comp — one comparison key string (e.g. "diabetic_empty_42-nondiabetic_empty_42").
    Returns:  list of up to TOP_N pathway ID strings, ranked by significance.

    Selection logic:
      1. Filter to rows that are significant (p.adj < threshold) for this comparison.
      2. Pivot to a matrix: rows = pathway IDs, columns = the 3 protein sets,
         values = neg_log_padj.
      3. Drop any pathway that is missing from at least one set (dropna) — only
         pathways significant in ALL 3 protein sets simultaneously are kept.
      4. Score each remaining pathway by its mean neg_log_padj across the 3 sets.
      5. Return the IDs of the top TOP_N pathways by score.
    """
    sig = data[(data["Comparison"] == comp) & (data["Significant"] == True)]
    piv = (
        sig.pivot_table(index="ID", columns="Set", values="neg_log_padj", aggfunc="mean")
        .reindex(columns=SETS)
        .dropna()                          # require significance in ALL 3 sets
    )
    piv["score"] = piv[SETS].mean(axis=1)  # rank by mean −log10(p.adj) across sets
    return list(piv.nlargest(TOP_N, "score").index)


def build_dot_df(data: pd.DataFrame, bone_ids: set) -> pd.DataFrame:
    """
    Receives: data     — the full stacked DataFrame from load_all().
              bone_ids — set of pathway ID strings from the bone-healing reference.
    Returns:  a DataFrame with one row per selected pathway, ready for the dotplot.
              Columns:
                - ID, Description, Category, go_depth — pathway metadata
                - one column per (comparison × set) combination, named
                  "<comp_key>__<set_key>", containing neg_log_padj if the pathway
                  is significant there, NaN otherwise (→ grey dot in the figure)
                - in_bone_meta: True if this pathway is in the bone-healing reference

    How pathways are selected:
      Call top_ids_for_comp() for each comparison to get up to TOP_N pathway IDs
      each, then merge into one deduplicated list (up to 20 total).
    """
    # Step 1: collect top-N IDs per comparison into one deduplicated ordered list
    all_ids: list = []
    for ck in COMP_KEYS:
        for pid in top_ids_for_comp(data, ck):
            if pid not in all_ids:
                all_ids.append(pid)

    # Step 2: look up Description, Category, go_depth from set1.
    # These fields are pathway-intrinsic (same across all sets), so set1 is sufficient.
    meta = (
        data[data["Set"] == SETS[0]]
        [["ID", "Description", "Category", "go_depth"]]
        .drop_duplicates("ID")
        .set_index("ID")
    )

    # Step 3: for each pathway, store neg_log_padj for every (comparison × set) cell.
    # NaN = pathway not significant in that cell → will be drawn as a grey dot.
    rows = []
    for pid in all_ids:
        if pid not in meta.index:
            continue
        row = {"ID": pid, **meta.loc[pid].to_dict()}
        for ck in COMP_KEYS:
            for s in SETS:
                cell = data[(data["ID"] == pid) & (data["Comparison"] == ck) & (data["Set"] == s)]
                row[f"{ck}__{s}"] = (
                    float(cell.iloc[0]["neg_log_padj"])
                    if len(cell) and bool(cell.iloc[0]["Significant"])
                    else np.nan
                )
        rows.append(row)

    df = pd.DataFrame(rows)
    # Sort by mean neg_log_padj across all 6 (comparison × set) columns, descending,
    # so the most enriched pathway appears at the top of the figure.
    # NaN cells (not significant) are excluded from the mean.
    val_cols = [f"{ck}__{s}" for ck, s in COL_KEYS]
    df["mean_score"] = df[val_cols].mean(axis=1)
    df = df.sort_values("mean_score", ascending=False).reset_index(drop=True)
    df = df.drop(columns="mean_score")
    df["in_bone_meta"] = df["ID"].map(normalize_pathway_id).isin(bone_ids)
    return df


# ── figure helpers ────────────────────────────────────────────────────────────

def _set_col_ticks(ax: plt.Axes, rotation: int = 45) -> None:
    """Sets the 6 x-axis tick labels (protein set names, repeated for each comparison)."""
    ax.set_xticks(COL_POS)
    ax.set_xticklabels(COL_TICK_LBLS, fontsize=6.5, rotation=rotation, ha="right",
                       rotation_mode="anchor")


def _add_comp_headers(ax: plt.Axes, y: float) -> None:
    """Draws the comparison group labels ("Empty defect", "PCL scaffold") above the columns."""
    for cx, cl in zip(COMP_CTR, COMP_LBLS):
        ax.text(cx, y, cl, ha="center", va="bottom", fontsize=7,
                fontweight="bold", color="0.25")


def _wrap_label(text: str, width: int = 26) -> str:
    """Wraps a pathway description to at most 2 lines for use as a y-axis tick label."""
    lines = textwrap.wrap(text, width=width, max_lines=2, placeholder="…")
    return "\n".join(lines)


# ── figure ────────────────────────────────────────────────────────────────────

def make_dotplot_figure(dot_df: pd.DataFrame) -> plt.Figure:
    """
    Receives: dot_df — the DataFrame produced by build_dot_df().
    Returns:  a matplotlib Figure object (not yet saved).

    Layout:
      - Y-axis: one row per pathway (pathway name as tick label)
      - X-axis: 6 columns — 3 protein sets × 2 comparisons
      - Left margin: coloured depth stripe (colour = namespace, alpha = GO depth)
                     + brown indicator for bone-healing reference pathways
      - Dots: coloured and sized by neg_log_padj if significant; small grey if not
      - Bottom: two legends — one for namespace colours, one for dot size scale
    """
    val_cols = [f"{ck}__{s}" for ck, s in COL_KEYS]
    max_val  = dot_df[val_cols].max().max()  # largest neg_log_padj → largest dot size
    n_rows   = len(dot_df)

    # figure height scales with row count so pathway labels never overlap
    fig, ax = plt.subplots(figsize=(4.6, max(3.2, 0.28 * n_rows + 1.3)))

    # ── left margin: depth stripe (colour = namespace, number = GO depth) ───────
    STRIPE_X = min(COL_POS) - 0.36
    STRIPE_W = 0.20

    for ri, row in dot_df.iterrows():
        d = row["go_depth"]
        ax.add_patch(plt.Rectangle(
            (STRIPE_X - STRIPE_W / 2, ri - 0.42), STRIPE_W, 0.84,
            facecolor=CAT_COLORS[row["Category"]], alpha=0.82, lw=0, zorder=2,
        ))
        ax.text(STRIPE_X, ri, str(int(d)) if not np.isnan(d) else "–",
                ha="center", va="center", fontsize=5.5, color="white",
                fontweight="bold", zorder=3)

    ax.text(STRIPE_X, -0.58, "GO\ndepth", ha="center", va="bottom",
            fontsize=5.5, linespacing=1.2)

    # ── dots: coloured + sized if significant, small grey if not ─────────────
    for ri, row in dot_df.iterrows():
        color = CAT_COLORS[row["Category"]]
        for xi, col in zip(COL_POS, val_cols):
            v = row[col]
            if np.isnan(v):
                ax.scatter(xi, ri, s=16, color="0.87", zorder=2, linewidths=0)
            else:
                ax.scatter(xi, ri, s=(v / max_val) * 150 + 11,
                           color=color, alpha=0.88, zorder=3, linewidths=0)

    ax.set_yticks(range(n_rows))
    ax.set_yticklabels(
        [
            _wrap_label(row["Description"]) + (BONE_STAR if row["in_bone_meta"] else "")
            for _, row in dot_df.iterrows()
        ],
        fontsize=6.5,
    )
    ax.set_xlim(min(COL_POS) - 0.58, max(COL_POS) + 0.18)
    ax.set_ylim(n_rows - 0.45, -0.55)

    _set_col_ticks(ax)
    ax.set_xlabel("Protein set", fontsize=8, labelpad=6)
    ax.axvline(COMP_SEP, color="0.45", lw=0.75, linestyle="--")

    _add_comp_headers(ax, y=-0.45)
    ax.grid(axis="x", lw=0.3, color="0.92", zorder=0)
    ax.spines[["top", "right"]].set_visible(False)

    # ── legends on the right ──────────────────────────────────────────────────
    # tight_layout constrains axes to the left 80 % of the figure width;
    # the remaining strip on the right holds the two legends and the star note.
    fig.tight_layout(rect=[0, 0, 0.80, 1])

    ns_handles = [mpatches.Patch(facecolor=CAT_COLORS[c], label=CAT_LABELS[c]) for c in CAT_ORDER]
    leg1 = fig.legend(handles=ns_handles, fontsize=6, frameon=False, ncol=1,
                      title="Term type", title_fontsize=6.5,
                      bbox_to_anchor=(0.81, 0.97), loc="upper left")
    leg1.get_title().set_fontweight("bold")

    # ── size legend: manual inset so blob-edge gaps are exactly equal ────────
    # s values (pt²) and radii (pt) mirror the main-plot formula exactly
    _fracs   = [1.0, 0.60, 0.25]
    _lbls    = [f"{f * max_val:.0f}" for f in _fracs]
    _dot_s   = [f * 150 + 11 for f in _fracs]
    _radii   = [np.sqrt(s / np.pi) for s in _dot_s]

    _EDGE_GAP_PT = 6   # equal gap between successive blob edges, in display points
    _TITLE_H_PT  = 9   # vertical space for the title

    # Cumulative centre positions (pt from top of inset, downward)
    _centers = [_TITLE_H_PT + _radii[0]]
    for _i in range(1, len(_radii)):
        _centers.append(_centers[-1] + _radii[_i - 1] + _EDGE_GAP_PT + _radii[_i])
    _total_h_pt = _centers[-1] + _radii[-1] + 4  # 4 pt bottom padding

    # Size the inset so that 1 axes unit == 1 display point → scatter s is correct
    _fh_in = fig.get_figheight()
    _inset_h = _total_h_pt / (_fh_in * 72)   # figure fraction
    _size_ax = fig.add_axes([0.81, 0.58 - _inset_h, 0.17, _inset_h])
    _size_ax.set_xlim(0, 1)
    _size_ax.set_ylim(_total_h_pt, 0)  # 0 at top
    _size_ax.axis("off")

    _size_ax.text(0.05, 0, "−log₁₀(q-value)",
                  fontsize=6.5, fontweight="bold", va="top", ha="left")
    for _c, _s, _lbl in zip(_centers, _dot_s, _lbls):
        _size_ax.scatter([0.22], [_c], s=_s, color="0.55",
                         zorder=3, linewidths=0, clip_on=False)
        _size_ax.text(0.44, _c, _lbl, fontsize=6, va="center", ha="left")

    fig.text(0.81, 0.24, "★ Bone-healing\nreference\npathway",
             fontsize=6.5, fontweight="bold", va="top", ha="left", color="0.35",
             transform=fig.transFigure, linespacing=1.4)

    return fig


# ── main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    go_depths = parse_go_depths(ensure_obo())
    bone_ids = set(pd.read_csv(BONE_FILE)["term_id"].dropna())

    for fc in FC_THRESHOLDS:
        for stab in STAB_THRESHOLDS:
            base = os.path.join("valid_DEPs", f"FC{fc}_Stab{stab}", "enrichment")

            # Check all required enrichment files exist before attempting to load
            missing = [
                os.path.join(base, s, f"{ck}_{s}_enrichment_complete.csv")
                for ck in COMP_KEYS for s in SETS
                if not os.path.exists(os.path.join(base, s, f"{ck}_{s}_enrichment_complete.csv"))
            ]
            if missing:
                print(f"\nSkipping FC{fc}_Stab{stab} — {len(missing)} enrichment file(s) missing")
                continue

            print(f"\n{'='*60}")
            print(f"Processing FC{fc}_Stab{stab} …")
            print(f"{'='*60}")

            data = load_all(go_depths, base)

            print("Selecting top shared pathways …")
            dot_df = build_dot_df(data, bone_ids)

            if dot_df.empty:
                print("  No pathways shared across all 3 sets in both comparisons — skipping figure")
                continue

            print(f"  {len(dot_df)} unique pathways in union (top {TOP_N} per comparison)")
            print(f"  {dot_df['in_bone_meta'].sum()} of these are in the bone-healing reference")

            out = os.path.join(base, "pathway_hierarchy_dotplot.png")
            make_dotplot_figure(dot_df).savefig(out, dpi=300, bbox_inches="tight")
            plt.close("all")
            print(f"Saved → {out}")