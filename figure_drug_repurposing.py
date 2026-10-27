"""
Scatter plot of the drug repurposing candidates: score in the empty defect vs. score with the PCL scaffold.

One marker per drug of Supplementary Table S7 (the candidates of either comparison): x = mean TrustRank score in
the empty defect, y = with the PCL scaffold, each the mean of the drug's scores in the three protein-set runs of
the comparison (1 = top-ranked drug of a run). In the comparison in which a drug is a candidate this is the
"Mean score" of S7, read from the intersection table of drug_repurposing.py; in the other comparison it is
computed from the per-run tables with a run that does not return the drug counted as 0. Drugs with identical
scores (drugs that bind the same target proteins) stack on one marker; the marker faces are translucent, so a
stack is darker than a single drug. Labelled by name are the top TOP_N drugs of each comparison, drugs tied with
the TOP_N-th included (like the union of the top-10 terms in figure_pathway_hierarchy.py); a stack of more than
MAX_NAMES_INLINE drugs is labelled by its count, and labelled markers that touch share one label.
Reads the drug_repurposing/ outputs of drug_repurposing.py for every threshold combination.
"""

import itertools
import os
import textwrap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba

FC_THRESHOLDS   = [4.04, 5.59, 9.32]
STAB_THRESHOLDS = [1.1, 1.2, 1.3]

COMP_KEYS = [
    "diabetic_empty_42-nondiabetic_empty_42",
    "diabetic_PCL_42-nondiabetic_PCL_42",
]
COMP_LBLS = ["Empty defect", "PCL scaffold"]
SETS = ["set1_tissue_level", "set2_tissue_plus_network", "set3_first_level"]

TOP_N = 3              # drugs labelled per comparison (ranked by mean score; drugs tied with the TOP_N-th are included)
MAX_NAMES_INLINE = 4   # a stack of more drugs with identical scores is labelled by its count
LABEL_WRAP = 54        # characters per label line
LABEL_GAP = 5.5        # distance between a marker centre and its label (pt)
LABEL_PAD = 1.2        # clearance between a label and markers or other labels (pt)
LABEL_DIRECTIONS = [(1, 0), (-1, 0), (0, 1), (0, -1),            # positions tried for a label, in this order:
                    (1, 1), (-1, 1), (1, -1), (-1, -1)]          # right, left, above, below, then the diagonals
FIGSIZE = (6.3, 2.1)   # inches; the text width of a Google Docs / A4 page, the y axis is not to scale with the x axis
MARKER_AREA = 24       # scatter marker size (pt²)
MARKER_EDGE = 0.6      # marker edge width (pt)
FACE_ALPHA = 0.35      # opacity of a marker face; stacked markers appear darker
MARKER_COLOR = "#4F6D8F"
INK = "0.15"
GRID_COLOR = "0.90"
AXIS_LIMITS = (-0.02, 1.02)


# ── data ──────────────────────────────────────────────────────────────────────

def drug_table(drug_dir: str) -> pd.DataFrame:
    """
    Receives: drug_repurposing/ directory of one threshold combination.
    Returns:  one row per candidate drug of either comparison with, per comparison label, candidate_<label> (bool)
              and score_<label>: the "Mean score" of the intersection table written by drug_repurposing.py
              (Supplementary Table S7) where the drug is a candidate, else the same mean over the three runs with a
              run that does not return the drug counted as 0; and labelled: among the top TOP_N candidates of either
              comparison, drugs tied with the TOP_N-th included.
    """
    rows = {}
    for comparison, label in zip(COMP_KEYS, COMP_LBLS):
        candidates = pd.read_csv(os.path.join(drug_dir, f"{comparison}_top_drugs_intersection.csv"))
        for drug, mean_score in zip(candidates["Drug"], candidates["Mean score"]):
            row = rows.setdefault(drug, {"Drug": drug})
            row[f"candidate_{label}"] = True
            row[f"score_{label}"] = mean_score
    for comparison, label in zip(COMP_KEYS, COMP_LBLS):
        runs = [pd.read_csv(os.path.join(drug_dir, f"{comparison}_{protein_set}_drugs.csv")) for protein_set in SETS]
        score_in_run = [dict(zip(run["Drug"], run["Score"])) for run in runs]
        for drug, row in rows.items():
            row.setdefault(f"candidate_{label}", False)
            row.setdefault(f"score_{label}", sum(scores.get(drug, 0.0) for scores in score_in_run) / len(SETS))
    table = pd.DataFrame(rows.values())
    table["labelled"] = False
    for label in COMP_LBLS:
        candidate_scores = table[f"score_{label}"].where(table[f"candidate_{label}"])
        table["labelled"] |= candidate_scores.rank(method="min", ascending=False) <= TOP_N
    return table


def marker_table(table: pd.DataFrame) -> pd.DataFrame:
    """One row per marker, i.e. per distinct coordinate (scores rounded to 4 decimals): x, y, drugs (names,
    alphabetical) and labelled; ordered from the bottom up and, at equal height, from the right."""
    table = table.assign(x=table[f"score_{COMP_LBLS[0]}"].round(4), y=table[f"score_{COMP_LBLS[1]}"].round(4))
    markers = table.groupby(["x", "y"], as_index=False).agg(
        drugs=("Drug", lambda names: sorted(names, key=str.lower)), labelled=("labelled", "any"))
    return markers.sort_values(["y", "x"], ascending=[True, False]).reset_index(drop=True)


def label_text(marker: pd.Series) -> str:
    """The drug names of the marker ("Bortezomib, Quercetin") or the count of a large stack ("29 drugs"); a name is
    never broken across lines."""
    names = marker["drugs"]
    if len(names) > MAX_NAMES_INLINE:
        return f"{len(names)} drugs"
    return textwrap.fill(", ".join(name.replace(" ", " ") for name in names), width=LABEL_WRAP)


# ── labels ────────────────────────────────────────────────────────────────────
# The layout is computed in points (display coordinates × 72 / dpi), so it does not depend on the data ranges or
# the aspect ratio of the axes. Boxes are (x0, y0, x1, y1) tuples.

def label_box(centre: tuple, size: tuple, direction: tuple, frame: tuple) -> tuple:
    """Box of a label of the given size (width, height) next to the marker at centre: LABEL_GAP away in the given
    direction (dx, dy with values -1, 0, 1) and centred on the marker along an axis with 0, where it is shifted
    into the frame (axes box) if necessary; the marker then still lies within the label's extent."""
    width, height = size
    dx, dy = direction
    x0 = centre[0] + (LABEL_GAP if dx > 0 else -LABEL_GAP - width if dx < 0 else -width / 2)
    y0 = centre[1] + (LABEL_GAP if dy > 0 else -LABEL_GAP - height if dy < 0 else -height / 2)
    if dx == 0:
        x0 = min(max(x0, frame[0]), frame[2] - width)
    if dy == 0:
        y0 = min(max(y0, frame[1]), frame[3] - height)
    return (x0, y0, x0 + width, y0 + height)


def boxes_overlap(a: tuple, b: tuple, pad: float) -> bool:
    return a[0] - pad < b[2] and a[2] + pad > b[0] and a[1] - pad < b[3] and a[3] + pad > b[1]


def inside(box: tuple, frame: tuple) -> bool:
    return box[0] >= frame[0] and box[1] >= frame[1] and box[2] <= frame[2] and box[3] <= frame[3]


def label_groups(ax: plt.Axes, markers: pd.DataFrame, centres: list, radius: float) -> list:
    """
    Receives: the axes (drawn), the marker table, the marker centres (points) and the marker radius (points).
    Returns:  one dict per label: text (the Text artist, not yet positioned), size (width, height in points) and
              centre (points). Labelled markers that touch or overlap (centres closer than two radii, in chains)
              cannot be told apart and share one label, one line per marker in the order of their score, at their
              mean centre.
    """
    renderer = ax.figure.canvas.get_renderer()
    points_per_pixel = 72 / ax.figure.dpi
    groups = [[index] for index in markers.index[markers["labelled"]]]
    merged = True
    while merged:
        merged = False
        for a, b in itertools.combinations(range(len(groups)), 2):
            if any(np.hypot(*np.subtract(centres[i], centres[j])) < 2 * radius for i in groups[a] for j in groups[b]):
                groups[a].extend(groups.pop(b))
                merged = True
                break
    labels = []
    for members in groups:
        by_score = markers.loc[members].sort_values(["x", "y"], ascending=False)
        text = ax.text(0, 0, "\n".join(label_text(marker) for _, marker in by_score.iterrows()), fontsize=6.5,
                       color=INK, linespacing=1.2, ha="left", va="bottom", zorder=6)
        extent = text.get_window_extent(renderer)
        labels.append({"text": text, "size": (extent.width * points_per_pixel, extent.height * points_per_pixel),
                       "centre": tuple(np.mean([centres[index] for index in members], axis=0))})
    return labels


def place_labels(ax: plt.Axes, markers: pd.DataFrame) -> None:
    """
    Labels the labelled markers, lowest and right-most first. Each label takes the first position of
    LABEL_DIRECTIONS at which it lies inside the axes and overlaps neither a marker nor a label placed before;
    if there is none it is put to the right of its marker. Must be called after the axes geometry is final.
    """
    fig = ax.figure
    fig.canvas.draw()
    points_per_pixel = 72 / fig.dpi
    centres = [tuple(ax.transData.transform((x, y)) * points_per_pixel) for x, y in zip(markers["x"], markers["y"])]
    radius = np.sqrt(MARKER_AREA) / 2 + MARKER_EDGE / 2
    marker_boxes = [(x - radius, y - radius, x + radius, y + radius) for x, y in centres]
    axes_box = tuple(np.asarray(ax.get_window_extent().extents) * points_per_pixel)
    to_data = ax.transData.inverted()

    placed = []
    for label in label_groups(ax, markers, centres, radius):
        for direction in LABEL_DIRECTIONS:
            box = label_box(label["centre"], label["size"], direction, axes_box)
            free = not any(boxes_overlap(box, other, LABEL_PAD) for other in marker_boxes + placed)
            if inside(box, axes_box) and free:
                break
        else:
            box = label_box(label["centre"], label["size"], LABEL_DIRECTIONS[0], axes_box)
            print(f"  WARNING: no free position for the label {label['text'].get_text()!r}")
        placed.append(box)
        label["text"].set_position(to_data.transform(np.asarray(box[:2]) / points_per_pixel))


# ── figure ────────────────────────────────────────────────────────────────────

def make_figure(table: pd.DataFrame, markers: pd.DataFrame) -> plt.Figure:
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.scatter(table[f"score_{COMP_LBLS[0]}"], table[f"score_{COMP_LBLS[1]}"], s=MARKER_AREA,
               facecolors=to_rgba(MARKER_COLOR, FACE_ALPHA), edgecolors=MARKER_COLOR, linewidths=MARKER_EDGE, zorder=5,
               clip_on=False)
    ax.set_xlim(*AXIS_LIMITS)
    ax.set_ylim(*AXIS_LIMITS)
    ticks = [0, 0.25, 0.5, 0.75, 1.0]
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.tick_params(labelsize=6.5, colors=INK, length=2.5, width=0.5)
    # two lines each: the one-line title would be longer than the short y axis
    ax.set_xlabel(f"Mean TrustRank score\n{COMP_LBLS[0]}", fontsize=7.5, color=INK, labelpad=4, linespacing=1.1)
    ax.set_ylabel(f"Mean TrustRank score\n{COMP_LBLS[1]}", fontsize=7.5, color=INK, labelpad=4, linespacing=1.1)
    ax.grid(lw=0.4, color=GRID_COLOR, zorder=0)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_color("0.75")
        spine.set_linewidth(0.5)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    place_labels(ax, markers)
    return fig


# ── main ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    for fc in FC_THRESHOLDS:
        for stab in STAB_THRESHOLDS:
            drug_dir = os.path.join("valid_DEPs", f"FC{fc}_Stab{stab}", "drug_repurposing")
            missing = [
                comparison for comparison in COMP_KEYS
                if not os.path.exists(os.path.join(drug_dir, f"{comparison}_top_drugs_intersection.csv"))
            ]
            if missing:
                print(f"\nSkipping FC{fc}_Stab{stab} — no intersection table for {missing}")
                continue

            print(f"\n{'=' * 60}\nProcessing FC{fc}_Stab{stab} …\n{'=' * 60}")
            table = drug_table(drug_dir)
            markers = marker_table(table)
            print("Candidates: " + ", ".join(f"{int(table[f'candidate_{label}'].sum())} {label}" for label in COMP_LBLS)
                  + f"; {len(markers)} markers, labelled (top {TOP_N} of either comparison, ties included):")
            for _, marker in markers[markers["labelled"]].iterrows():
                print(f"  x={marker['x']:.3f} y={marker['y']:.3f} n={len(marker['drugs']):>2}  "
                      + " ".join(label_text(marker).split()))

            out = os.path.join(drug_dir, "drug_candidates_scatter.png")
            make_figure(table, markers).savefig(out, dpi=300, bbox_inches="tight", pad_inches=0.05)
            plt.close("all")
            print(f"Saved → {out}")
