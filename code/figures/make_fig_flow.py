#!/usr/bin/env python3
"""REMARK-style flow diagram: biopsies accrued through to the analysed matched pairs.

Every count is read from the analysis files, never typed here — run it and the numbers move
with the data. The SAP states the report follows REMARK, which requires this diagram; without
it a reader has to reconstruct the funnel from prose scattered across the Results.
"""
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from cns_style import apply_style, palette, save_figure, MM_PER_INCH, plt

apply_style("Nature")
GREY = "0.45"

t1 = pd.read_csv("table1_source_deid.csv", comment="#")
p = pd.read_csv("pairs_final.csv")
a = p[~p.flag]
grp = t1.group.value_counts()

N_BX = len(t1)
N_LIVER, N_PERI, N_PDAC = int(grp["Liver"]), int(grp["Peritoneum"]), int(grp["Primary PDAC"])
N_COURSE, N_PT_ALL = len(p), p.study_id.nunique()
# biopsies, not patients — a patient can contribute more than one biopsy, so the two
# counts are not interchangeable and the exclusion box must use biopsies.
N_BX_MATCHED = p.sample_label.nunique()
N_FLAG, N_FLAG_PT = int(p.flag.sum()), p[p.flag].study_id.nunique()
N_PAIR, N_PT = len(a), a.study_id.nunique()
N_PRE = int((a.day < 0).sum())
N_POST, N_POST_PT = int((a.day >= 0).sum()), a[a.day >= 0].study_id.nunique()

print(f"{N_BX} biopsies -> {N_COURSE} courses / {N_PT_ALL} patients -> "
      f"{N_PAIR} pairs / {N_PT} patients -> {N_POST} prospective")

fig = plt.figure(figsize=(120 / MM_PER_INCH, 132 / MM_PER_INCH))
ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)

MAIN_X, MAIN_W = 0.055, 0.60
SIDE_X, SIDE_W = 0.695, 0.285


def box(x, y, w, h, lines, bold_first=True, face="white", edge="0.35", fs=5.6):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.006,rounding_size=0.008",
                                facecolor=face, edgecolor=edge, linewidth=0.7, zorder=2))
    n = len(lines)
    for k, ln in enumerate(lines):
        ax.text(x + w / 2, y + h - (h / (n + 0.6)) * (k + 0.75), ln, ha="center", va="center",
                fontsize=(6.4 if (k == 0 and bold_first) else fs),
                fontweight=("bold" if (k == 0 and bold_first) else "normal"),
                color="black" if k == 0 else "0.25", zorder=3)


def down(y0, y1, x=MAIN_X + MAIN_W / 2):
    ax.add_patch(FancyArrowPatch((x, y0), (x, y1), arrowstyle="-|>", mutation_scale=7,
                                 linewidth=0.8, color="0.35", zorder=1))


def branch(y, x0=MAIN_X + MAIN_W / 2, x1=SIDE_X):
    ax.plot([x0, x1], [y, y], "-", color=GREY, lw=0.7, zorder=1)
    ax.add_patch(FancyArrowPatch((x1 - 0.02, y), (x1, y), arrowstyle="-|>", mutation_scale=6,
                                 linewidth=0.7, color=GREY, zorder=1))


H = 0.098
rows = [
    (0.828, [f"Biopsies cultured as organotypic slices (n = {N_BX})",
             f"{N_LIVER} liver, {N_PERI} peritoneal, {N_PDAC} primary pancreatic"]),
    (0.632, [f"Courses matched to an ex vivo result (n = {N_COURSE})",
             f"from {N_PT_ALL} patients",
             "same agent or agent class, scorable RECIST best response"]),
    (0.436, [f"Analysis set (n = {N_PAIR} drug–patient pairs)",
             f"from {N_PT} patients"]),
    (0.240, [f"Prospective subset (n = {N_POST})",
             f"from {N_POST_PT} patients",
             "treatment begun after the biopsy"]),
]
for y, lines in rows:
    box(MAIN_X, y, MAIN_W, H, lines)
for i in range(len(rows) - 1):
    down(rows[i][0], rows[i + 1][0] + H)

SIDE = [
    (0.668, [f"No matched clinical course",
             f"{N_BX - N_BX_MATCHED} of {N_BX} biopsies",
             "no scorable response to a tested agent"]),
    (0.472, [f"Tissue-exhausted specimen",
             f"{N_FLAG} courses, {N_FLAG_PT} patient",
             "excluded from all statistics; shown in grey"]),
    (0.276, [f"Course began before the biopsy",
             f"{N_PRE} of {N_PAIR} pairs",
             "assay could not have informed the decision"]),
]
for y, lines in SIDE:
    box(SIDE_X, y, SIDE_W, 0.086, lines, bold_first=False, face="0.965", edge=GREY, fs=5.0)
    branch(y + 0.043)

ax.text(MAIN_X, 0.985, "Flow of biopsies through the analysis", fontsize=8,
        fontweight="bold", va="top", ha="left")
ax.text(MAIN_X, 0.955,
        "Boxes on the right are exclusions and qualifications, not further groups.",
        fontsize=5.2, va="top", ha="left", color="0.35")

ax.text(MAIN_X, 0.196,
        f"Of the {N_PAIR} analysed pairs, {N_PRE} began before the biopsy and {N_POST} after. "
        f"The {N_POST} prospective\npairs, from {N_POST_PT} patients, are the only subset in "
        f"which the assay result preceded the\ntreatment decision, and are the effect size a "
        f"prospective study should be designed around.",
        fontsize=5.2, va="top", ha="left", color="0.25", linespacing=1.55)

save_figure(fig, "fig_flow_diagram", outdir="figures")
