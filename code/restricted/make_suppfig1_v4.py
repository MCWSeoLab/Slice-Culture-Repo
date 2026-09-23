#!/usr/bin/env python3
"""Supplementary Figure 1 - primary PDAC and peritoneal slices in extended culture.

A  primary PDAC, MTS absorbance, days 1-7
B  primary PDAC, nuclear density on H&E, days 0-7
C  primary PDAC, H&E micrographs, days 2-7
D  primary PDAC, four-marker IHC micrographs, day 2 and day 7
E  primary PDAC, IHC quantification, days 0-7
F  metastatic peritoneal core, H&E at day 3 and day 7, 4x and 10x

Data: deid/suppfig1_pdac_*.csv, built and cross-checked by build_supp_data.py.
Panels A and B carry NO significance test. Both series come from a single specimen with
technical replicates (2 slices per day in A, 6 fields per day in B), so a one-way ANOVA across
days treats technical replication as biological and its p-value does not support a claim about
the platform. The Prism F statistics that used to be printed here - F(3, 4) = 13.10, p = 0.0155
and F(4, 25) = 1.23, p = 0.32 - were hard-coded strings; both reproduce exactly from the deid
CSVs, and both were removed. Panel B's was also being read as evidence that density was
"unchanged", which is accepting the null. The panels state the design only
("One specimen, N ... per day"); the figure caption is what states that no test was applied.

Panel F was previously a separate supplementary figure. That figure's drug-exposure
panels (MTS, nuclear density and IHC by agent) were cut: the experiment used
staurosporine at 5 uM, which the project's own positive-control dose-response
shows is not an adequate control in this assay - 10 uM gives 0.6 to 14.8% of
control across four biopsies, while 5 uM left this core at 78%. With an
underdosed positive control and no separation between agents on any endpoint the
panels were not interpretable. Their source data stay in
deid/suppfig2_peritoneal_*.csv.

Laid out in explicit millimetres rather than a GridSpec: micrograph tiles are
4:3 and imshow preserves aspect, so a row sized by ratio leaves the tiles
floating in a taller cell.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cns_style import MM_PER_INCH, apply_style, palette, save_figure  # noqa: E402
import imgprep  # noqa: E402

DEID = Path("deid")
OUT = Path("figures")

apply_style(base_fontsize=7, small_fontsize=6, linewidth=0.5)
COL = palette("Nature")

mts = pd.read_csv(DEID / "suppfig1_pdac_mts_deid.csv")
nd = pd.read_csv(DEID / "suppfig1_pdac_nuclear_density_deid.csv")
ihc = pd.read_csv(DEID / "suppfig1_pdac_ihc_deid.csv")

# ------------------------------------------------------------------ layout
W, H = 183.0, 266.0
L, R = 13.0, 3.0
CW = W - L - R                      # content width
TILE_GAP = 3.0
TILE_W = (CW - 3 * TILE_GAP) / 4
TILE_H = TILE_W * 3 / 4
TITLE_H = 4.0

fig = plt.figure(figsize=(W / MM_PER_INCH, H / MM_PER_INCH))


def ax_mm(x, y_top, w, h):
    """Axes placed in millimetres from the top-left of the canvas."""
    return fig.add_axes([x / W, 1 - (y_top + h) / H, w / W, h / H])


def panel_letter(x_mm, y_mm, letter):
    fig.text(x_mm / W, 1 - y_mm / H, letter, fontsize=9, fontweight="bold",
             va="top", ha="left")


def row_label(y_mm, text, size=6.5):
    fig.text((L - 2.5) / W, 1 - y_mm / H, text, fontsize=size, rotation=90,
             va="center", ha="center")


def strip(ax, df, order, color, size=8, jitter=0.13):
    """Mean bar with SD whisker and every replicate shown."""
    rng = np.random.default_rng(20260818)
    for i, g in enumerate(order):
        v = df.loc[df.group == g, "value"].to_numpy()
        ax.bar(i, v.mean(), width=0.62, color=color, alpha=0.28,
               edgecolor=color, linewidth=0.5, zorder=1)
        if len(v) > 1:
            ax.errorbar(i, v.mean(), yerr=v.std(ddof=1), color=color,
                        linewidth=0.6, capsize=1.8, capthick=0.6, zorder=2)
        ax.scatter(i + rng.uniform(-jitter, jitter, len(v)), v, s=size,
                   color=color, edgecolor="white", linewidth=0.25, zorder=3,
                   clip_on=False)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)
    ax.set_xlim(-0.65, len(order) - 0.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(which="both", direction="out")


DAYS_MTS = ["Day 1", "Day 3", "Day 5", "Day 7"]
DAYS_ALL = ["Day 0", "Day 2", "Day 3", "Day 5", "Day 7"]

y = 8.0
PLOT_H = 42.0
half = (CW - 17.0) / 2

# ---------------------------------------------------------------- A  MTS
axA = ax_mm(L, y, half, PLOT_H)
strip(axA, mts, DAYS_MTS, COL[3], size=10)
axA.set_ylabel("MTS absorbance (490 nm)")
axA.set_title("Primary PDAC: metabolic activity", fontsize=7)
axA.set_ylim(0, 1.45)
axA.text(0.03, 0.97, "One specimen, 2 slices per day",
         transform=axA.transAxes, fontsize=5.4, va="top", ha="left",
         linespacing=1.35)
panel_letter(1.5, y - 2.5, "A")

# --------------------------------------------------- B  nuclear density
axB = ax_mm(L + half + 17.0, y, half, PLOT_H)
strip(axB, nd, DAYS_ALL, COL[2], size=7)
axB.set_ylabel("Nuclear density (nuclei/mm$^2$)")
axB.set_title("Primary PDAC: cellularity on H&E", fontsize=7)
axB.set_ylim(0, 3300)
axB.text(0.03, 0.97, "One specimen, 6 fields per day",
         transform=axB.transAxes, fontsize=5.4, va="top", ha="left",
         linespacing=1.35)
panel_letter(L + half + 6.0, y - 2.5, "B")

# ------------------------------------------------------ C  H&E micrographs
y = 58.0
panel_letter(1.5, y - 1.0, "C")
row_label(y + TITLE_H + TILE_H / 2, "Primary PDAC, H&E")
for i, (name, lab) in enumerate(
        [("PDAC_HE_Day2", "Day 2"), ("PDAC_HE_Day3", "Day 3"),
         ("PDAC_HE_Day5", "Day 5"), ("PDAC_HE_Day7", "Day 7")]):
    ax = ax_mm(L + i * (TILE_W + TILE_GAP), y + TITLE_H, TILE_W, TILE_H)
    imgprep.show(ax, name, lab)

# ------------------------------------------------------ D  IHC micrographs
y = 100.0
panel_letter(1.5, y - 1.0, "D")
MARKERS = [("CC3", "Cleaved caspase-3"), ("Ki67", "Ki-67"),
           ("CD45", "CD45"), ("EpCAM", "EpCAM")]
for row, (day, lab) in enumerate([("D2", "Day 2"), ("D7", "Day 7")]):
    y_row = y + TITLE_H + row * (TILE_H + 2.5)
    row_label(y_row + TILE_H / 2, lab)
    for i, (key, mlab) in enumerate(MARKERS):
        ax = ax_mm(L + i * (TILE_W + TILE_GAP), y_row, TILE_W, TILE_H)
        imgprep.show(ax, f"PDAC_{day}_{key}", mlab if row == 0 else None)

# --------------------------------------------------- E  IHC quantification
y = 174.0
axE = ax_mm(L, y, CW, 44.0)
ORDER = ["Cleaved caspase-3", "Ki-67", "CD45", "EpCAM"]
mcol = dict(zip(ORDER, [COL[0], COL[1], COL[2], COL[5]]))
rng = np.random.default_rng(20260818)
width = 0.19
for j, m in enumerate(ORDER):
    sub = ihc[ihc.marker == m]
    for i, d in enumerate(DAYS_ALL):
        v = sub.loc[sub.group == d, "value"].to_numpy()
        x = i + (j - 1.5) * width
        axE.bar(x, v.mean(), width=width * 0.88, color=mcol[m], alpha=0.28,
                edgecolor=mcol[m], linewidth=0.45, zorder=1,
                label=m if i == 0 else None)
        if len(v) > 1:
            axE.errorbar(x, v.mean(), yerr=v.std(ddof=1), color=mcol[m],
                         linewidth=0.5, capsize=1.2, capthick=0.5, zorder=2)
        axE.scatter(x + rng.uniform(-0.035, 0.035, len(v)), v, s=3.5,
                    color=mcol[m], edgecolor="white", linewidth=0.15, zorder=3,
                    clip_on=False)
axE.set_xticks(range(len(DAYS_ALL)))
axE.set_xticklabels(DAYS_ALL)
axE.set_xlim(-0.6, len(DAYS_ALL) - 0.4)
axE.set_ylabel("Positive cells (% of total)")
axE.set_title("Primary PDAC: immunohistochemistry in culture", fontsize=7)
axE.set_ylim(0, 45)
axE.spines["top"].set_visible(False)
axE.spines["right"].set_visible(False)
axE.tick_params(which="both", direction="out")
axE.legend(frameon=False, fontsize=6, ncol=4, loc="upper left",
           handlelength=1.0, handletextpad=0.5, columnspacing=1.4,
           borderaxespad=0.15)
axE.text(0.995, 0.97, "4 fields per marker per day", transform=axE.transAxes,
         fontsize=5.4, va="top", ha="right")
panel_letter(1.5, y - 2.5, "E")

# ------------------------------------- F  peritoneal core H&E (was Supp Fig 2)
y = 226.0
panel_letter(1.5, y - 1.0, "F")
row_label(y + TITLE_H + TILE_H / 2, "Peritoneal core, H&E")
# These are digital fields cropped from a 20x whole-slide scan, not captures through a 4x or
# a 10x objective, and after the 4 Sep 2026 rematch to Figure 4 the detail tiles are 250 um
# fields - about 2.5x tighter than the old "10x" label implied.  Label them by what they are
# and let the scale bars carry the magnification.
PERI = [("PERI_HE_Day3_4x", "Day 3, overview"), ("PERI_HE_Day3_10x", "Day 3, detail"),
        ("PERI_HE_Day7_4x", "Day 7, overview"), ("PERI_HE_Day7_10x", "Day 7, detail")]
for i, (name, lab) in enumerate(PERI):
    ax = ax_mm(L + i * (TILE_W + TILE_GAP), y + TITLE_H, TILE_W, TILE_H)
    imgprep.show(ax, name, lab)

OUT.mkdir(parents=True, exist_ok=True)
save_figure(fig, "figS1_extended_culture_v4", outdir=OUT, formats=("png",), dpi=500)
# Vector output: matplotlib resamples embedded rasters to savefig.dpi, which the
# house style sets to 600. The micrograph tiles print ~40 mm wide, so 350 dpi is
# already above the 300 dpi journals ask for and keeps the SVG under 20 MB.
import matplotlib as mpl  # noqa: E402

with mpl.rc_context({"savefig.dpi": 350}):
    save_figure(fig, "figS1_extended_culture_v4", outdir=OUT, formats=("pdf", "svg"))
