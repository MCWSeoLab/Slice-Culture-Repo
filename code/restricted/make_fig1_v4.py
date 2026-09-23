#!/usr/bin/env python3
"""Figure 1 v4 — v3 with panel B brought to a single magnification.

Only change from v3: the day-3 and day-5 H&E tiles are now 250 um fields, matching day 1
and day 2.  Nothing else about the figure moves.


Replaces the Illustrator artboard `SRC_core_workflow_viability`. The schematic is kept as
artwork (imported as one image); everything else is redrawn so Figure 1 matches Figures 2-4
in font, palette, axis weight and panel-label convention.

  A  core biopsy to plated slice to readout — workflow schematic
  B  H&E across seven days in culture
  C  immunophenotype at day 3
  D  slice surface area over seven days
  E  MTS absorbance over seven days

Changes from the artboard, all deliberate:
  * panel letters added (there were none)
  * the MTS-vs-nuclear-density scatter is dropped — Figure 2A already carries it, with a
    bootstrap CI and honest n, so keeping it here duplicated the weaker version
  * D and E now show one line per biopsy instead of a mean-only bar/line, because the three
    biopsies have unequal follow-up and the apparent day-3 dip is one biopsy, not a trend
  * rainbow fills removed; colour no longer implies grouping that does not exist
  * both panels state the ANOVA over days, which is the actual feasibility claim
"""
import numpy as np, pandas as pd
from scipy import stats
import matplotlib.image as mpimg
from matplotlib.gridspec import GridSpec
from cns_style import apply_style, palette, finalize_axes, save_figure, MM_PER_INCH, plt

apply_style("Nature")

TC = pd.read_csv("fig1_timecourse_deid.csv", comment="#")
BXCOL = dict(zip(["Bx1", "Bx2", "Bx3"], palette("Nature", n=3)))
MEAN = "0.20"

W, H = 183 / MM_PER_INCH, 206 / MM_PER_INCH
fig = plt.figure(figsize=(W, H))
# Panels B and C enlarged 4 Sep 2026 - same magnification, more page.  The bands were
# height-bound at 0.40, which capped the tiles at 21 mm; 0.60 lets the width bind
# instead.  27 mm is the ceiling: the kept day-3/day-5 and IHC tiles are 320 px
# originals, and past 27 mm they fall below 300 ppi.
gs = GridSpec(4, 2, figure=fig, height_ratios=[1.80, 0.50, 0.50, 1.02],
              hspace=0.36, wspace=0.30, left=0.095, right=0.975, top=0.985, bottom=0.065)


def letter(x, y, L):
    fig.text(x, y, L, fontsize=9, fontweight="bold", va="top", ha="left")


# Field width in micrometres for every tile in panels B and C.  As of v4 ALL FOUR panel B
# tiles are 250 um fields, so the H&E row is at one magnification - which is what changed
# from v3.  Day 1 and day 2 are digital crops from the whole-slide scans at 0.5 um/px; day 3
# and day 5 are 1471 px crops taken from the ORIGINAL 2560 x 1920 microscope captures at 20x
# (0.17 um/px, so 2560 px = a 435 um frame and 1471 px = 250 um).  Cropping the source frame
# rather than the 320 px figure tile means this is a real optical zoom, not an upscale.
#
# The 250 um whole-slide crops of `Liver BX2 A3B3 Day 3` and `A4B4 Day 4` are NOT usable for
# this - the day-3 section is out of focus across the whole slide and the day-4 section is
# acellular (see wsi-field-selection).  The deck TIFFs are the same tissue as the tiles they
# replace.
#
# Panel C is unchanged: the four IHC tiles are still the original microscope captures resized
# from a 2560 px frame at 10x, 320 px x 2.72 um/px = an 870 um field.  So panel C alone is at
# a wider field width, and every tile still carries its own bar.
FIELD_UM = {
    "HE_Day1": 250.0, "HE_Day2": 250.0,
    "HE_Day3": 250.0, "HE_Day5": 250.0,
    "IHC_CC3": 320 * 2.72, "IHC_Ki67": 320 * 2.72,
    "IHC_CD45": 320 * 2.72, "IHC_EpCAM": 320 * 2.72,
}


def scale_bar(ax, stem):
    fum = FIELD_UM.get(stem)
    if not fum:
        return
    bar = 50.0 if fum < 600.0 else 100.0
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()          # imshow inverts y
    w, h = abs(x1 - x0), abs(y0 - y1)
    bx1 = min(x0, x1) + w * 0.955
    bx0 = bx1 - w * (bar / fum)
    by = max(y0, y1) - h * 0.085
    ax.plot([bx0, bx1], [by, by], color="black", lw=1.3, solid_capstyle="butt",
            clip_on=False, zorder=6)
    ax.text((bx0 + bx1) / 2, by - h * 0.04, f"{bar:.0f} \u00b5m", color="black",
            fontsize=4.6, ha="center", va="bottom", zorder=6)


def tile_row(cell, items, folder, x0, x1, titles=True, fs=6.0):
    """Lay a row of micrographs out by hand so they fill the band at one common size."""
    box = cell.get_position(fig)
    n = len(items)
    gap = 0.007
    tw = (x1 - x0 - gap * (n - 1)) / n
    imh, imw = mpimg.imread(f"{folder}/{items[0][1]}.png").shape[:2]
    th = tw * (imh / imw) * (W / H)
    if th > box.height:
        th = box.height; tw = th / ((imh / imw) * (W / H))
    bx = x0 + ((x1 - x0) - (tw * n + gap * (n - 1))) / 2
    by = box.y0 + (box.height - th) / 2
    for k, (lab, stem) in enumerate(items):
        b = fig.add_axes([bx + k * (tw + gap), by, tw, th])
        b.imshow(mpimg.imread(f"{folder}/{stem}.png"))
        b.set_xticks([]); b.set_yticks([])
        for sp in b.spines.values():
            sp.set_visible(True); sp.set_linewidth(0.4); sp.set_color("0.45")
        if titles:
            b.set_title(lab, fontsize=fs, pad=2)
        scale_bar(b, stem)
    return bx, by, th


# ---- A: workflow schematic, kept as artwork
cellA = gs[0, :].get_position(fig)
sch = mpimg.imread("fig1_images/schematic.png")
sh = (sch.shape[0] / sch.shape[1]) * (cellA.x1 - cellA.x0) * (W / H)
axA = fig.add_axes([cellA.x0, cellA.y1 - sh, cellA.x1 - cellA.x0, sh])
axA.imshow(sch); axA.axis("off")
letter(0.012, cellA.y1 + 0.004, "A")

# ---- B: H&E across seven days
bx, by, th = tile_row(gs[1, :], [("Day 1", "HE_Day1"), ("Day 2", "HE_Day2"),
                                 ("Day 3", "HE_Day3"), ("Day 5", "HE_Day5")],
                      "fig1_images", 0.115, 0.885)
fig.text(bx - 0.018, by + th / 2, "H&E", fontsize=6.4, rotation=90, va="center", ha="right")
letter(0.012, by + th + 0.020, "B")

# ---- C: immunophenotype at day 3
bx, by, th = tile_row(gs[2, :], [("Cleaved casp-3", "IHC_CC3"), ("Ki-67", "IHC_Ki67"),
                                 ("CD45", "IHC_CD45"), ("EpCAM", "IHC_EpCAM")],
                      "fig1_images", 0.115, 0.885)
fig.text(bx - 0.018, by + th / 2, "Day 3", fontsize=6.4, rotation=90, va="center", ha="right")
letter(0.012, by + th + 0.020, "C")


def timecourse(ax, key, ylab, note):
    d = TC[TC.panel == key]
    for b, g in d.groupby("biopsy"):
        g = g.sort_values("day")
        ax.plot(g.day, g.value, "-o", ms=3.0, lw=0.8, color=BXCOL[b],
                markerfacecolor="white", markeredgecolor=BXCOL[b], markeredgewidth=0.9,
                label=b, zorder=3)
    m = d.groupby("day").value.agg(["mean", "count"]).reset_index()
    ax.plot(m.day, m["mean"], "-", lw=1.6, color=MEAN, zorder=4)
    ax.plot(m.day, m["mean"], "s", ms=3.6, color=MEAN, zorder=5)
    groups = [g.value.values for _, g in d.groupby("day")]
    F, pv = stats.f_oneway(*[g for g in groups if len(g) > 0])
    dfn = len(groups) - 1
    dfd = int(sum(len(g) for g in groups)) - len(groups)
    ax.set_xticks(sorted(d.day.unique()))
    ax.set_xlabel("Day in culture", fontsize=6.4)
    ax.set_ylabel(ylab, fontsize=6.4)
    ax.text(0.02, 0.025,
            f"One-way ANOVA over days: F({dfn},{dfd}) = {F:.2f}, p = {pv:.2f}\n"
            f"n = {', '.join(str(int(c)) for c in m['count'])} biopsies at day "
            f"{', '.join(str(int(x)) for x in m.day)}\n" + note,
            transform=ax.transAxes, va="bottom", fontsize=5.2, linespacing=1.45)
    finalize_axes(ax, tight=False)


# ---- D / E: the two seven-day series
axD = fig.add_subplot(gs[3, 0]); axE = fig.add_subplot(gs[3, 1])
timecourse(axD, "surface", "Slice surface area",
           "black = mean across the biopsies measured that day")
timecourse(axE, "mts", "MTS absorbance (490 nm)",
           "black = mean across the biopsies measured that day")
axD.set_ylim(0, None); axE.set_ylim(0, None)
axD.legend(loc="upper left", frameon=False, fontsize=5.4, handletextpad=0.5,
           borderaxespad=0.2, labelspacing=0.25,
           title="core biopsy", title_fontsize=5.4)
pos = axD.get_position(); letter(0.012, pos.y1 + 0.022, "D")
pos = axE.get_position(); letter(pos.x0 - 0.083, pos.y1 + 0.022, "E")

fig.text(0.095, 0.030,
         "B, C: representative slices; scale bars as marked. "
         "D, E: three liver core biopsies with unequal follow-up — neither surface area nor "
         "MTS changes significantly across seven days in culture.",
         fontsize=5.4, va="top")

save_figure(fig, "fig1_workflow_feasibility_v4", outdir="figures")
for k in ["surface", "mts"]:
    d = TC[TC.panel == k]
    print(k, d.groupby("day").value.agg(["mean", "count"]).round(3).to_dict("index"))
