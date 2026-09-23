#!/usr/bin/env python3
"""ERK figure in the Figure-4 idiom: two cases side by side, with the ERK readouts
appended as their own bands.

Cases here are the two patients who actually have ERK staining — PT-014 (KRAS G12D)
and PT-025 (KRAS G12R) — so unlike the merged Figure 4, the alleles in the ERK bands
belong to the same columns as the MTS and IHC above them.

  A / E  ex vivo viability at 48 h by agent
  B / F  representative stains, DMSO vs RMC-7977 — the four IHC markers and the two
         ERK stains in one band (the ERK stains were their own panel, F / L, before)
  C / G  IHC fold change vs DMSO
  D / H  total ERK, phospho-ERK, p-ERK/T-ERK ratio (DMSO vs RMC-7977)

The clinical-course panels (D / J) were removed on request; the courses for these
patients are shown in Figure 3C and, for PT-015, in Figure 4.
"""
import pandas as pd, numpy as np
import matplotlib.image as mpimg
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from course_panel import draw_course
from cns_style import apply_style, palette, finalize_axes, save_figure, MM_PER_INCH, plt

apply_style("Nature")

CASES = [dict(sid="PT-014", ihc="Liver Bx.6", allele="KRAS G12D", img="G12D",
              title="Case A — KRAS G12D liver metastasis", key="Gem/nab-Pac"),
         dict(sid="PT-025", ihc="Liver Bx. [DATE]", allele="KRAS G12R", img="G12R",
              title="Case B — KRAS G12R liver metastasis", key="Gem/nab-Pac")]

DRUGCOL = dict(zip(["FOLFIRINOX", "Gem/nab-Pac", "RMC-7977", "Trametinib", "Vemurafenib"],
                   palette("Nature", n=5)))
RESPCOL = {"CR": "0.20", "PR": "0.35", "SD": "0.62", "PD": "0.88", "": "1.0"}
ARMS, ACOL = ["DMSO", "RMC-7977"], {"DMSO": "0.62", "RMC-7977": DRUGCOL["RMC-7977"]}
MARKERS = ["Cleaved caspase-3", "Ki-67", "CD45", "EpCAM"]


def canon(s):
    s = str(s).lower().replace(" ", "").replace("-", "").replace("/", "")
    if "folfirinox" in s: return "FOLFIRINOX"
    if "gem" in s and ("nab" in s or "pac" in s): return "Gem/nab-Pac"
    if "rmc" in s: return "RMC-7977"
    if "trametinib" in s: return "Trametinib"
    if "vemuraf" in s: return "Vemurafenib"
    return None


mllink = pd.read_csv("masterlist_linked_deid.csv"); mllink["drug"] = mllink.treatment.map(canon)
ihc = pd.read_csv("ihc_tidy.csv")
gg = ihc.groupby(["biopsy", "marker", "treatment"]).value.mean().reset_index()
ref = gg[gg.treatment == "DMSO"].set_index(["biopsy", "marker"]).value.rename("dmso")
gg = gg.join(ref, on=["biopsy", "marker"])
gg = gg[(gg.treatment != "DMSO") & gg.dmso.gt(0)].copy()
gg["fold"] = gg.value / gg.dmso; gg["drug"] = gg.treatment.map(canon)

rg = pd.read_csv("regimens_merged.csv").fillna("")
for c in ["day_start_rel_biopsy", "day_end_rel_biopsy"]:
    rg[c] = pd.to_numeric(rg[c], errors="coerce")
EXCL = r"whipple|\bbx\b|biopsy|sbrt|chemort|cgy|\bfx\b|resect|NED|LAB"

erk = pd.read_csv("erk_tidy.csv")
ew = erk.pivot_table(index=["allele", "treatment", "replicate"],
                     columns="marker", values="value").reset_index()
ew["ratio"] = ew["p-ERK"] / ew["T-ERK"]

W, H = 183 / MM_PER_INCH, 244 / MM_PER_INCH
fig = plt.figure(figsize=(W, H))
# Row 1 carries a 2 x 6 stain band now (4 IHC + 2 ERK), so the tiles are narrower than
# the old 2 x 4 band and the row needs proportionally less height.
# Row 1 holds six 4:3 tiles across half the width, so its content is only ~24 mm tall;
# give it a matching share or the band floats in a tall empty cell.
# The stain bands were 6 tiles across a HALF-page column, which capped the tiles at 11 mm
# and left the row floating in a tall empty cell.  Giving each case's band the full page
# width, stacked, converts that dead space into tile size.
gs = GridSpec(5, 2, figure=fig, height_ratios=[1.00, 1.25, 1.25, 1.00, 0.96],
              hspace=0.58, wspace=0.13, left=0.070, right=0.992, top=0.950, bottom=0.105)
letters = iter("ABCDEFGH")

for col, C in enumerate(CASES):
    sid = C["sid"]

    # ---- ex vivo MTS
    a = fig.add_subplot(gs[0, col])
    d48 = mllink[(mllink.study_id == sid) & (mllink.timepoint_hr == 48) & mllink.drug.notna()]
    order = [x for x in DRUGCOL if x in set(d48.drug)]
    for i, dr in enumerate(order):
        v = d48[d48.drug == dr].mts_pct_of_control.mean()
        a.bar(i, v, width=0.64, color=DRUGCOL[dr], edgecolor="black", linewidth=0.5,
              zorder=2, hatch="//" if dr == C["key"] else None)
        a.text(i, v + 3, f"{v:.0f}", ha="center", va="bottom", fontsize=5.6)
    a.axhline(100, color="0.35", lw=0.5, ls="--", zorder=1)
    a.axhline(70, color="0.35", lw=0.5, ls=":", zorder=1)
    a.set_xticks(range(len(order)))
    a.set_xticklabels(order, rotation=28, ha="right", fontsize=5.6)
    a.set_xlim(-0.62, len(order) - 0.38); a.set_ylim(0, 150)
    a.set_ylabel("Ex vivo viability at 48 h\n(% of DMSO control)", fontsize=6.2)
    a.set_title(C["title"], fontsize=7, pad=6)
    if col == 1:
        a.text(len(order) - 0.4, 103, "DMSO control", fontsize=5.2, ha="right",
               va="bottom", color="0.35")
        a.text(len(order) - 0.4, 73, "30% inhibition", fontsize=5.2, ha="right",
               va="bottom", color="0.35")
    finalize_axes(a)
    a.text(-0.20, 1.22, next(letters), transform=a.transAxes, fontsize=9,
           fontweight="bold", va="top")

    # ---- representative stains, DMSO vs RMC-7977: the four IHC markers followed by the
    # two ERK stains, which used to be a separate panel. All six are the same magnification
    # and the same two arms, so one band reads better than two.
    STAINS = [("ihc_images", "CC3", "Cleaved casp-3"), ("ihc_images", "Ki67", "Ki-67"),
              ("ihc_images", "CD45", "CD45"), ("ihc_images", "EpCAM", "EpCAM"),
              ("erk_images", "TERK", "T-ERK"), ("erk_images", "pERK", "p-ERK")]
    # A narrow empty column between EpCAM and T-ERK separates the IHC markers from the ERK
    # stains, so the band reads as two groups without needing a rule drawn across it.
    SPACER = 4
    wr = [1.0] * (len(STAINS) + 1); wr[SPACER] = 0.16
    band_cell = gs[1 + col, :]
    gsm = band_cell.subgridspec(2, len(STAINS) + 1, width_ratios=wr,
                                 wspace=0.035, hspace=0.035)
    for rr, arm in enumerate(["DMSO", "RMC7977"]):
        for cc, (folder, mk, lab) in enumerate(STAINS):
            b = fig.add_subplot(gsm[rr, cc if cc < SPACER else cc + 1])
            b.imshow(mpimg.imread(f"{folder}/{C['img']}_{arm}_{mk}.png"))
            b.set_xticks([]); b.set_yticks([])
            for sp in b.spines.values():
                sp.set_linewidth(0.4); sp.set_color("0.4")
            if rr == 0:
                b.set_title(lab, fontsize=5.0, pad=2)
            if cc == 0:
                b.set_ylabel({"DMSO": "DMSO", "RMC7977": "RMC-7977"}[arm], fontsize=5.4)
    ax_m = fig.add_subplot(band_cell); ax_m.axis("off")
    ax_m.text(-0.052, 1.10, next(letters), transform=ax_m.transAxes, fontsize=9,
              fontweight="bold", va="top")

    # ---- IHC fold change
    a = fig.add_subplot(gs[3, col])
    sub = gg[gg.biopsy == C["ihc"]]
    drugs = [x for x in order if x in set(sub.drug)]
    wdt = 0.8 / max(len(drugs), 1)
    for j, dr in enumerate(drugs):
        xs = np.arange(len(MARKERS)) + (j - (len(drugs) - 1) / 2) * wdt
        vals = [sub[(sub.drug == dr) & (sub.marker == m)].fold.mean() for m in MARKERS]
        for mi, mv in enumerate(vals):
            if not np.isfinite(mv) and j == 0:
                a.text(mi, 0.02, "n.d.", ha="center", va="bottom", fontsize=4.8,
                       color="0.45", transform=a.get_xaxis_transform())
        a.bar(xs, vals, width=wdt * 0.9, color=DRUGCOL[dr], edgecolor="black",
              linewidth=0.4, zorder=2, hatch="//" if dr == C["key"] else None)
    a.axhline(1.0, color="0.35", lw=0.5, ls="--", zorder=1)
    a.set_xticks(range(len(MARKERS)))
    a.set_xticklabels(MARKERS, rotation=22, ha="right", fontsize=5.6)
    a.set_ylabel("IHC fold change\nvs DMSO", fontsize=6.2)
    if not len(sub):
        a.text(0.5, 0.5, "no IHC recorded for this biopsy", transform=a.transAxes,
               ha="center", va="center", fontsize=5.6, color="0.45")
    finalize_axes(a)
    a.text(-0.20, 1.16, next(letters), transform=a.transAxes, fontsize=9,
           fontweight="bold", va="top")

    # ---- ERK quantification for THIS patient's allele
    gsq = gs[4, col].subgridspec(1, 3, wspace=0.75)
    for k, (meas, ylab) in enumerate([("T-ERK", "Total ERK\n(% positive)"),
                                      ("p-ERK", "Phospho-ERK\n(% positive)"),
                                      ("ratio", "p-ERK / T-ERK")]):
        b = fig.add_subplot(gsq[k])
        for j, arm in enumerate(ARMS):
            v = ew[(ew.allele == C["allele"]) & (ew.treatment == arm)][meas].values
            b.bar(j, np.mean(v), width=0.5, color=ACOL[arm], alpha=0.40,
                  edgecolor="none", zorder=1)
            b.plot([j - 0.25, j + 0.25], [np.mean(v)] * 2, color=ACOL[arm], lw=1.0, zorder=2)
            b.plot(np.full(len(v), j) + np.linspace(-0.10, 0.10, len(v)), v, "o", ms=2.8,
                   markerfacecolor=ACOL[arm], markeredgecolor="black",
                   markeredgewidth=0.35, zorder=3)
        b.set_xticks([0, 1]); b.set_xticklabels(ARMS, rotation=20, ha="right", fontsize=5.4)
        b.set_xlim(-0.6, 1.6)
        b.set_ylabel(ylab, fontsize=5.9)
        finalize_axes(b)
    ax_q = fig.add_subplot(gs[4, col]); ax_q.axis("off")
    ax_q.text(-0.115, 1.16, next(letters), transform=ax_q.transAxes, fontsize=9,
              fontweight="bold", va="top")


h1 = [Line2D([], [], marker="s", ls="none", ms=5, mfc=DRUGCOL[d], mec="black", mew=0.4, label=d)
      for d in ["FOLFIRINOX", "Gem/nab-Pac", "RMC-7977", "Trametinib"]]
h2 = [Line2D([], [], marker='s', ls='none', ms=5, mfc='white', mec='black', mew=0.6,
             label='Hatched = agent the patient received'),
      Line2D([], [], ls='--', color='0.35', lw=0.6, label='DMSO control (100%)'),
      Line2D([], [], ls=':', color='0.35', lw=0.6, label='30% inhibition threshold')]
fig.legend(handles=h1 + h2, loc="lower center", bbox_to_anchor=(0.5, 0.006), frameon=False,
           fontsize=5.8, ncol=4, handletextpad=0.5, columnspacing=1.3)
fig.text(0.085, 0.060,
         "Stain bands: the four leftmost columns are the IHC markers, the two to the right "
         "of the gap are total and phospho-ERK.\n"
         "ERK panels: each point is one slice (4 technical replicates per condition). "
         "Images 10x.",
         fontsize=5.2, va="top", linespacing=1.5)

save_figure(fig, "figS_erk_cases_wide", outdir="figures")
print("  cases:", [c["sid"] for c in CASES])
