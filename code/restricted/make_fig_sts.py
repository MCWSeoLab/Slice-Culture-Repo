#!/usr/bin/env python3
"""Staurosporine positive-control dose-response.

Two IC50 dilution series exist in the lab workbooks and no figure has ever shown them. They are
the standard credential for a viability assay: it has dynamic range, and it detects a known
killer of cells in a dose-ordered way. Absence of a positive control is conspicuous in a
feasibility paper.

Both sheets store blank-corrected absorbance per concentration ("T-C" columns). Viability is
expressed against the same sheet's DMSO control, matching the SAP 4.1 definition used
everywhere else.
"""
import numpy as np
import openpyxl
import pandas as pd

from cns_style import apply_style, palette, finalize_axes, save_figure, MM_PER_INCH, plt

apply_style("Nature")
COLS = palette("Nature", n=3)

LIVER = ("/root/.claude/uploads/5ca89d3c-111f-578e-bab9-226ffe80bccc/"
         "272a5f6c-Liver_Bx__Diaphgram_Bx__MTS_assay.xlsx")
PDACF = ("/root/.claude/uploads/5ca89d3c-111f-578e-bab9-226ffe80bccc/"
         "8afd034f-PDAC_Bx__Peritoneum_Bx__MTS_assay.xlsx")


def conc_to_um(s):
    t = str(s).strip().lower().replace("μ", "u")
    if not t.endswith("um"):
        return None
    try:
        return float(t[:-2])
    except ValueError:
        return None


def scrape(path, sheet):
    """Pull (concentration, value_col_1, value_col_2) triples plus the DMSO control row.

    Both sheets lay the summary block out the same way: a label column holding "STS",
    "DMSO Control" or "Untreated", the concentration beside it, then one value per timepoint.
    Found by scanning for those labels rather than by fixed coordinates, because the block
    starts at a different column on each sheet.
    """
    grid = [[c.value for c in r] for r in openpyxl.load_workbook(path, data_only=True)[sheet]
            .iter_rows()]
    doses, ctrl = [], {}
    for r in grid:
        for k, cell in enumerate(r):
            lab = str(cell).strip().lower() if cell is not None else ""
            if lab == "sts":
                um = conc_to_um(r[k + 1]) if k + 1 < len(r) else None
                vals = [v for v in r[k + 2:k + 4] if isinstance(v, (int, float))]
                if um is not None and vals:
                    doses.append((um, vals))
            elif conc_to_um(cell) is not None:
                # The other sheet omits the "STS" prefix and puts the concentration in the
                # label column. Its RAW blocks look the same except that the next cell is a
                # well id, so require the two following cells to both be numeric — that is
                # only true in the summary block.
                nxt = r[k + 1:k + 3]
                if len(nxt) == 2 and all(isinstance(v, (int, float)) for v in nxt):
                    doses.append((conc_to_um(cell), list(nxt)))
            elif lab in ("dmso control", "untreated"):
                vals = [v for v in r[k + 1:k + 3] if isinstance(v, (int, float))]
                if vals:
                    ctrl.setdefault(lab, vals)
    return doses, ctrl


SERIES = [("Liver Bx 3", LIVER, "L.Bx 3 STS IC50"),
          ("PDAC Bx 2", PDACF, "PDAC 2 STS IC-50")]

rows = []
for label, path, sheet in SERIES:
    doses, ctrl = scrape(path, sheet)
    base = ctrl.get("dmso control")
    print(f"\n{label}  ({sheet})")
    print(f"  DMSO control: {base}   untreated: {ctrl.get('untreated')}")
    if not base or not doses:
        print("  !! could not resolve a control or any dose — skipped")
        continue
    seen = set()
    for um, vals in sorted(doses):
        if um in seen:      # both scrape rules can match the same summary row
            continue
        seen.add(um)
        pct = 100 * vals[0] / base[0]
        rows.append(dict(series=label, um=um, pct=pct))
        print(f"    {um:>6.2f} uM  ->  {vals[0]:.4f} / {base[0]:.4f}  =  {pct:6.1f}% of control")

d = pd.DataFrame(rows)
if d.empty:
    raise SystemExit("no dose-response data resolved")

# ---------------------------------------------------------------- single-dose STS
# Staurosporine was ALSO run as a single-concentration positive control on four further
# cohort biopsies. Reading only the two IC50 sheets understates the positive control badly.
t = pd.read_csv("mts_treated_pct_deid.csv")
sd = (t[t.drug == "Staurosporine"]
      .groupby(["sample", "timepoint", "concentration", "qc_fail"])
      .pct_of_control.median().reset_index())
# PDAC Bx 1's STS sits on a day-series sheet the extractor excludes by design; taken from the
# sheet's own day-7 summary, DMSO 1.352 vs STS 0.0235.
sd = pd.concat([sd, pd.DataFrame([dict(sample="PDAC Bx 1", timepoint="day 7",
                                       concentration="10uM", qc_fail=False,
                                       pct_of_control=100 * 0.0235 / 1.352)])],
               ignore_index=True)
# the top dose of each IC50 series is itself a 10 uM single-dose point
for lab in d.series.unique():
    top = d[d.series == lab].sort_values("um").iloc[-1]
    sd = pd.concat([sd, pd.DataFrame([dict(sample=lab, timepoint="48h", concentration="10uM",
                                           qc_fail=False, pct_of_control=top.pct)])],
                   ignore_index=True)
sd["um"] = sd.concentration.str.lower().str.replace("mol", "", regex=False).str.rstrip("um ")
sd["um"] = pd.to_numeric(sd.um, errors="coerce")
print("\nsingle-dose staurosporine, median % of control per biopsy and timepoint")
print(sd.sort_values(["um", "sample"]).to_string(index=False,
      float_format=lambda x: f"{x:.1f}"))

TPX = {"24h": 0, "48h": 1, "72h": 2, "day 7": 3}
fig, (axA, axB) = plt.subplots(1, 2, figsize=(183 / MM_PER_INCH, 74 / MM_PER_INCH))
fig.subplots_adjust(left=0.075, right=0.988, top=0.855, bottom=0.275, wspace=0.28)

# ---- A: dose-response
for i, (label, g) in enumerate(d.groupby("series")):
    g = g.sort_values("um")
    axA.plot(g.um, g.pct, "-o", ms=3.6, lw=1.0, color=COLS[i], label=label, zorder=3)
axA.axhline(100, color="0.35", lw=0.6, ls="--", zorder=1)
axA.axhline(70, color="0.35", lw=0.6, ls=":", zorder=1)
axA.set_xscale("log")
axA.set_xlabel("Staurosporine (\u00b5M, log scale)", fontsize=6.4)
axA.set_ylabel("Ex vivo viability at 48 h\n(% of DMSO control)", fontsize=6.4)
axA.set_title("Dose-response, two biopsies", fontsize=7, pad=5)
axA.set_ylim(-5, max(140, d.pct.max() * 1.12))
axA.text(axA.get_xlim()[1], 101, "DMSO control", fontsize=5.2, ha="right", va="bottom",
         color="0.35")
axA.text(axA.get_xlim()[1], 71, "30% inhibition", fontsize=5.2, ha="right", va="bottom",
         color="0.35")
axA.legend(loc="lower left", frameon=False, fontsize=5.6, handlelength=1.6,
           handletextpad=0.5, borderaxespad=0.5)
finalize_axes(axA, tight=False)

# ---- B: single-dose, all biopsies that had one
CC = {5.0: COLS[1], 10.0: COLS[0]}
LBL = []
for (samp, um), g in sd.groupby(["sample", "um"]):
    g = g.assign(x=g.timepoint.map(TPX)).sort_values("x")
    col = CC.get(um, "0.5")
    axB.plot(g.x, g.pct_of_control, "-", lw=0.9, color=col, zorder=2, alpha=0.9)
    for _, r in g.iterrows():
        axB.plot(r.x, r.pct_of_control, "o", ms=4.0, color=col, zorder=3,
                 markerfacecolor=("white" if r.qc_fail else col))
    LBL.append((g.x.iloc[-1], float(g.pct_of_control.iloc[-1]), samp, col))
# several series end near zero; push labels apart vertically so none sits on another
for k, (x, y, samp, col) in enumerate(sorted(LBL, key=lambda z: (z[0], z[1]))):
    near = [o for o in LBL if o[0] == x and abs(o[1] - y) < 9 and o[2] != samp]
    dy = 4.5 if (near and y <= min(o[1] for o in near)) else (-4.5 if near else 0)
    axB.text(x + 0.09, y + dy, samp, fontsize=4.6, va="center", ha="left", color="0.3")
axB.axhline(100, color="0.35", lw=0.6, ls="--", zorder=1)
axB.axhline(70, color="0.35", lw=0.6, ls=":", zorder=1)
axB.set_xticks(list(TPX.values())); axB.set_xticklabels(list(TPX.keys()), fontsize=6.0)
axB.set_xlim(-0.3, 4.3); axB.set_ylim(-5, 125)
axB.set_xlabel("Timepoint", fontsize=6.4)
axB.set_ylabel("Ex vivo viability (% of DMSO control)", fontsize=6.4)
axB.set_title("Single-dose positive control, six biopsies", fontsize=7, pad=5)
hb = [plt.Line2D([], [], color=COLS[0], marker="o", ms=4, lw=1, label="10 \u00b5M"),
      plt.Line2D([], [], color=COLS[1], marker="o", ms=4, lw=1, label="5 \u00b5M"),
      plt.Line2D([], [], color="0.4", marker="o", ms=4, lw=0, markerfacecolor="white",
                 label="experiment fails the 20% QC gate")]
axB.legend(handles=hb, loc="upper right", frameon=False, fontsize=5.2, handlelength=1.4,
           handletextpad=0.5, borderaxespad=0.5)
finalize_axes(axB, tight=False)

for ax, L in [(axA, "A"), (axB, "B")]:
    ax.text(-0.13, 1.09, L, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")

ten = sd[sd.um == 10]
five = sd[sd.um == 5]
TOP = d[d.um == d.um.max()]
CAP = (
    "Staurosporine, expressed against the same-plate DMSO control.  "
    "A  Six-point dilution series on two biopsies at 48 h; viability falls from above "
    "control at the lowest concentration to "
    f"{TOP.pct.min():.0f}\u2013{TOP.pct.max():.0f}% at 10 \u00b5M.\n"
    "B  Median viability where staurosporine was given at a single concentration, one line "
    "per biopsy; open symbols mark experiments failing the DMSO quality-control gate. "
    f"At 10 \u00b5M the\ncontrol reaches {ten.pct_of_control.min():.0f}\u2013"
    f"{ten.pct_of_control.max():.0f}% of control; at 5 \u00b5M it is weaker and slower, one "
    f"biopsy remaining at {five[five.timepoint == '48h'].pct_of_control.max():.0f}% of "
    "control at 48 h."
)
fig.text(0.075, 0.155, CAP, fontsize=4.9, va="top", linespacing=1.55)

save_figure(fig, "figS_positive_control", outdir="figures")
