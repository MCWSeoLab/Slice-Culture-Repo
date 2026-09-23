#!/usr/bin/env python3
"""Do two slices of the same condition give the same answer?

A "condition" is one drug, on one biopsy, at one timepoint. Most were cultured as several
slices, and the analysis uses the MEAN of those slices. This asks how much the individual
slices move around that mean, and whether a single slice would have given the same
sensitive/resistant call.

The headline metric is deliberately the mild one: the fraction of individual slices whose own
call matches their condition's mean call. An earlier version counted any condition whose slices
straddled the boundary, which scored a condition with eight slices and one marginal outlier the
same as a genuine 50/50 split, and overstated the problem.
"""
import numpy as np
import pandas as pd

from cns_style import apply_style, palette, finalize_axes, save_figure, MM_PER_INCH, plt

apply_style("Nature")
COLS = palette("Nature", n=3)
CUT = 70.0

t = pd.read_csv("mts_treated_pct_deid.csv")
t["call"] = t.pct_of_control < CUT

cond = (t.groupby(["sample", "timepoint", "drug"])
        .pct_of_control.agg(n="size", mean="mean", lo="min", hi="max").reset_index())
cond = cond[cond.n >= 2].copy()
cond["mean_call"] = cond["mean"] < CUT
cond["spread"] = cond.hi - cond.lo
qc = t[["sample", "timepoint", "qc_fail"]].drop_duplicates()
cond = cond.merge(qc, on=["sample", "timepoint"])

j = t.merge(cond[["sample", "timepoint", "drug", "mean_call"]],
            on=["sample", "timepoint", "drug"])
j["agrees"] = j.call == j.mean_call
agree_all = j.agrees.mean()
agree_ok = j[~j.qc_fail].agrees.mean()
per = j.groupby(["sample", "timepoint", "drug"]).agrees.mean().rename("frac").reset_index()
cond = cond.merge(per, on=["sample", "timepoint", "drug"])
cond["unanimous"] = cond.frac == 1

print(f"conditions with >1 slice : {len(cond)}   slices: {len(j)}")
print(f"single slice reproduces its condition's mean call: {100*agree_all:.0f}% of slices "
      f"({100*agree_ok:.0f}% in QC-passing experiments)")
print(f"conditions where every slice agrees with the mean call: "
      f"{int(cond.unanimous.sum())} of {len(cond)}")
near = cond[(cond["mean"] - CUT).abs() <= 25]
far = cond[(cond["mean"] - CUT).abs() > 25]
print(f"  within 25 pp of the boundary : {int(near.unanimous.sum())}/{len(near)} unanimous")
print(f"  further than 25 pp away      : {int(far.unanimous.sum())}/{len(far)} unanimous")

fig, (axA, axB) = plt.subplots(1, 2, figsize=(183 / MM_PER_INCH, 78 / MM_PER_INCH),
                               gridspec_kw={"width_ratios": [1.25, 1.0]})
fig.subplots_adjust(left=0.075, right=0.985, top=0.855, bottom=0.275, wspace=0.26)

# ---------------------------------------------------------------- A: the slices themselves
srt = cond.sort_values("mean").reset_index(drop=True)
for i, r in srt.iterrows():
    v = t[(t["sample"] == r["sample"]) & (t.timepoint == r.timepoint)
          & (t.drug == r.drug)].pct_of_control
    c = "0.78" if r.qc_fail else (COLS[0] if not r.unanimous else COLS[1])
    axA.plot([r.lo, r.hi], [i, i], "-", color=c, lw=0.6, zorder=2, alpha=0.9)
    axA.plot(v, [i] * len(v), "o", ms=1.9, color=c, zorder=3, alpha=0.9)
    axA.plot(r["mean"], i, "|", ms=4.6, color="0.15", markeredgewidth=0.9, zorder=4)
axA.axvline(CUT, color="0.2", lw=1.0, ls=":", zorder=1)
axA.set_yticks([]); axA.set_ylim(-1.5, len(srt) + 0.5)
axA.set_xlim(-10, 210)
axA.set_xlabel("Ex vivo viability of one slice (% of DMSO control)", fontsize=6.4)
axA.set_ylabel(f"One row per condition (n = {len(srt)}), ordered by mean", fontsize=6.4)
axA.set_title("Every slice, and the mean the analysis uses", fontsize=7, pad=5)
axA.text(CUT - 4, len(srt) - 0.5,
         "left of this line = sensitive\nright = resistant", fontsize=5.0, ha="right",
         va="top", color="0.35", linespacing=1.5)
h = [plt.Line2D([], [], color=COLS[1], marker="o", ms=3, lw=0.8,
                label="every slice agrees with the mean call"),
     plt.Line2D([], [], color=COLS[0], marker="o", ms=3, lw=0.8,
                label="at least one slice disagrees"),
     plt.Line2D([], [], color="0.78", marker="o", ms=3, lw=0.8,
                label="experiment fails the 20% QC gate"),
     plt.Line2D([], [], color="0.15", marker="|", ms=5, lw=0, markeredgewidth=0.9,
                label="condition mean")]
axA.legend(handles=h, loc="lower right", frameon=False, fontsize=4.9, handlelength=1.4,
           handletextpad=0.5, borderaxespad=0.4)
finalize_axes(axA, tight=False)

# ---------------------------------------------------------------- B: where it goes wrong
for lab, sub, col, mfc in [("every slice agrees", cond[cond.unanimous], COLS[1], None),
                           ("at least one disagrees", cond[~cond.unanimous], COLS[0], "white")]:
    axB.plot(sub["mean"], sub.spread, "o", ms=4.2, color=col,
             markerfacecolor=(mfc or col), markeredgewidth=0.9, label=lab, zorder=3)
axB.axvline(CUT, color="0.2", lw=1.0, ls=":", zorder=1)
axB.set_xlim(-10, 210); axB.set_ylim(-4, cond.spread.max() * 1.10)
axB.set_xlabel("Condition mean viability (% of DMSO control)", fontsize=6.4)
axB.set_ylabel("Spread between slices\n(highest minus lowest, pp)", fontsize=6.4)
axB.set_title("Disagreement is concentrated near the call boundary", fontsize=7, pad=5)
axB.text(CUT + 3, cond.spread.max() * 1.04, "call boundary", fontsize=5.0, ha="left",
         va="top", color="0.35")
axB.legend(loc="upper right", frameon=False, fontsize=5.2, handlelength=1.2,
           handletextpad=0.5, borderaxespad=0.5)
finalize_axes(axB, tight=False)

for ax, L in [(axA, "A"), (axB, "B")]:
    ax.text(-0.10, 1.09, L, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")

CAP = (
    "A condition is one drug, on one biopsy, at one timepoint. Only conditions cultured as "
    "more than one slice are shown.\n"
    "A  Each row is one condition: dots are the individual slices, the tick is the mean the "
    f"analysis uses. A single slice reproduces its own condition's mean call in "
    f"{100*agree_all:.0f}% of slices "
    f"({100*agree_ok:.0f}% among experiments passing the quality-control gate).\n"
    f"B  The same conditions, plotted as mean against spread. Conditions whose slices disagree "
    f"sit close to the boundary: {int(near.unanimous.sum())} of {len(near)} are unanimous "
    f"within 25 pp of it, against {int(far.unanimous.sum())} of {len(far)} further away."
)
fig.text(0.075, 0.165, CAP, fontsize=4.9, va="top", linespacing=1.55)

save_figure(fig, "figS_replicate_agreement", outdir="figures")
