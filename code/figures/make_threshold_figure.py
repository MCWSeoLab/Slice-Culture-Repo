#!/usr/bin/env python3
"""Where the >30% inhibition threshold comes from.

Rebuilt 20 Aug 2026 against the WELL-LEVEL lab MTS workbooks (`scripts/extract_mts_wells.py`
-> `mts_dmso_cv_wells_deid.csv`, `mts_treated_pct_deid.csv`). The earlier version used the
Prism extract, which held only a fraction of the control wells and made both derivations look
unsupportable. With the real wells, both reproduce.

  A  DMSO replicate CV per experiment — the measured control-well variability
  B  the two independent noise estimates, and where each puts the cutoff
  C  how the classification moves as the cutoff is varied
"""
import numpy as np
import pandas as pd

from cns_style import apply_style, palette, finalize_axes, save_figure, MM_PER_INCH, plt

apply_style("Nature")
COLS = palette("Nature", n=3)
QC = 0.20

cv = pd.read_csv("mts_dmso_cv_wells_deid.csv")
tp = pd.read_csv("mts_treated_pct_deid.csv")

rep = cv[cv.n_control_wells >= 2].sort_values("dmso_cv").reset_index(drop=True)
pas = rep[rep.qc == "PASS"]
MED = float(pas.dmso_cv.median())
SD_CV = 100 * MED * np.sqrt(2)
CUT_CV = 2 * SD_CV

# half-normal fit to treated wells above 100% of control, QC-passing experiments only
d = (tp.loc[(~tp.qc_fail) & (tp.pct_of_control > 100), "pct_of_control"] - 100).to_numpy(float)
rng = np.random.default_rng(20260818)
hn = lambda x: float(np.sqrt(np.mean(np.asarray(x, float) ** 2)))
SD_HN = hn(d)
BS = [hn(d[rng.integers(0, len(d), len(d))]) for _ in range(4000)]
HN_LO, HN_HI = np.percentile(BS, 2.5), np.percentile(BS, 97.5)
CUT_HN = 2 * SD_HN

print(f"A  experiments with >=2 DMSO wells {len(rep)}, passing QC {len(pas)} "
      f"from {pas['sample'].nunique()} biopsies; max control wells {int(cv.n_control_wells.max())}")
print(f"B  CV route          median {MED*100:.2f}%  -> SD {SD_CV:.2f} pp -> {CUT_CV:.1f}%")
print(f"   half-normal route sigma {SD_HN:.2f} pp (CI {HN_LO:.1f}-{HN_HI:.1f}), "
      f"{len(d)} wells >100% -> {CUT_HN:.1f}%")

fig, axes = plt.subplots(1, 3, figsize=(183 / MM_PER_INCH, 66 / MM_PER_INCH))
fig.subplots_adjust(left=0.075, right=0.988, top=0.855, bottom=0.265, wspace=0.44)
axA, axB, axC = axes

# ---------------------------------------------------------------- A
for i, r in rep.iterrows():
    fail = str(r.qc).startswith("FAIL")
    axA.plot(r.dmso_cv * 100, i, "o", ms=3.4,
             markerfacecolor=("0.78" if fail else COLS[0]),
             markeredgecolor=("0.5" if fail else COLS[0]), markeredgewidth=0.5, zorder=3)
axA.axvline(MED * 100, color="0.25", lw=1.0, zorder=2)
axA.axvline(QC * 100, color="0.35", lw=0.6, ls=":", zorder=1)
axA.set_yticks([]); axA.set_ylim(-1.2, len(rep) + 0.4); axA.set_xlim(-3, 90)
axA.set_ylabel(f"Experiments (n = {len(rep)}), ranked", fontsize=6.4)
axA.set_xlabel("DMSO replicate CV (%)", fontsize=6.4)
axA.set_title("Control-well variability\nper experiment", fontsize=7, pad=5)
axA.text(MED * 100 + 2.5, len(rep) - 0.5, f"median of the {len(pas)}\npassing {MED*100:.1f}%",
         fontsize=5.0, va="top", ha="left", color="0.25")
axA.text(QC * 100 + 1.5, 0.6, "QC gate 20%", fontsize=5.0, ha="left", va="bottom",
         color="0.45", rotation=90)
axA.text(0.985, 0.035, f"{len(rep) - len(pas)} of {len(rep)} fail the gate",
         transform=axA.transAxes, fontsize=4.8, ha="right", va="bottom", color="0.45")
finalize_axes(axA, tight=False)

# ---------------------------------------------------------------- B
xs = np.linspace(30, 170, 400)
dens = np.exp(-0.5 * ((xs - 100) / SD_CV) ** 2) / (SD_CV * np.sqrt(2 * np.pi))
axB.plot(xs, dens, "-", color="0.30", lw=0.9, zorder=3)
band = (xs >= 100 - 2 * SD_CV) & (xs <= 100 + 2 * SD_CV)
axB.fill_between(xs[band], 0, dens[band], color="0.30", alpha=0.12, linewidth=0, zorder=1)
axB.axvline(100, color="0.35", lw=0.6, ls="--", zorder=2)
axB.axvline(100 - CUT_CV, color=COLS[0], lw=1.1, zorder=4)
axB.axvline(100 - CUT_HN, color=COLS[2], lw=1.1, ls=(0, (4, 2)), zorder=4)
axB.set_xlim(30, 170); axB.set_ylim(0, dens.max() * 1.62)
axB.set_xlabel("Viability (% of DMSO control)", fontsize=6.4)
axB.set_ylabel("Expected density under\nno drug effect", fontsize=6.4)
axB.set_title("Two independent estimates\nof the same noise floor", fontsize=7, pad=5)
axB.text(0.03, 0.985,
         f"control-well CV route\n  SD {SD_CV:.1f} pp → cutoff {CUT_CV:.0f}%\n\n"
         f"wells reading >100% route\n  SD {SD_HN:.1f} pp ({HN_LO:.0f}–{HN_HI:.0f}) → "
         f"cutoff {CUT_HN:.0f}%\n\nthreshold set at 30%",
         transform=axB.transAxes, fontsize=4.9, va="top", ha="left", linespacing=1.5)
finalize_axes(axB, tight=False)

# ---------------------------------------------------------------- C
p = pd.read_csv("pairs_final.csv")
a = p[~p.flag].copy()
a["pos"] = a.response.isin(["CR", "PR"])
cuts = np.arange(10, 61, 1)
g = pd.DataFrame([dict(cut=c,
                       TP=int(((a.mts_pct_of_control < 100 - c) & a.pos).sum()),
                       FP=int(((a.mts_pct_of_control < 100 - c) & ~a.pos).sum()),
                       FN=int((~(a.mts_pct_of_control < 100 - c) & a.pos).sum()),
                       TN=int((~(a.mts_pct_of_control < 100 - c) & ~a.pos).sum()))
                  for c in cuts])
g["correct"] = g.TP + g.TN
obs = np.sort(a.mts_pct_of_control.values)
edges = np.unique(np.round(np.concatenate(([0.0], 100 - obs[::-1], [100.0])), 4))
edges = edges[(edges >= 0) & (edges <= 100)]
FLAT_LO = float(edges[edges <= 30].max()); FLAT_HI = float(edges[edges > 30].min())
axC.axvspan(FLAT_LO, FLAT_HI, color="0.30", alpha=0.11, linewidth=0, zorder=1)
for lab, col, c in [("Correctly classified", "correct", COLS[0]),
                    ("False positives", "FP", COLS[1]), ("False negatives", "FN", COLS[2])]:
    axC.step(g.cut, g[col], where="post", lw=1.1, color=c, label=lab, zorder=3)
axC.axvline(30, color="0.25", lw=1.0, zorder=2)
axC.set_xlim(10, 60); axC.set_ylim(-0.4, len(a) + 0.6)
axC.set_xlabel("Cutoff (% inhibition)", fontsize=6.4)
axC.set_ylabel(f"Courses (of {len(a)})", fontsize=6.4)
axC.set_title("Classification against cutoff\n(CR/PR vs SD/PD)", fontsize=7, pad=5)
axC.text(30.6, len(a) + 0.35, "30%", fontsize=5.2, ha="left", va="top", color="0.25")
axC.text((FLAT_LO + FLAT_HI) / 2, 0.35, f"no call changes\n{FLAT_LO:.0f}–{FLAT_HI:.0f}%",
         fontsize=4.8, ha="center", va="bottom", color="0.42", linespacing=1.5)
axC.legend(loc="center left", frameon=False, fontsize=5.2, handlelength=1.4,
           handletextpad=0.5, borderaxespad=0.4)
finalize_axes(axC, tight=False)

for ax, L, dx in [(axA, "A", -0.20), (axB, "B", -0.235), (axC, "C", -0.215)]:
    ax.text(dx, 1.10, L, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom")

fig.text(0.075, 0.135,
         f"A  Coefficient of variation among same-plate DMSO control wells, one point per "
         f"experiment and timepoint, ranked; grey points fail the 20% QC gate. Median of the "
         f"{len(pas)} passing experiments {MED*100:.1f}%.\n"
         f"B  A treated well normalised to a control well carries that error twice, giving "
         f"SD {SD_CV:.1f} pp. Independently, treated wells reading above 100% of control carry "
         f"no drug effect by construction; a half-normal\n"
         f"    fitted to that tail ({len(d)} wells) gives SD {SD_HN:.1f} pp. Two standard "
         f"deviations below control is {CUT_CV:.0f}% and {CUT_HN:.0f}% inhibition respectively.\n"
         f"C  Counts over the {len(a)} matched courses as the cutoff is varied; shading marks "
         f"the {FLAT_LO:.0f}–{FLAT_HI:.0f}% window within which no call changes.",
         fontsize=4.9, va="top", linespacing=1.55)

save_figure(fig, "figS_threshold_derivation", outdir="figures")
