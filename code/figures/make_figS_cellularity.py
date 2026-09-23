#!/usr/bin/env python3
"""Supplementary figure — tumour content of the cultured slices.

A reviewer's first question about a whole-slice viability readout is what fraction of the
slice is actually tumour. EpCAM staining on the DMSO control slices answers it directly:
in a liver metastasis the surrounding hepatocytes are EpCAM-low, so EpCAM-positive area is
a reasonable proxy for tumour content.

  A  EpCAM-positive fraction of each biopsy's DMSO control slices
  B  ex vivo viability against tumour content, for the matched pairs

This is descriptive. With six patients in the concordance set it cannot support a claim that
cellularity does or does not drive the assay result, and the figure does not make one.
"""
import numpy as np, pandas as pd
from scipy import stats
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from cns_style import apply_style, palette, finalize_axes, save_figure, MM_PER_INCH, plt

apply_style("Nature")

# lab IHC sheet name -> study id, where the link is unambiguous. The two date-redacted
# sheets already carry their study_id in ihc_tidy.csv; the rest match on the lab's own
# biopsy number, which the de-identified masterlist confirms.
LINK = {"Liver Bx.5": "PT-012", "Liver Bx.6": "PT-014", "Liver Bx.7": "PT-015",
        "Liver bx. [DATE]": "PT-024", "Liver Bx. [DATE]": "PT-025",
        "PDAC 5": "PT-011", "Umbilicus": "PT-013"}
SITE = {"Liver": "Liver metastasis", "PDAC": "Primary pancreas",
        "Umbilicus": "Peritoneal metastasis"}
SITECOL = dict(zip(["Liver metastasis", "Peritoneal metastasis", "Primary pancreas"],
                   palette("Nature", n=3)))


def site_of(lbl):
    if lbl.lower().startswith("liver"): return "Liver metastasis"
    if lbl.lower().startswith("pdac"):  return "Primary pancreas"
    return "Peritoneal metastasis"


ihc = pd.read_csv("ihc_tidy.csv")
ep = ihc[(ihc.marker == "EpCAM") & (ihc.treatment == "DMSO")].copy()
g = (ep.groupby("biopsy").value.agg(["mean", "std", "count"])
     .reset_index().sort_values("mean"))
g["site"] = g.biopsy.map(site_of)
g["sid"] = g.biopsy.map(LINK)

pairs = pd.read_csv("pairs_final.csv")
inpairs = set(pairs[~pairs.flag].study_id)

fig = plt.figure(figsize=(183 / MM_PER_INCH, 82 / MM_PER_INCH))
gs = GridSpec(1, 2, figure=fig, width_ratios=[1.35, 1.0], wspace=0.30,
              left=0.085, right=0.985, top=0.90, bottom=0.30)

# ---- A: tumour content per biopsy
axA = fig.add_subplot(gs[0])
for i, r in enumerate(g.itertuples()):
    axA.bar(i, r.mean, width=0.66, color=SITECOL[r.site], edgecolor="black",
            linewidth=0.4, zorder=2)
    if np.isfinite(r.std):
        axA.plot([i, i], [r.mean, r.mean + r.std], "-", color="0.25", lw=0.7, zorder=3)
        axA.plot([i - 0.13, i + 0.13], [r.mean + r.std] * 2, "-", color="0.25", lw=0.7, zorder=3)
    if r.sid in inpairs:
        axA.plot(i, 1.5, "^", ms=3.4, color="0.20", clip_on=False, zorder=4)
axA.set_xticks(range(len(g)))
axA.set_xticklabels([sid if isinstance(sid, str) else b.replace(" [DATE]", "")
                     for b, sid in zip(g.biopsy, g.sid)],
                    rotation=32, ha="right", fontsize=5.6)
axA.set_ylabel("EpCAM-positive area in DMSO control\n(% of slice)", fontsize=6.4)
axA.set_ylim(0, 85)
axA.text(0.02, 0.97,
         f"{len(g)} biopsies, {int(g['count'].sum())} DMSO control slices\n"
         f"median {g['mean'].median():.0f}%, range {g['mean'].min():.0f}–{g['mean'].max():.0f}%\n"
         "▲ = biopsy contributes to the concordance analysis",
         transform=axA.transAxes, va="top", fontsize=5.2, linespacing=1.45)
finalize_axes(axA, tight=False)
axA.text(-0.135, 1.06, "A", transform=axA.transAxes, fontsize=9, fontweight="bold", va="bottom")

# ---- B: viability against tumour content, matched pairs only
axB = fig.add_subplot(gs[1])
cell = g.dropna(subset=["sid"]).set_index("sid")["mean"]
q = pairs[~pairs.flag].copy()
q["cellularity"] = q.study_id.map(cell)
q = q.dropna(subset=["cellularity"])
for _, r in q.iterrows():
    axB.plot(r.cellularity, r.mts_pct_of_control, "o" if r.benefit else "s", ms=5,
             markerfacecolor=("white" if r.benefit else "0.35"),
             markeredgecolor="0.20", markeredgewidth=0.9, zorder=3)
axB.axhline(70, color="0.35", lw=0.5, ls=":", zorder=1)
axB.set_xlabel("Tumour content of the biopsy (% EpCAM)", fontsize=6.4)
axB.set_ylabel("Ex vivo viability at 48 h\n(% of DMSO control)", fontsize=6.4)
axB.set_xlim(0, 80); axB.set_ylim(0, 130)
if len(q) > 2:
    # same standard as Figure 3: patients are the resampling unit, and the naive p-value
    # from 10 non-independent pairs is not reported
    rng = np.random.default_rng(20260818)
    ids = q.study_id.unique()
    idx = [np.flatnonzero((q.study_id == i).values) for i in ids]
    X, Y = q.cellularity.to_numpy(float), q.mts_pct_of_control.to_numpy(float)
    boot = []
    for b in range(4000):
        sel = np.concatenate([idx[k] for k in rng.integers(0, len(ids), len(ids))])
        if len(np.unique(X[sel])) < 2 or len(np.unique(Y[sel])) < 2: continue
        boot.append(stats.spearmanr(X[sel], Y[sel])[0])
    rho = stats.spearmanr(X, Y)[0]
    lo, hi = np.percentile(boot, [2.5, 97.5])
    axB.text(0.03, 0.03,
             f"Spearman ρ = {rho:+.2f}\n"
             f"cluster-bootstrap 95% CI {lo:+.2f} to {hi:+.2f}\n"
             f"{len(q)} pairs, {q.study_id.nunique()} patients",
             transform=axB.transAxes, va="bottom", fontsize=5.2, linespacing=1.45)
    print(f"\ncellularity vs viability: rho {rho:.3f}  CI {lo:.3f} to {hi:.3f}")
finalize_axes(axB, tight=False)
axB.text(-0.20, 1.06, "B", transform=axB.transAxes, fontsize=9, fontweight="bold", va="bottom")

h = [Line2D([], [], marker="s", ls="none", ms=5, mfc=SITECOL[s], mec="black", mew=0.4, label=s)
     for s in SITECOL]
h += [Line2D([], [], marker="o", ls="none", ms=5, mfc="white", mec="0.20", mew=0.9,
             label="Clinical benefit"),
      Line2D([], [], marker="s", ls="none", ms=5, mfc="0.35", mec="0.20", mew=0.9,
             label="Progressive disease")]
fig.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, 0.005), frameon=False,
           fontsize=5.8, ncol=5, handletextpad=0.5, columnspacing=1.4)

save_figure(fig, "figS_tumour_content", outdir="figures")

out = g[["biopsy", "site", "sid", "mean", "std", "count"]].copy()
out.columns = ["biopsy", "site", "study_id", "epcam_pct_mean", "epcam_pct_sd", "n_slices"]
out.round(2).to_csv("tumour_content_deid.csv", index=False)
print(out.round(1).to_string(index=False))
print(f"\nmedian {g['mean'].median():.1f}%  range {g['mean'].min():.1f}-{g['mean'].max():.1f}%")
print(f"below 10% EpCAM: {int((g['mean'] < 10).sum())} of {len(g)} biopsies")
