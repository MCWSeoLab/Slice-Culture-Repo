#!/usr/bin/env python3
"""Figure 3 — concordance between the ex vivo slice assay and the matched clinical course.

Rebuilt on `pairs_final.csv`, the adjudicated pair set (SAP v1.2 §13.3): chart abstraction
supplies PT-014, PT-015, PT-024 and PT-025; the workbook supplies the rest. Tissue-exhausted
specimens (PT-023) are drawn but held out of every statistic, exactly as they were before.

  A  ex vivo viability by binary clinical outcome (benefit vs progressive disease)
  B  ex vivo viability against real time on that regimen
  C  every matched course on one timeline, ordered by ex vivo viability
  D  2x2 at the 30%-inhibition threshold, with exact confidence intervals

Colour = agent class. Marker = biopsy site. Open face = the course began before the biopsy,
so the assay could not have informed it.
"""
import numpy as np, pandas as pd
from scipy import stats
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize, to_rgb
import itertools
import math
import re
import os
from matplotlib.gridspec import GridSpec
from cns_style import (apply_style, palette, finalize_axes, save_figure, MM_PER_INCH, plt)

apply_style("Nature")
SEED = 20260818
rng = np.random.default_rng(SEED)
B = 4000
CUT = 70.0        # % of DMSO control; >30% inhibition is the pre-specified call
# FIG3_OUTCOME=benefit (default) groups panel A as CR/PR/SD vs PD, the endpoint the SAP
# pre-specified. FIG3_OUTCOME=response groups it as CR/PR vs SD/PD.
#
# The two are NOT interchangeable statistically. Under `response` every patient's courses
# fall on the same side, so the comparison is 3 patients vs 3, not 5 pairs vs 5, and the
# naive pair-level Mann-Whitney (p = 0.032) is an artifact of treating 10 clustered courses
# as independent. `response` therefore reports an exact permutation of the outcome label
# across PATIENTS instead. See scripts/clustered_p.py for the full working.
# FIG3_LABELS=biopsy names panel C's rows the way Table 2 does - "Liver Bx 8", "PDAC Bx 5",
# "Peritoneum Bx 2" - instead of by study id. Each patient here contributed exactly one
# cultured biopsy, so nothing about the clustering is lost. FIG3_LABELS=patient (default)
# keeps the PT identifiers.
BXLAB = os.environ.get("FIG3_LABELS", "patient") == "biopsy"

OUTCOME = os.environ.get("FIG3_OUTCOME", "benefit")
RESPONSE = OUTCOME == "response"
if RESPONSE:
    BEN = "Objective response\n(CR/PR)"
    PDS = "No response\n(SD/PD)"
else:
    BEN = "Clinical benefit\n(CR/PR/SD)"
    PDS = "Progressive\ndisease"

# ---------------------------------------------------------------- data
p = pd.read_csv("pairs_final.csv")
site = (pd.read_csv("masterlist_linked_deid.csv")[["sample_label", "tissue_type",
                                                   "biopsy_site", "kras"]]
        .drop_duplicates("sample_label").set_index("sample_label"))
SITELAB = {("Metastatic", "Liver"): "Liver metastasis",
           ("Metastatic", "Peritoneum"): "Peritoneal metastasis",
           ("Primary", "Pancreas"): "Primary pancreas"}
p["site"] = p.sample_label.map(lambda s: SITELAB[(site.loc[s].tissue_type,
                                                  site.loc[s].biopsy_site)])
SITEMK = {"Liver metastasis": "o", "Peritoneal metastasis": "^", "Primary pancreas": "D"}

# `sample_label` in the analysis CSVs is the lab's name for the specimen; Table 2 numbers the
# same specimens per site. The join below is asserted against the KRAS allele in
# table2_source_deid.csv, so a wrong pairing fails the build instead of mislabelling a row.
TABLE2_SAMPLE = {"PDAC 5": "PDAC Bx 5",
                 "Umbilicus  1": "Peritoneum Bx 2",
                 "Liver 6 resection 2": "Liver Bx 6",
                 "Liver 7": "Liver Bx 7",
                 "Liver 8": "Liver Bx 8",
                 "Liver 9": "Liver Bx 9",
                 "Liver 10": "Liver Bx 10"}


def _allele(s):
    """KRAS allele as a bare token: 'G12D', 'WT', or '' when none is recorded."""
    s = (s or "").strip()
    if not s:
        return ""
    if "wild-type" in s.lower() or s.upper() == "WT":
        return "WT"
    m = re.search(r"KRAS\s+([A-Z]\d+[A-Z])", s)
    if m:
        return m.group(1)
    m = re.match(r"^([A-Z]\d+[A-Z])$", s)
    return m.group(1) if m else ""


def check_table2_join():
    """Fail loudly if a sample_label has been paired with the wrong Table 2 row."""
    t2 = pd.read_csv("table2_source_deid.csv", comment="#", keep_default_na=False)
    t2 = t2.set_index("sample")
    for lab, t2name in TABLE2_SAMPLE.items():
        assert t2name in t2.index, f"{t2name!r} is not a row in Table 2"
        want = _allele(site.loc[lab].kras if lab in site.index else "")
        got = _allele(t2.loc[t2name, "genomic"])
        if lab == "PDAC 5":
            # masterlist_linked_deid.csv still records G12R for this biopsy. Chart review
            # resolved the allele as KRAS G12D, TP53, GNAS, which is what Table 2 prints, so
            # the masterlist is not the arbiter for this one row.
            assert got == "G12D", f"PDAC Bx 5 allele moved: {got}"
            continue
        assert want == got, f"{lab} -> {t2name}: KRAS {want} vs Table 2 {got}"
    return t2


DRUGS = ["FOLFIRINOX", "Gem/nab-Pac", "MEK inhibitor"]
COL = dict(zip(DRUGS, palette("Nature", n=3)))
GREY = "0.72"
if BXLAB:
    check_table2_join()
    p["row_label"] = p.sample_label.map(TABLE2_SAMPLE)
    assert p.row_label.notna().all(), "a matched pair has no Table 2 name"
    # one biopsy per patient, so the relabelling is one-to-one and the clustering is intact
    assert p.groupby("study_id").row_label.nunique().eq(1).all()
    assert p.row_label.nunique() == p.study_id.nunique()
else:
    p["row_label"] = p.study_id

p["sens"] = p.mts_pct_of_control < CUT           # ex vivo sensitive
p["pre"] = p.day < 0                             # course started before the biopsy
p["dur"] = p.day_stop - p.day                    # real time on that regimen
p["censored"] = p.day_stop.isna()

p["responder"] = p.response.isin(["CR", "PR"])
p["pos"] = p.responder if RESPONSE else p.benefit.astype(bool)

a = p[~p.flag].copy()                            # the analysis set
NP, NPT = len(a), a.study_id.nunique()

TP = int((a.sens & a.pos).sum()); FP = int((a.sens & ~a.pos).sum())
FN = int((~a.sens & a.pos).sum()); TN = int((~a.sens & ~a.pos).sum())


def cp(k, n):
    """Clopper-Pearson exact interval; returns (point, lo, hi) as percentages."""
    if n == 0: return (np.nan, np.nan, np.nan)
    lo = 0.0 if k == 0 else stats.beta.ppf(0.025, k, n - k + 1)
    hi = 1.0 if k == n else stats.beta.ppf(0.975, k + 1, n - k)
    return (100 * k / n, 100 * lo, 100 * hi)


METRICS = [("Sensitivity", cp(TP, TP + FN)), ("Specificity", cp(TN, TN + FP)),
           ("PPV", cp(TP, TP + FP)), ("NPV", cp(TN, TN + FN)),
           ("Accuracy", cp(TP + TN, NP))]
FISHER = stats.fisher_exact([[TP, FP], [FN, TN]])[1]


def cluster_boot(df, x, y):
    """Spearman with a patient-level cluster bootstrap — pairs from one patient are not
    independent, so patients, not pairs, are the resampling unit."""
    ids = df.study_id.unique()
    idx = [np.flatnonzero((df.study_id == i).values) for i in ids]
    X, Y = df[x].to_numpy(float), df[y].to_numpy(float)
    # own generator: a shared stateful rng would make this depend on call order
    picks = np.random.default_rng(SEED).integers(0, len(ids), size=(B, len(ids)))
    out = []
    for b in range(B):
        sel = np.concatenate([idx[k] for k in picks[b]])
        xx, yy = X[sel], Y[sel]
        if len(np.unique(xx)) < 2 or len(np.unique(yy)) < 2: continue
        out.append(stats.spearmanr(xx, yy)[0])
    out = np.asarray(out)
    return stats.spearmanr(X, Y)[0], np.percentile(out, 2.5), np.percentile(out, 97.5)


def boot_band(df, x, y, xs):
    """95% band for a least-squares fit, resampling PATIENTS with replacement.

    Same resampling unit as `cluster_boot` and the rest of the SAP: courses from one patient
    are not independent, so a textbook OLS confidence band would be too narrow here. With only
    a handful of patients this band is wide, which is the honest picture, not a defect.
    """
    ids = df.study_id.unique()
    idx = [np.flatnonzero((df.study_id == i).values) for i in ids]
    X, Y = df[x].to_numpy(float), df[y].to_numpy(float)
    # own generator: a shared stateful rng would make this depend on call order
    picks = np.random.default_rng(SEED).integers(0, len(ids), size=(B, len(ids)))
    lines = []
    for b in range(B):
        sel = np.concatenate([idx[k] for k in picks[b]])
        xx, yy = X[sel], Y[sel]
        if len(np.unique(xx)) < 2: continue
        sl_b, ic_b = np.polyfit(xx, yy, 1)
        lines.append(ic_b + sl_b * xs)
    L = np.asarray(lines)
    return np.percentile(L, 2.5, axis=0), np.percentile(L, 97.5, axis=0), len(L)



def block_perm_p(df, x, y):
    """Exact p for the rank correlation, permuting each patient's y-block between patients.

    Methods describes "all 720 permutations of the six patients", but a patient contributing 2
    courses cannot receive a 1-course block, so only the size-preserving assignments are valid.
    Returns (p, n_valid, n_total) so the figure and the paper can state the real denominator.
    """
    ids = list(df.study_id.unique())
    blocks = {i: df.loc[df.study_id == i, y].to_numpy(float) for i in ids}
    X = df[x].to_numpy(float)
    obs = abs(stats.spearmanr(X, df[y].to_numpy(float))[0])
    cnt = valid = total = 0
    for perm in itertools.permutations(ids):
        total += 1
        if any(len(blocks[a]) != len(blocks[b]) for a, b in zip(ids, perm)):
            continue
        yy = np.concatenate([blocks[b] for b in perm])
        r = stats.spearmanr(X, yy)[0]
        valid += 1
        cnt += abs(r) >= obs - 1e-12
    return cnt / valid, valid, total

RHO, RLO, RHI = cluster_boot(a, "mts_pct_of_control", "rank")

# ---------------------------------------------------------------- canvas
# FIG3_LAYOUT=stacked (default) puts A and B on top with the timeline spanning the width;
# FIG3_LAYOUT=row puts all three side by side. The timeline carries 13 long y-tick labels and
# an 1100-day x range, so in `row` it gets a third of the width and the labels have to be
# shortened — compare the two before choosing.
LAYOUT = os.environ.get("FIG3_LAYOUT", "stacked")
ROW = LAYOUT == "row"

# FIG3_PANELC=heat replaces the "(139%)" suffix on panel C's tick labels with a colour-coded
# cell carrying the same number, in its own column to the left of the timeline, and rules a
# line between the last resistant and the first sensitive course. FIG3_PANELC=text (default)
# leaves panel C exactly as it was.
HEAT = os.environ.get("FIG3_PANELC", "text") == "heat"
HEAT_CMAP = plt.get_cmap("Purples")
# Deliberately NOT one of the three agent colours (red / cyan / teal) and not RdBu_r, so the
# viability column cannot be mistaken for a drug. Floor below 0 so a 1% cell is still tinted.
HEAT_NORM = Normalize(vmin=-45.0, vmax=152.0)

# FIG3_STATS=asis reproduces the annotation strings that are in the Illustrator-edited
# `fig3_concordance_response_row.pdf` currently pasted into the manuscript, so that a rebuild
# changes nothing but panel C. Every number in panel A is still derived from the data; the
# panel B "Exact p" is the one value that is carried through as a literal because nothing in
# this script computes it. FIG3_STATS=computed (default) prints the cluster-aware statistics
# this script derives, which for panel A is the exact permutation across patients, not the
# pair-level Mann-Whitney. See scripts/clustered_p.py.
STATS_ASIS = os.environ.get("FIG3_STATS", "computed") == "asis"

ASIS_PANEL_B_P = "Exact p=0.015"   # NOT recomputed here - carried from the edited PDF

if ROW:
    fig = plt.figure(figsize=(183 / MM_PER_INCH, 92 / MM_PER_INCH))
    # Panel C's tick labels are long and hang to the left of its axes box, so the gap
    # before it has to be wide enough to hold them or they run into panel B.
    gs = GridSpec(1, 3, figure=fig,
                  width_ratios=[1.0, 1.0, 1.62 if HEAT else 1.46],
                  wspace=(1.05 if BXLAB else 0.88) if HEAT else 1.02,
                  left=0.066, right=0.992, top=0.915, bottom=0.235)
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
    if HEAT:
        gsC = gs[0, 2].subgridspec(1, 2, width_ratios=[0.115, 1.0], wspace=0.055)
        axH = fig.add_subplot(gsC[0, 0]); axC = fig.add_subplot(gsC[0, 1])
    else:
        axH, axC = None, fig.add_subplot(gs[0, 2])
else:
    fig = plt.figure(figsize=(183 / MM_PER_INCH, 176 / MM_PER_INCH))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[0.88, 1.32], hspace=0.34, wspace=0.30,
                  left=0.115, right=0.975, top=0.950, bottom=0.130)
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
    if HEAT:
        gsC = gs[1, :].subgridspec(1, 2, width_ratios=[0.055, 1.0], wspace=0.03)
        axH = fig.add_subplot(gsC[0, 0]); axC = fig.add_subplot(gsC[0, 1])
    else:
        axH, axC = None, fig.add_subplot(gs[1, :])   # the timeline spans the full width


def draw(ax, r, x, y, ms=5.0, z=3):
    ax.plot(x, y, SITEMK[r.site], markersize=ms,
            markerfacecolor=(GREY if r.flag else ("white" if r.pre else COL[r.drug])),
            markeredgecolor=("0.45" if r.flag else COL[r.drug]),
            markeredgewidth=(0.6 if r.flag else 1.0), zorder=z)


# ---- A: viability by binary outcome
for k, (lab, mask) in enumerate([(BEN, p.pos), (PDS, ~p.pos)]):
    d = p[mask]
    off = np.linspace(-0.20, 0.20, len(d)) if len(d) > 1 else [0.0]
    for o, (_, r) in zip(off, d.sort_values("mts_pct_of_control").iterrows()):
        draw(axA, r, k + o, r.mts_pct_of_control)
    kept = d[~d.flag].mts_pct_of_control
    if len(kept):
        axA.plot([k - 0.30, k + 0.30], [kept.median()] * 2, "-", color="0.25", lw=1.4, zorder=2)
axA.axhline(100, color="0.35", lw=0.5, ls="--", zorder=1)
axA.axhline(CUT, color="0.35", lw=0.5, ls=":", zorder=1)
axA.text(1.42, 101, "DMSO control", fontsize=5.2, ha="right", va="bottom", color="0.35")
axA.text(1.42, CUT + 1, "30% inhibition", fontsize=5.2, ha="right", va="bottom", color="0.35")
# The response labels are wider than the benefit ones and collide in the one-row layout.
# Guarded so the already-generated `benefit` files are not disturbed.
axA.set_xticks([0, 1])
# "Objective response (CR/PR)" is too wide for a third-page panel even at 5.4 pt, so in the
# one-row layout it wraps onto three lines instead of shrinking further.
SHORTTICKS = os.environ.get("FIG3_TICKS", "auto") == "short"
if STATS_ASIS or SHORTTICKS:
    TICKS = ["Responder\n(CR/PR)", "Non-responder\n(SD/PD)"]
elif RESPONSE and ROW:
    TICKS = ["Objective\nresponse\n(CR/PR)", "No\nresponse\n(SD/PD)"]
else:
    TICKS = [BEN, PDS]
axA.set_xticklabels(TICKS, fontsize=(5.5 if (BXLAB and ROW) else 6.0) if (STATS_ASIS or SHORTTICKS)
                    else (5.6 if (RESPONSE and ROW) else 6.0))
# extra headroom in the one-row layout, where the annotation block is four lines deep in a
# narrower panel and would otherwise sit on the highest points
axA.set_xlim(-0.55, 1.55)
axA.set_ylim(0, 178 if ROW else 152)
axA.set_ylabel("Ex vivo viability at 48 h\n(% of DMSO control)", fontsize=6.4)
mb = a[a.pos].mts_pct_of_control; mp = a[~a.pos].mts_pct_of_control


def cluster_perm_p():
    """Exact p from permuting the outcome label across PATIENTS, not courses.

    Only defined when the outcome is constant within every patient, which it is under the
    `response` grouping and is not under `benefit` (PT-013 and PT-014 each contribute one SD
    and one PD, so they sit in both groups at once).
    """
    byp = a.groupby("study_id").pos.agg(nuniq="nunique", first="first")
    if (byp.nuniq > 1).any():
        return None, None, None
    ids = byp.index.to_list()
    k = int(byp["first"].sum())
    n1, n2 = len(mb), len(mp)
    obs = abs(stats.mannwhitneyu(mb, mp, alternative="two-sided")[0] - n1 * n2 / 2)
    cnt = tot = 0
    for sel in itertools.combinations(ids, k):
        m = a.study_id.isin(sel)
        pp, nn = a[m].mts_pct_of_control, a[~m].mts_pct_of_control
        if not len(pp) or not len(nn):
            continue
        u = stats.mannwhitneyu(pp, nn, alternative="two-sided")[0]
        tot += 1
        cnt += abs(u - len(pp) * len(nn) / 2) >= obs - 1e-9
    return cnt / tot, k, len(ids) - k


PERM_P, NPOS_PT, NNEG_PT = cluster_perm_p()
if STATS_ASIS:
    # The pair-level Mann-Whitney. It treats 10 clustered courses as 10 independent
    # observations and is an artifact; kept only so a rebuild matches the manuscript.
    _U = stats.mannwhitneyu(mb, mp, alternative="two-sided")[1]
    stat_lines = (f"n = {len(mb)} vs {len(mp)} courses, from {NPOS_PT} vs {NNEG_PT} patients\n"
                  f"median {mb.median():.0f}% vs {mp.median():.0f}%\n"
                  f"Mann-Whitney p={_U:.3f}")
elif PERM_P is not None:
    floor = 2 / math.comb(NPOS_PT + NNEG_PT, NPOS_PT)
    # The attainable-floor line is not drawn on the panel. `floor` is still computed and
    # printed to stdout, since it is the answer to why no smaller p is reportable here.
    stat_lines = (f"n = {len(mb)} vs {len(mp)} courses, from {NPOS_PT} vs {NNEG_PT} patients\n"
                  f"median {mb.median():.0f}% vs {mp.median():.0f}%\n"
                  f"exact p = {PERM_P:.2f}, outcome permuted across patients")
    print(f"  panel A: exact p={PERM_P:.4f}  smallest attainable at "
          f"{NPOS_PT} vs {NNEG_PT} patients = {floor:.4f} (not printed on the panel)")
else:
    U = stats.mannwhitneyu(mb, mp, alternative="two-sided")[1]
    stat_lines = (f"n = {len(mb)} vs {len(mp)} matched pairs\n"
                  f"median {mb.median():.0f}% vs {mp.median():.0f}%\n"
                  f"Mann–Whitney p = {U:.3f}\n"
                  "tissue-exhausted specimens shown, not counted")
axA.text(0.02, 0.985, stat_lines,
         transform=axA.transAxes, fontsize=5.0 if ROW else 5.2, va="top", linespacing=1.45)
finalize_axes(axA, tight=False)

# ---- B: viability vs real time on treatment
d = p[p.dur.notna()]
for _, r in d.iterrows():
    draw(axB, r, r.mts_pct_of_control, r.dur)
kept = d[~d.flag]
if len(kept) > 2:
    sl, ic, rv, pv, se = stats.linregress(kept.mts_pct_of_control, kept.dur)
    # Draw only over the observed x range. The line used to run 0-150, extrapolating past the
    # data; a confidence band over an extrapolated region would be worse than no band.
    xs = np.linspace(kept.mts_pct_of_control.min(), kept.mts_pct_of_control.max(), 80)
    blo, bhi, nfit = boot_band(kept, "mts_pct_of_control", "dur", xs)
    axB.fill_between(xs, blo, bhi, color="0.30", alpha=0.13, linewidth=0, zorder=1)
    axB.plot(xs, ic + sl * xs, "-", color="0.30", lw=0.9, zorder=2)
    rho2, lo2, hi2 = cluster_boot(kept, "mts_pct_of_control", "dur")
    if STATS_ASIS:
        _tail = ASIS_PANEL_B_P
    else:
        _bp, _nv, _nt = block_perm_p(kept, "mts_pct_of_control", "dur")
        # the {_nv}/{_nt} denominator belongs in Methods; spelled out here it overruns the
        # panel in the one-row layout and runs into panel C.
        _tail = (f"exact p = {_bp:.2f}, durations permuted across patients\n"
                 f"{len(kept)} courses with a recorded stop date, "
                 f"{kept.study_id.nunique()} patients")
        print(f"  panel B: rho={rho2:.4f} CI {lo2:.4f} to {hi2:.4f}  "
              f"exact p={_bp:.4f} ({_nv}/{_nt} valid)  n={len(kept)} courses "
              f"{kept.study_id.nunique()} patients")
    axB.text(0.02, 0.985,
             f"slope {sl:.2f} d per % viability,  R² = {rv**2:.2f}\n"
             f"Spearman ρ = {rho2:.2f} (95% CI {lo2:.2f} to {hi2:.2f})\n"
             + _tail,
             transform=axB.transAxes, ha="left", va="top",
             fontsize=5.0 if ROW else 5.2, linespacing=1.45)
axB.axvline(CUT, color="0.35", lw=0.5, ls=":", zorder=1)
axB.axvline(100, color="0.35", lw=0.5, ls="--", zorder=1)
# headroom for the four-line annotation block, which sits top-left
axB.set_xlim(-6, 152); axB.set_ylim(0, float(d.dur.max()) * 1.34)
axB.set_xlabel("Ex vivo viability at 48 h (% of DMSO control)", fontsize=6.4)
axB.set_ylabel("Time on that regimen (days)", fontsize=6.4)
finalize_axes(axB, tight=False)

# ---- C: every matched course on one timeline
# All 13 courses, so the row count matches A and B. The three tissue-exhausted PT-023 rows
# stay grey and out of the statistics; one of them ran five years before the biopsy and is
# shown as an off-scale arrow rather than being allowed to stretch the axis.
d = p.sort_values("mts_pct_of_control", ascending=False).reset_index(drop=True)
XL, XR = -680.0, 420.0
axC.axvspan(XL, 0, color="0.93", zorder=0)        # everything left of the line is retrospective
for i, r in d.iterrows():
    c = GREY if r.flag else COL[r.drug]
    end = r.day_stop if pd.notna(r.day_stop) else r.day + 22
    if end < XL:                                   # entirely off-scale to the left
        axC.annotate("", xy=(XL + 4, i), xytext=(XL + 40, i), zorder=3,
                     arrowprops=dict(arrowstyle="-|>", color=c, lw=1.1, shrinkA=0, shrinkB=0))
        axC.text(XL + 48, i, f"{r.day:.0f} to {end:.0f} d", ha="left", va="center",
                 fontsize=4.8, color="0.45")
        continue
    axC.barh(i, end - max(r.day, XL), left=max(r.day, XL), height=0.58, color=c,
             edgecolor="black", linewidth=0.3, zorder=2)
    if pd.isna(r.day_stop):
        axC.annotate("", xy=(end + 28, i), xytext=(end, i), zorder=3,
                     arrowprops=dict(arrowstyle="-|>", color=c, lw=0.8, shrinkA=0, shrinkB=0))
    else:
        axC.plot([end, end], [i - 0.30, i + 0.30], "-", color="black", lw=0.8, zorder=3)
    axC.text(max(r.day, XL) - 12, i, r.response, ha="right", va="center", fontsize=5.0,
             fontweight="bold", color="0.25")
axC.set_xlim(XL, XR)
axC.axvline(0, color="black", lw=0.8, zorder=4)
axC.set_yticks(range(len(d)))
# In the one-row layout the timeline gets a third of the page, so the viability value is
# dropped from the tick label and the type comes down to keep the labels from eating the axes.
ABBR = {"Gem/nab-Pac": "Gem/nab-P", "MEK inhibitor": "MEK inh", "FOLFIRINOX": "FOLFIRINOX"}
# Abbreviated in the row layout so the type can stay at 5 pt; 4 pt is below most journals'
# minimum. With HEAT the viability moves out of the label and into its own cell.
if ROW:
    TICKLAB = [f"{r.row_label} {ABBR.get(r.drug, r.drug)}"
               + ("" if HEAT else f" ({r.mts_pct_of_control:.0f}%)") for _, r in d.iterrows()]
    TICKFS = 5.0
else:
    TICKLAB = [f"{r.row_label} · {r.drug}"
               + ("" if HEAT else f" ({r.mts_pct_of_control:.0f}%)") for _, r in d.iterrows()]
    TICKFS = 5.2
if HEAT:
    axC.set_yticklabels([])
    axC.tick_params(axis="y", length=0)
else:
    axC.set_yticklabels(TICKLAB, fontsize=TICKFS)
axC.set_ylim(-0.7, len(d) - 0.3)
axC.set_xlabel("Days relative to biopsy", fontsize=6.4)
axC.set_title("Matched drug\u2013patient courses" if STATS_ASIS else
              (f"Matched drug\u2013patient courses\n({NP} pairs, {NPT} patients; "
               f"{len(d) - NP} tissue-exhausted, excluded)" if ROW else
               f"Matched drug\u2013patient courses ({NP} pairs, {NPT} patients; "
               f"{len(d) - NP} tissue-exhausted, excluded)"),
              fontsize=6.0 if ROW else 7, pad=4)
if ROW:
    axC.text(XL + 14, -0.64, "Before biopsy", fontsize=4.4, va="bottom", color="0.42")
    axC.text(16, -0.64, "After", fontsize=4.4, va="bottom", color="0.42")
else:
    axC.text(XL + 14, -0.64, "Course began before the biopsy", fontsize=5.2, va="bottom",
             color="0.42")
    axC.text(16, -0.64, "Began after", fontsize=5.2, va="bottom", color="0.42")
axC.invert_yaxis()
finalize_axes(axC, tight=False)

# The rule sits between the last course above the threshold and the first below it. Derived
# from the data, never from a row number: `d` is sorted by viability, so the first sensitive
# row is the boundary. It lands above PT-024 FOLFIRINOX (61%) for the current pair set.
_below = np.flatnonzero((d.mts_pct_of_control < CUT).to_numpy())
CUT_ROW = float(_below[0]) - 0.5 if len(_below) else None

if HEAT:
    axH.set_xlim(0, 1)
    axH.set_ylim(*axC.get_ylim())          # already inverted by axC
    axH.set_xticks([])
    axH.set_yticks(range(len(d)))
    axH.set_yticklabels(TICKLAB, fontsize=TICKFS)
    axH.tick_params(axis="y", length=0, pad=2)
    for sp in axH.spines.values():
        sp.set_visible(False)
    for i, r in d.iterrows():
        f = HEAT_CMAP(HEAT_NORM(r.mts_pct_of_control))
        axH.add_patch(Rectangle((0.0, i - 0.40), 1.0, 0.80, facecolor=f,
                                edgecolor="0.55", linewidth=0.3, zorder=2))
        # relative luminance of the fill decides black or white type, so every cell is legible
        lum = np.dot(to_rgb(f), (0.2126, 0.7152, 0.0722))
        axH.text(0.5, i, f"{r.mts_pct_of_control:.0f}", ha="center", va="center",
                 fontsize=4.8, color=("white" if lum < 0.55 else "0.12"), zorder=3)
    # tissue-exhausted specimens keep their measured value but are named in grey, matching
    # the grey bars and the legend entry
    for lab, (_, r) in zip(axH.get_yticklabels(), d.iterrows()):
        if r.flag:
            lab.set_color("0.45")
    # Captioned below the column, not above it: the space over panel C is already taken by
    # the panel title and the "Before biopsy" note, and a header wide enough to read would
    # have to run back across panel B.
    axH.text(0.42, -0.030, "Ex vivo\nviability at 48 h\n(% of DMSO control)",
             transform=axH.transAxes, ha="center", va="top",
             fontsize=4.6, color="0.25", linespacing=1.35)
    if CUT_ROW is not None:
        # runs past the right edge of the column so it meets panel C's rule across the gap
        axH.plot([-0.02, 1.32], [CUT_ROW] * 2, "-", color="0.15", lw=0.8,
                 clip_on=False, zorder=5)

if CUT_ROW is not None:
    axC.axhline(CUT_ROW, color="0.15", lw=0.8, zorder=5)
    axC.text(XR - 14, CUT_ROW + 0.14,
             "\u226530% ex vivo inhibition below",
             ha="right", va="top", fontsize=4.6, color="0.25", zorder=5)

# Panel D (the 2x2 and its exact intervals) was removed on request; the same numbers
# are reported in the Results text and the Statistical Analysis Plan.

# ---------------------------------------------------------------- legend
h1 = [Line2D([], [], marker="o", ls="none", ms=5, mfc=COL[d], mec=COL[d], label=d) for d in DRUGS]
h2 = [Line2D([], [], marker=SITEMK[s], ls="none", ms=5, mfc="0.35", mec="black", mew=0.4,
             label=f"{s} (n = {p[p.site == s].study_id.nunique()})") for s in SITEMK]
h3 = [Line2D([], [], marker="o", ls="none", ms=5, mfc="white", mec="0.25", mew=1.0,
             label="Open symbol (A, B) = course began before the biopsy"),
      Line2D([], [], marker="o", ls="none", ms=5, mfc=GREY, mec="0.45", mew=0.6,
             label="Tissue-exhausted specimen (not in statistics)"),
      Line2D([], [], marker="|", ls="none", ms=6, mec="black", label="Regimen stopped"),
      Line2D([], [], marker=">", ls="none", ms=5, mfc="0.35", mec="0.35",
             label="Still on regimen / no stop date recorded")]
fig.legend(handles=h1 + h2 + h3, loc="lower center", bbox_to_anchor=(0.5, 0.004),
           frameon=False, fontsize=5.0 if ROW else 5.6, ncol=4 if ROW else 3,
           handletextpad=0.5, columnspacing=1.2 if ROW else 1.4)

_axCL = axH if HEAT else axC
letters_pos = ([(axA, "A", -0.235), (axB, "B", -0.215),
                (_axCL, "C", -2.00 if HEAT else -0.400)] if ROW else
               [(axA, "A", -0.155), (axB, "B", -0.135),
                (_axCL, "C", -0.30 if HEAT else -0.075)])
for ax, L, dx in letters_pos:
    ax.text(dx, 1.045, L, transform=ax.transAxes, fontsize=9, fontweight="bold",
            va="bottom", ha="left")
name = ("fig3_concordance" + ("_response" if RESPONSE else "") + ("_row" if ROW else "")
        + ("_heat" if HEAT else "") + ("_bx" if BXLAB else "")
        + ("_asis" if STATS_ASIS else ""))
save_figure(fig, name, outdir="figures")
print(f"  outcome={OUTCOME} layout={LAYOUT}  {NP} pairs / {NPT} patients   "
      f"TP {TP} FP {FP} FN {FN} TN {TN}")
print(f"  rho {RHO:.3f}  CI {RLO:.3f} to {RHI:.3f}   Fisher p {FISHER:.4f}")
for n, (v, lo, hi) in METRICS:
    print(f"  {n:<12s} {v:5.1f}%  ({lo:.0f}-{hi:.0f})")
