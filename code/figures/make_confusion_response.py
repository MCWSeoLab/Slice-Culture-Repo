#!/usr/bin/env python3
"""Confusion matrix for the CR/PR vs SD/PD classification.

Ex vivo call: sensitive = >30% inhibition = viability < 70% of DMSO control at 48 h.
Clinical call: responder = CR or PR; non-responder = SD or PD.

Every count here is over drug-patient COURSES. Six patients contribute ten courses, and under
this grouping each patient's courses share one outcome, so the effective comparison is 3
patients vs 3. The intervals below do not know that and are therefore optimistic; the
patient-level test is in clustered_p.py.
"""
import numpy as np
import pandas as pd
from scipy import stats
from matplotlib.patches import Rectangle

from cns_style import apply_style, palette, finalize_axes, save_figure, MM_PER_INCH, plt

apply_style("Nature")
CUT = 70.0

p = pd.read_csv("pairs_final.csv")
a = p[~p.flag].copy()                       # tissue-exhausted specimens excluded
a["sens"] = a.mts_pct_of_control < CUT
a["pos"] = a.response.isin(["CR", "PR"])

TP = a[a.sens & a.pos]; FP = a[a.sens & ~a.pos]
FN = a[~a.sens & a.pos]; TN = a[~a.sens & ~a.pos]
n = {"TP": len(TP), "FP": len(FP), "FN": len(FN), "TN": len(TN)}


def cp(k, m):
    """Clopper-Pearson exact interval, as percentages."""
    if m == 0: return (np.nan, np.nan, np.nan)
    lo = 0.0 if k == 0 else stats.beta.ppf(0.025, k, m - k + 1)
    hi = 1.0 if k == m else stats.beta.ppf(0.975, k + 1, m - k)
    return (100 * k / m, 100 * lo, 100 * hi)


METRICS = [
    ("Sensitivity", n["TP"], n["TP"] + n["FN"]),
    ("Specificity", n["TN"], n["TN"] + n["FP"]),
    ("PPV",         n["TP"], n["TP"] + n["FP"]),
    ("NPV",         n["TN"], n["TN"] + n["FN"]),
    ("Accuracy",    n["TP"] + n["TN"], len(a)),
]
FISHER = stats.fisher_exact([[n["TP"], n["FP"]], [n["FN"], n["TN"]]])[1]

print(f"courses {len(a)}  patients {a.study_id.nunique()}")
print(f"TP {n['TP']}  FP {n['FP']}  FN {n['FN']}  TN {n['TN']}")
for lab, k, m in METRICS:
    v, lo, hi = cp(k, m)
    print(f"  {lab:<12s} {k}/{m} = {v:5.0f}%   (95% CI {lo:.0f}-{hi:.0f})")
print(f"  Fisher exact p = {FISHER:.3f}")
print("\ncell membership:")
for tag, d in [("TP", TP), ("FP", FP), ("FN", FN), ("TN", TN)]:
    for _, r in d.iterrows():
        print(f"  {tag}  {r.study_id} {r.drug:<14s} {r.mts_pct_of_control:6.1f}%  {r.response}")

# ------------------------------------------------------------------ figure
fig = plt.figure(figsize=(120 / MM_PER_INCH, 62 / MM_PER_INCH))
ax = fig.add_axes([0.02, 0.02, 0.96, 0.96]); ax.axis("off")

x0, y0, w, h = 0.16, 0.46, 0.44, 0.44
COLS = ["Objective response\n(CR/PR)", "No response\n(SD/PD)"]
ROWS = ["Ex vivo sensitive\n(>30% inhibition)", "Ex vivo resistant\n(≤30% inhibition)"]

ax.text(x0 + w / 2, y0 + h + 0.115, "Clinical outcome", ha="center", va="bottom",
        fontsize=6.6, fontweight="bold")
for j, lab in enumerate(COLS):
    ax.text(x0 + w * (0.25 + 0.5 * j), y0 + h + 0.028, lab, ha="center", va="bottom",
            fontsize=5.8)
for i, lab in enumerate(ROWS):
    ax.text(x0 - 0.022, y0 + h * (0.75 - 0.5 * i), lab, ha="right", va="center", fontsize=5.8)

grid = [[("TP", n["TP"]), ("FP", n["FP"])], [("FN", n["FN"]), ("TN", n["TN"])]]
for i, row in enumerate(grid):
    for j, (tag, k) in enumerate(row):
        cx = x0 + w * (0.25 + 0.5 * j)
        cy = y0 + h * (0.75 - 0.5 * i)
        ax.add_patch(Rectangle((x0 + w * 0.5 * j, y0 + h * (0.5 - 0.5 * i)), w * 0.5, h * 0.5,
                               facecolor="none", edgecolor="0.55", lw=0.6))
        ax.text(cx, cy + 0.035, str(k), ha="center", va="center", fontsize=15)
        ax.text(cx, cy - 0.062, tag, ha="center", va="center", fontsize=5.4, color="0.45")

lines = [f"{lab:<12s} {k}/{m} = {cp(k, m)[0]:3.0f}%   (95% CI {cp(k, m)[1]:.0f}–{cp(k, m)[2]:.0f})"
         for lab, k, m in METRICS]
lines += ["", f"{'Fisher exact':<12s} p = {FISHER:.3f}"]
ax.text(0.655, y0 + h, "\n".join(lines), va="top", ha="left", fontsize=5.4,
        family="monospace", linespacing=1.62)

ax.text(0.0, 0.335,
        f"{len(a)} drug–patient courses from {a.study_id.nunique()} patients; three "
        "tissue-exhausted courses excluded.\n"
        "Counts and intervals are over courses. Under this grouping every patient's courses "
        "share one\noutcome, so the comparison is 3 patients versus 3; the intervals shown do "
        "not account for that.\n"
        "Ex vivo sensitive is <70% of same-plate DMSO control at 48 h.",
        va="top", ha="left", fontsize=5.0, linespacing=1.55)

save_figure(fig, "fig3_confusion_response", outdir="figures")
