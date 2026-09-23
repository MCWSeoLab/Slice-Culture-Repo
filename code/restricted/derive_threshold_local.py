#!/usr/bin/env python3
"""Regenerate BOTH derivations of the response threshold from well-level Masterlist data.

RUN THIS LOCALLY. It reads Masterlist.xlsx, which contains PHI and must not leave the machine.
Everything it writes is de-identified aggregate: counts, CVs, sigma, and the cutoffs. Paste the
printed report back into the session, or hand over the two CSVs it writes.

    python3 derive_threshold_local.py
    python3 derive_threshold_local.py --path "/some/other/Masterlist.xlsx"

Dependencies: pandas, numpy, openpyxl. Deliberately NO scipy — the analysis VM does not have
it. The half-normal MLE and the bootstrap are a few lines of numpy.

WHY THIS EXISTS
The manuscript says two independent routes converged on ~30% inhibition. Only one of them
(propagated DMSO CV) can currently be reproduced from the de-identified files, and it rests on
three CVs from two patients. The other (half-normal fit to treated wells reading >100% of
control) needs well-level MTS data, which lives only here. If "two independent routes" is going
to stay in the Methods, this script is what backs it up.

WHAT IT DOES
  1. Prints a schema report so column guesses can be checked before trusting any number.
  2. Derivation B — DMSO replicate CV per experiment, and a full census of how many
     experiments could and could not contribute one.
  3. Derivation A — half-normal fit to treated wells above 100% of control, with a bootstrap
     interval on sigma, computed both before and after the QC gate.
  4. Writes deid/threshold_cv_census_deid.csv and deid/threshold_halfnormal_deid.csv.
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

# --------------------------------------------------------------------------- config
SHEET = "Metadata"
QC_GATE = 0.20          # SAP 4.3: exclude any experiment whose DMSO replicate CV exceeds 20%
N_BOOT = 4000
SEED = 20260818
CONTROL_WORDS = ("dmso", "control", "untreated", "vehicle", "ctrl")

DATE_RE = re.compile(r"\b\d{1,2}[-/.]\d{1,2}[-/.]\d{2,4}\b|\b\d{4}[-/.]\d{1,2}[-/.]\d{1,2}\b")


def scrub(s):
    """Redact anything date-shaped from a label before it is written out."""
    return DATE_RE.sub("[DATE]", str(s)).strip()


def pick(cols, *patterns, required=True, what=""):
    """Find the one column matching any pattern. Refuse on ambiguity rather than guess."""
    hits = [c for c in cols if any(re.search(p, str(c), re.I) for p in patterns)]
    if len(hits) == 1:
        return hits[0]
    if not hits:
        if required:
            sys.exit(f"\nFAILED: no column found for {what}. Looked for {patterns}.\n"
                     f"Columns present: {list(cols)}\n"
                     f"Edit the pattern at the top of this script and re-run.")
        return None
    sys.exit(f"\nFAILED: {len(hits)} columns match {what}: {hits}\n"
             f"Narrow the pattern in this script so exactly one matches, then re-run.")


def halfnormal_sigma(d):
    """MLE for a half-normal folded at zero: sigma = sqrt(mean(d^2)). No scipy needed."""
    d = np.asarray(d, float)
    return float(np.sqrt(np.mean(d ** 2))) if len(d) else float("nan")


def boot_sigma(d, rng, n=N_BOOT):
    d = np.asarray(d, float)
    if len(d) < 2:
        return (float("nan"), float("nan"))
    out = [halfnormal_sigma(d[rng.integers(0, len(d), len(d))]) for _ in range(n)]
    return (float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5)))


# --------------------------------------------------------------------------- load
ap = argparse.ArgumentParser()
ap.add_argument("--path", default=None, help="path to Masterlist.xlsx")
ap.add_argument("--outdir", default="deid", help="where the de-identified CSVs go")
args = ap.parse_args()

path = args.path
if path is None:
    here = os.path.dirname(os.path.abspath(__file__))
    for cand in [os.path.join(here, "..", "Manuscript Drafts with Results", "Correlation analysis",
                              "Masterlist.xlsx"),
                 os.path.join(here, "..", "Masterlist.xlsx"),
                 os.path.join(here, "Masterlist.xlsx")]:
        if os.path.exists(cand):
            path = cand
            break
if path is None or not os.path.exists(path):
    sys.exit("Could not find Masterlist.xlsx. Pass it explicitly:\n"
             '  python3 derive_threshold_local.py --path "/full/path/Masterlist.xlsx"')

m = pd.read_excel(path, sheet_name=SHEET)
m.columns = [str(c).strip() for c in m.columns]

print("=" * 78)
print("0. SCHEMA — check these guesses before trusting anything below")
print("=" * 78)
print(f"file  : {os.path.basename(path)}")
print(f"sheet : {SHEET}   rows: {len(m)}   columns: {len(m.columns)}")
print("\ncolumns:")
for c in m.columns:
    nn = m[c].notna().sum()
    print(f"    {c:<28s} non-null {nn:>5d}   e.g. {str(m[c].dropna().iloc[0])[:34] if nn else '-'}")

C_SAMPLE = pick(m.columns, r"^sample", r"sample.*id", what="the sample identifier")
C_TREAT = pick(m.columns, r"^treat", r"drug", r"condition", what="the treatment")
C_MTS = pick(m.columns, r"mts", what="the MTS percent-of-control")
C_TIME = pick(m.columns, r"timepoint", r"time.*h", r"^hr", required=False, what="the timepoint")
C_RAW = pick(m.columns, r"absorb", r"^od\b", r"raw", r"490", required=False,
             what="raw absorbance")

print(f"\nusing:  sample={C_SAMPLE!r}  treatment={C_TREAT!r}  mts={C_MTS!r}  "
      f"timepoint={C_TIME!r}  raw={C_RAW!r}")
if C_RAW is None:
    print("  NOTE: no raw-absorbance column found. DMSO CVs will be computed from the\n"
          "        percent-of-control column instead, which is valid only if control wells\n"
          "        are normalized to the group mean rather than to themselves. Check the\n"
          "        DMSO values printed below: if they are all exactly 100, this is wrong and\n"
          "        the CV route needs the raw plate readings.")

m[C_MTS] = pd.to_numeric(m[C_MTS], errors="coerce")
if C_RAW:
    m[C_RAW] = pd.to_numeric(m[C_RAW], errors="coerce")
m["_is_control"] = m[C_TREAT].astype(str).str.lower().str.replace(r"[^a-z]", "", regex=True) \
                    .apply(lambda s: any(w in s for w in CONTROL_WORDS))
m["_expt"] = (m[C_SAMPLE].astype(str) + " | " +
              (m[C_TIME].astype(str) if C_TIME else "all"))

print(f"\ncontrol rows: {int(m._is_control.sum())}   treated rows: {int((~m._is_control).sum())}")
ctrl_vals = m.loc[m._is_control, C_MTS].dropna()
if len(ctrl_vals):
    print(f"control {C_MTS} values: min {ctrl_vals.min():.1f}  median "
          f"{ctrl_vals.median():.1f}  max {ctrl_vals.max():.1f}")

# --------------------------------------------------------------------- derivation B
print()
print("=" * 78)
print("B. DMSO REPLICATE CV PER EXPERIMENT")
print("=" * 78)
val = C_RAW if C_RAW else C_MTS
rows = []
for expt, g in m[m._is_control].groupby("_expt"):
    v = g[val].dropna()
    rows.append(dict(experiment=scrub(expt), n_control_wells=len(v),
                     mean=float(v.mean()) if len(v) else np.nan,
                     sd=float(v.std(ddof=1)) if len(v) > 1 else np.nan,
                     dmso_cv=float(v.std(ddof=1) / v.mean()) if len(v) > 1 and v.mean() else np.nan))
cv = pd.DataFrame(rows).sort_values("n_control_wells", ascending=False)
cv["qc"] = np.where(cv.dmso_cv.isna(), "UNKNOWN (single control well)",
                    np.where(cv.dmso_cv > QC_GATE, "FAIL (CV > 20%)", "PASS"))

n_total = len(cv)
n_multi = int((cv.n_control_wells >= 2).sum())
n_pass = int((cv.qc == "PASS").sum())
print(f"experiments (sample x timepoint)          : {n_total}")
print(f"  with a single control well (no CV)      : {n_total - n_multi}")
print(f"  with >= 2 control wells                 : {n_multi}")
print(f"  of those, passing the {QC_GATE:.0%} QC gate      : {n_pass}")
print(f"largest number of control wells anywhere  : {int(cv.n_control_wells.max())}")
print()
print(cv[["experiment", "n_control_wells", "dmso_cv", "qc"]].to_string(index=False))

passing = cv[cv.qc == "PASS"]
if len(passing):
    MED = float(passing.dmso_cv.median())
    SD_RATIO = 100 * MED * np.sqrt(2)
    print(f"\nmedian CV of passing experiments : {MED*100:.2f}%   (n = {len(passing)})")
    print(f"  -> SD of a single-well ratio    : {SD_RATIO:.2f} pp")
    print(f"  -> 2 SD                         : {2*SD_RATIO:.2f}% inhibition")
    print(f"  SAP quotes 10.8% -> 15.3 pp -> >30%")
else:
    MED = SD_RATIO = float("nan")
    print("\nNo experiment passes the QC gate — the CV route cannot be derived.")

# --------------------------------------------------------------------- derivation A
print()
print("=" * 78)
print("A. HALF-NORMAL FIT TO TREATED WELLS ABOVE 100% OF CONTROL")
print("=" * 78)
print("No agent can truly raise a slice's viability above its own DMSO control, so readings")
print("above 100% carry no signal and estimate the assay's noise directly.")
print()
print("CAVEAT, worth stating in the Methods: only wells whose true effect is near zero can")
print("land above 100%, so this fit is estimated from a censored tail. If most treated wells")
print("in the dataset had a genuine effect, few are left to fit and sigma is both imprecise")
print("and biased. Read the bootstrap interval and the well count below before quoting it.")
rng = np.random.default_rng(SEED)
failed = set(cv.loc[cv.qc.str.startswith("FAIL"), "experiment"])
treated = m[~m._is_control].copy()
treated["_expt_s"] = treated._expt.map(scrub)

hn_rows = []
for label, sub in [("all treated wells", treated),
                   ("QC-passing experiments only", treated[~treated._expt_s.isin(failed)])]:
    over = sub.loc[sub[C_MTS] > 100, C_MTS].dropna()
    d = (over - 100).to_numpy(float)
    sig = halfnormal_sigma(d)
    lo, hi = boot_sigma(d, rng)
    print(f"\n  {label}")
    print(f"    treated wells total          : {int(sub[C_MTS].notna().sum())}")
    print(f"    wells reading > 100%         : {len(d)}")
    if len(d):
        print(f"    excess over 100%, median     : {np.median(d):.1f} pp   max {d.max():.1f} pp")
        print(f"    half-normal sigma            : {sig:.2f} pp   "
              f"(bootstrap 95% CI {lo:.1f} to {hi:.1f})")
        print(f"    2 sigma                      : {2*sig:.2f}% inhibition")
    hn_rows.append(dict(subset=label, n_treated=int(sub[C_MTS].notna().sum()),
                        n_over_100=len(d), sigma_pp=sig, sigma_lo=lo, sigma_hi=hi,
                        cutoff_pct_inhibition=2 * sig))
print("\n  SAP quotes sigma = 16.2 pp -> >32%")

# --------------------------------------------------------------------- write + verdict
os.makedirs(args.outdir, exist_ok=True)
p1 = os.path.join(args.outdir, "threshold_cv_census_deid.csv")
p2 = os.path.join(args.outdir, "threshold_halfnormal_deid.csv")
cv.to_csv(p1, index=False)
pd.DataFrame(hn_rows).to_csv(p2, index=False)

print()
print("=" * 78)
print("VERDICT")
print("=" * 78)
qc_only = [r for r in hn_rows if r["subset"].startswith("QC-passing")][0]
a_cut, b_cut = qc_only["cutoff_pct_inhibition"], 2 * SD_RATIO
if np.isfinite(a_cut) and np.isfinite(b_cut):
    print(f"  half-normal route : {a_cut:5.1f}% inhibition   (SAP: 32%)")
    print(f"  DMSO CV route     : {b_cut:5.1f}% inhibition   (SAP: 30%)")
    print(f"  spread between the two routes: {abs(a_cut - b_cut):.1f} pp")
    if abs(a_cut - b_cut) <= 6:
        print("  -> The two routes still agree. 'Two independent derivations converged'")
        print("     stands, and the numbers above should replace the quoted ones in SAP 4.2")
        print("     if they differ.")
    else:
        print("  -> The two routes now DISAGREE by more than 6 pp. Do not repeat the")
        print("     'converged' claim in the Methods without re-examining both. Send this")
        print("     report back before the Methods paragraph is finalised.")
print(f"\nwrote {p1}")
print(f"wrote {p2}")
print("\nBoth files are de-identified aggregates and are safe to share back.")
