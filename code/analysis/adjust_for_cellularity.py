#!/usr/bin/env python3
"""Can the concordance analysis be adjusted for tumour content? Four candidate approaches,
each run and reported honestly, including the ones that fail.

Throughout, patients are the resampling unit (4000 cluster bootstrap resamples, seed as in
the SAP) and no naive p-value is reported for anything computed over drug-patient pairs.
"""
import numpy as np, pandas as pd
from scipy import stats

SEED, B = 20260818, 4000
rng = np.random.default_rng(SEED)

tc = pd.read_csv("tumour_content_deid.csv").dropna(subset=["study_id"]).set_index("study_id")
p = pd.read_csv("pairs_final.csv")
p = p[~p.flag].copy()
p["cell"] = p.study_id.map(tc.epcam_pct_mean)
p = p.dropna(subset=["cell"])
p["sens"] = p.mts_pct_of_control < 70
print(f"{len(p)} pairs, {p.study_id.nunique()} patients, "
      f"{p.cell.nunique()} distinct cellularity values\n")


def cboot(df, fn):
    """cluster bootstrap over patients; returns (point, lo, hi, n_usable)"""
    ids = df.study_id.unique()
    idx = {i: np.flatnonzero((df.study_id == i).values) for i in ids}
    out = []
    for _ in range(B):
        sel = np.concatenate([idx[i] for i in rng.choice(ids, len(ids), replace=True)])
        v = fn(df.iloc[sel])
        if v is not None and np.isfinite(v):
            out.append(v)
    out = np.asarray(out)
    return fn(df), np.percentile(out, 2.5), np.percentile(out, 97.5), len(out)


def spear(d, a, b):
    if d[a].nunique() < 2 or d[b].nunique() < 2: return None
    return stats.spearmanr(d[a], d[b])[0]


print("=" * 74)
print("STEP 0.  Is tumour content actually a confounder?")
print("=" * 74)
print("A confounder must be associated with BOTH the exposure (viability) and the")
print("outcome (clinical benefit). If it is associated only with viability, it is a")
print("measurement-quality problem, not a confounder, and adjustment is the wrong tool.\n")

r, lo, hi, _ = cboot(p, lambda d: spear(d, "cell", "mts_pct_of_control"))
print(f"  cellularity ~ viability      rho {r:+.2f}   95% CI {lo:+.2f} to {hi:+.2f}")
r, lo, hi, _ = cboot(p, lambda d: spear(d, "cell", "rank"))
print(f"  cellularity ~ RECIST rank    rho {r:+.2f}   95% CI {lo:+.2f} to {hi:+.2f}")
ben, pd_ = p[p.benefit].cell, p[~p.benefit].cell
print(f"  cellularity by outcome       benefit median {ben.median():.1f}% (n={len(ben)}), "
      f"PD median {pd_.median():.1f}% (n={len(pd_)})")
print("  NOTE: only 2 pairs, from 2 patients, had progressive disease. Any statement")
print("        about the outcome side of this rests on those two observations.\n")

print("=" * 74)
print("APPROACH 1.  Partial rank correlation, holding cellularity constant")
print("=" * 74)


def partial(d):
    """Spearman partial correlation of viability with RECIST rank given cellularity."""
    if min(d.mts_pct_of_control.nunique(), d["rank"].nunique(), d.cell.nunique()) < 3:
        return None
    R = np.corrcoef(np.vstack([stats.rankdata(d.mts_pct_of_control),
                               stats.rankdata(d["rank"]), stats.rankdata(d.cell)]))
    try:
        P = np.linalg.inv(R)
    except np.linalg.LinAlgError:
        return None
    den = np.sqrt(P[0, 0] * P[1, 1])
    return -P[0, 1] / den if den > 0 else None


raw, rlo, rhi, _ = cboot(p, lambda d: spear(d, "mts_pct_of_control", "rank"))
par, plo, phi, nb = cboot(p, partial)
print(f"  unadjusted   rho {raw:+.2f}   95% CI {rlo:+.2f} to {rhi:+.2f}")
print(f"  partial      rho {par:+.2f}   95% CI {plo:+.2f} to {phi:+.2f}   "
      f"({nb}/{B} resamples estimable)")
print("  Verdict: estimable, but cellularity is constant within a patient, so partialling")
print("  it out removes most of the between-patient variation the correlation rests on.")
print("  The interval is uninformative. Report as a sensitivity analysis at most.\n")

print("=" * 74)
print("APPROACH 2.  Restrict to biopsies above a tumour-content floor")
print("=" * 74)
for thr in [0, 10, 20, 25, 30]:
    q = p[p.cell >= thr]
    if len(q) < 3: continue
    TP = int((q.sens & q.benefit).sum()); FP = int((q.sens & ~q.benefit).sum())
    FN = int((~q.sens & q.benefit).sum()); TN = int((~q.sens & ~q.benefit).sum())
    sens = f"{TP}/{TP+FN}" if TP + FN else "n/a"
    spec = f"{TN}/{TN+FP}" if TN + FP else "n/a"
    acc = (TP + TN) / len(q)
    print(f"  >={thr:2d}% tumour: {len(q):2d} pairs / {q.study_id.nunique()} patients   "
          f"TP {TP} FP {FP} FN {FN} TN {TN}   sens {sens}  spec {spec}  acc {acc:.0%}")
print("  Verdict: the floor improves every metric, which is exactly why it cannot be")
print("  presented as a primary result — the threshold was chosen after seeing the")
print("  confound. Defensible only as a declared post-hoc sensitivity analysis, and as")
print("  the pre-specified eligibility rule for the prospective study.\n")

print("=" * 74)
print("APPROACH 3.  Mechanistic correction to a tumour-compartment viability")
print("=" * 74)
print("  If measured signal = f*V_tumour + (1-f)*V_other and the drug spares non-tumour")
print("  tissue, then V_tumour = (V_measured - (1 - f)) / f.\n")
p["f"] = p.cell / 100.0
p["v"] = p.mts_pct_of_control / 100.0
p["v_tumour"] = (p.v - (1 - p.f)) / p.f
show = p[["study_id", "drug", "cell", "mts_pct_of_control", "v_tumour"]].copy()
show["v_tumour"] = (100 * show.v_tumour).round(0)
show.columns = ["patient", "drug", "tumour %", "measured %", "corrected %"]
print(show.sort_values("tumour %").to_string(index=False))
bad = int(((p.v_tumour < -0.5) | (p.v_tumour > 2)).sum())
print(f"\n  {bad} of {len(p)} corrected values fall outside -50% to 200%.")
print("  Verdict: the correction divides by f, so at 4-7% tumour it multiplies the noise")
print("  by 15-20x and returns impossible viabilities. It also assumes cytotoxics spare")
print("  stroma and hepatocytes, which is false. Not usable.\n")

print("=" * 74)
print("APPROACH 4.  Change the endpoint to a tumour-restricted one")
print("=" * 74)
print("  Rather than adjusting the whole-slice metabolic readout, score drug effect only")
print("  inside EpCAM-positive regions — cleaved caspase-3 or Ki-67 within tumour. This")
print("  removes the dilution at source instead of correcting for it afterwards.")
ihc = pd.read_csv("ihc_tidy.csv")
have = (ihc[ihc.marker.isin(["Cleaved caspase-3", "Ki-67"])]
        .groupby("biopsy").treatment.nunique())
print(f"\n  Biopsies with CC3/Ki-67 across >1 condition: {int((have > 1).sum())} of {len(have)}")
print("  Verdict: the right answer, but it needs region-restricted rescoring in QuPath")
print("  (co-registered EpCAM mask), which is new image analysis, not a reanalysis of")
print("  numbers already in hand.\n")

print("=" * 74)
print("RECOMMENDATION")
print("=" * 74)
print("""  Do not adjust the primary analysis. With 6 patients, 2 of whom supply every
  progressive-disease observation, no adjustment can separate tumour content from
  drug sensitivity, and each method above either destroys the estimate (1), is
  circular (2), is mathematically unstable (3), or requires new data (4).

  Report instead:
    - tumour content as a measured property of the cohort (Supplementary Figure 1)
    - the >=20% subset as one clearly-labelled post-hoc sensitivity analysis
    - a pre-specified minimum tumour content and a tumour-restricted endpoint as
      design requirements for the prospective study""")
