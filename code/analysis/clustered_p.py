#!/usr/bin/env python3
"""What happens to the CR/PR vs SD/PD p-value once patients, not pairs, are the unit.

The naive Mann-Whitney over 5 vs 5 drug-patient pairs gives p = 0.032. Ten pairs from six
patients are not ten independent observations, so that number is not interpretable on its own.
Three cluster-aware versions are computed here:

  1. one-pair-per-patient exact enumeration  — every way of picking one course per patient
  2. patient-level aggregation               — one summary value per patient, exact MWU
  3. cluster permutation                     — permute the outcome label across PATIENTS,
                                               recompute the pair-level statistic, exact

Run: python clustered_p.py
"""
import itertools
import math

import numpy as np
import pandas as pd
from scipy import stats

CUT = 70.0
p = pd.read_csv("pairs_final.csv")
a = p[~p.flag].copy()                       # analysis set: tissue-exhausted excluded
a["responder"] = a.response.isin(["CR", "PR"])
a["benefit"] = a.benefit.astype(bool)

print("=" * 78)
print("THE DATA")
print("=" * 78)
print(a[["study_id", "drug", "response", "mts_pct_of_control", "responder", "benefit"]]
      .sort_values(["study_id", "drug"]).to_string(index=False))

for LABEL, COLM in [("OBJECTIVE RESPONSE  (CR/PR vs SD/PD)", "responder"),
                    ("AS PUBLISHED  (CR/PR/SD vs PD)", "benefit")]:
    print()
    print("=" * 78)
    print(LABEL)
    print("=" * 78)

    # --- is the outcome constant within each patient? ------------------------------------
    byp = a.groupby("study_id")[COLM].agg(nuniq="nunique", size="size", first="first")
    mixed = byp[byp.nuniq > 1]
    npos_pt = int((byp.nuniq.eq(1) & byp["first"]).sum())
    nneg_pt = int((byp.nuniq.eq(1) & ~byp["first"]).sum())
    print(f"patients: {len(byp)}   pairs: {len(a)}")
    print(f"patients whose courses are ALL positive: {npos_pt}")
    print(f"patients whose courses are ALL negative: {nneg_pt}")
    print(f"patients with a mix of both outcomes:    {len(mixed)}"
          + (f"  ({', '.join(mixed.index)})" if len(mixed) else ""))

    pos = a[a[COLM]].mts_pct_of_control
    neg = a[~a[COLM]].mts_pct_of_control
    naive = stats.mannwhitneyu(pos, neg, alternative="two-sided")[1]
    print(f"\n  naive MWU over pairs ({len(pos)} vs {len(neg)}):  p = {naive:.4f}"
          "   <- NOT interpretable, pairs are clustered")

    # --- 1. one pair per patient, every combination ---------------------------------------
    groups = [g.index.to_list() for _, g in a.groupby("study_id")]
    combos = list(itertools.product(*groups))
    ps, ok = [], 0
    for c in combos:
        s = a.loc[list(c)]
        pp, nn = s[s[COLM]].mts_pct_of_control, s[~s[COLM]].mts_pct_of_control
        if len(pp) < 1 or len(nn) < 1:
            continue
        pv = stats.mannwhitneyu(pp, nn, alternative="two-sided")[1]
        ps.append(pv); ok += pv < 0.05
    ps = np.asarray(ps)
    print(f"\n  1. one pair per patient, {len(combos)} combinations")
    print(f"     p ranges {ps.min():.4f} to {ps.max():.4f}, median {np.median(ps):.4f}")
    print(f"     subsets reaching p < 0.05:  {ok} / {len(ps)}")

    # --- 2. patient-level aggregation ------------------------------------------------------
    # One value per patient. Patients with a mixed outcome cannot be assigned to a side, so
    # this is only defined when the outcome is constant within every patient.
    if len(mixed) == 0:
        agg = a.groupby(["study_id", COLM]).mts_pct_of_control.mean().reset_index()
        ap, an = agg[agg[COLM]].mts_pct_of_control, agg[~agg[COLM]].mts_pct_of_control
        pv = stats.mannwhitneyu(ap, an, alternative="two-sided")[1]
        floor = 2 / (math.comb(len(ap) + len(an), len(ap)))
        print(f"\n  2. patient-level means, exact MWU ({len(ap)} vs {len(an)} patients)")
        print(f"     medians {ap.median():.1f}% vs {an.median():.1f}%   p = {pv:.4f}")
        print(f"     smallest two-sided p ATTAINABLE at {len(ap)} vs {len(an)}: {floor:.4f}")
    else:
        print(f"\n  2. patient-level aggregation not defined — {len(mixed)} patient(s) "
              "contribute courses to both outcome groups")

    # --- 3. cluster permutation ------------------------------------------------------------
    if len(mixed) == 0:
        ids = byp.index.to_list()
        lab = {i: bool(byp.loc[i, "first"]) for i in ids}
        k = sum(lab.values())
        obs = stats.mannwhitneyu(pos, neg, alternative="two-sided")[0]
        # centre the statistic so "more extreme" is symmetric
        n1, n2 = len(pos), len(neg)
        obs_c = abs(obs - n1 * n2 / 2)
        cnt = tot = 0
        for sel in itertools.combinations(ids, k):
            m = a.study_id.isin(sel)
            pp, nn = a[m].mts_pct_of_control, a[~m].mts_pct_of_control
            if len(pp) < 1 or len(nn) < 1:
                continue
            u = stats.mannwhitneyu(pp, nn, alternative="two-sided")[0]
            tot += 1
            cnt += abs(u - len(pp) * len(nn) / 2) >= obs_c - 1e-9
        print(f"\n  3. cluster permutation: outcome label permuted across PATIENTS")
        print(f"     {tot} distinct assignments ({k} of {len(ids)} patients labelled positive)")
        print(f"     exact p = {cnt}/{tot} = {cnt / tot:.4f}")
        print(f"     smallest attainable p = 1/{tot} = {1 / tot:.4f}")
    else:
        print("\n  3. cluster permutation not defined — the outcome is not constant "
              "within every patient, so it cannot be permuted at the patient level")

print()
print("=" * 78)
print("BOTTOM LINE")
print("=" * 78)
print("""Under CR/PR vs SD/PD the outcome is perfectly nested inside patient: every patient's
courses share one outcome. The comparison is therefore 3 patients versus 3 patients, not
5 pairs versus 5. At 3 vs 3 the smallest two-sided p any rank test can return is 0.10, so
p < 0.05 is unreachable however cleanly the groups separate. The naive p = 0.032 is an
artifact of counting ten courses as ten independent observations.""")
