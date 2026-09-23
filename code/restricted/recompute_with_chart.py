#!/usr/bin/env python3
"""Recompute the concordance analysis using the completed chart abstraction.

Chart data (PT-024, PT-015, PT-014, PT-025) supersedes the workbook's regimen rows for
those four patients; the remaining patients keep their Masterlist-derived courses.
Matching rules are unchanged from SAP v1.1: drug-class level, FOLFOX/FOLFIRI NOT matched
to FOLFIRINOX, course closest in time to the biopsy where an agent was given twice.
"""
import pandas as pd, numpy as np
from scipy import stats
from statsmodels.stats.proportion import proportion_confint as ci

rng = np.random.default_rng(20260818); B = 4000
CUT = 30.0                       # >30% inhibition = ex vivo sensitive
BEN = {"CR", "PR", "SD"}         # binary clinical benefit (SAP v1.1, post hoc)
RANK = {"PD": 0, "SD": 1, "PR": 2, "CR": 3}

def canon(s):
    s = str(s).lower().replace(" ", "").replace("-", "").replace("/", "")
    if "folfiri" in s and "nox" not in s: return None      # FOLFIRI is NOT FOLFIRINOX
    if "folfirinox" in s or "folforinox" in s: return "FOLFIRINOX"
    if "gem" in s and ("nab" in s or "pac" in s or "abrax" in s): return "Gem/nab-Pac"
    if "rmc" in s: return "RMC-7977"
    if "trametinib" in s or "cobimetinib" in s: return "MEK inhibitor"
    if "vemuraf" in s: return "Vemurafenib"
    return None

# ---- ex vivo, 48 h, % of DMSO control
ml = pd.read_csv("masterlist_linked_deid.csv")
ml = ml[ml.timepoint_hr == 48].copy()
ml["drug"] = ml.treatment.map(canon)
ev = (ml.dropna(subset=["drug"])
        .groupby(["study_id", "sample_label", "drug"], as_index=False)
        .mts_pct_of_control.mean())
FLAG = {"PT-023"}                                   # tissue-exhausted (lab's own note)

# ---- clinical courses
CH = pd.read_csv("chart_therapy_deid.csv")
CH["drug"] = CH.regimen.map(canon)
CH = CH.rename(columns={"day_start": "day", "best_response": "response"})
CH["src"] = "chart"
chart_pts = set(CH.study_id)

OLD = pd.read_csv("regimens_deid.csv")
OLD = OLD[~OLD.study_id.isin(chart_pts)].copy()
OLD["drug"] = OLD.regimen.map(canon)
OLD = OLD.rename(columns={"day_start_rel_biopsy": "day", "response": "response"})
OLD["src"] = "workbook"; OLD["day_progression"] = np.nan; OLD["day_stop"] = OLD.day_end_rel_biopsy

cl = pd.concat([CH[["study_id","drug","regimen","day","day_stop","response","day_progression","src"]],
                OLD[["study_id","drug","regimen","day","day_stop","response","day_progression","src"]]],
               ignore_index=True)
cl = cl.dropna(subset=["drug", "day"])
cl = cl[cl.response.notna() & cl.response.astype(str).str.upper().isin(RANK)]
cl["response"] = cl.response.str.upper()

# course closest in time to the biopsy when an agent was given more than once
cl["absday"] = cl.day.abs()
cl = cl.sort_values("absday").drop_duplicates(subset=["study_id", "drug"], keep="first")

# ---- pair up
pairs = ev.merge(cl, on=["study_id", "drug"], how="inner")
pairs["flag"] = pairs.study_id.isin(FLAG)
pairs["sens"] = pairs.mts_pct_of_control < (100 - CUT)
pairs["benefit"] = pairs.response.isin(BEN)
pairs["rank"] = pairs.response.map(RANK)
pairs = pairs.sort_values(["study_id", "drug"]).reset_index(drop=True)
pairs.to_csv("pairs_chart_updated.csv", index=False)

print("MATCHED PAIRS AFTER CHART ABSTRACTION")
print(pairs[["study_id","drug","mts_pct_of_control","regimen","day","response",
             "benefit","sens","flag","src"]].to_string(index=False,
      float_format=lambda v: f"{v:.1f}"))

u = pairs[~pairs.flag]
print(f"\nall pairs {len(pairs)} from {pairs.study_id.nunique()} patients | "
      f"unflagged {len(u)} from {u.study_id.nunique()} patients")

def tab(d, label):
    TP=int((d.sens&d.benefit).sum()); FP=int((d.sens&~d.benefit).sum())
    FN=int((~d.sens&d.benefit).sum()); TN=int((~d.sens&~d.benefit).sum())
    def f(k,n):
        if n==0: return "n/a"
        lo,hi=ci(k,n,method="beta"); return f"{100*k/n:>3.0f}% ({100*lo:.0f}-{100*hi:.0f})"
    p=stats.fisher_exact([[TP,FP],[FN,TN]])[1]
    print(f"\n{label}: TP={TP} FP={FP} FN={FN} TN={TN}   n={len(d)} pairs, {d.study_id.nunique()} patients")
    print(f"   sensitivity {f(TP,TP+FN)}   specificity {f(TN,TN+FP)}   PPV {f(TP,TP+FP)}"
          f"   NPV {f(TN,TN+FN)}   accuracy {f(TP+TN,len(d))}   Fisher p = {p:.3f}")
    return TP,FP,FN,TN

tab(u, "2x2, binary clinical benefit, tissue-exhausted excluded")
tab(pairs, "2x2, all pairs including tissue-exhausted")

def boot(df, x, y):
    ids=df.study_id.unique(); idx=[np.flatnonzero((df.study_id==i).values) for i in ids]
    X=df[x].to_numpy(float); Y=df[y].to_numpy(float)
    picks=rng.integers(0,len(ids),size=(B,len(ids))); out=[]
    for b in range(B):
        sel=np.concatenate([idx[k] for k in picks[b]]); xx,yy=X[sel],Y[sel]
        if len(np.unique(xx))<2 or len(np.unique(yy))<2: continue
        out.append(stats.spearmanr(xx,yy)[0])
    out=np.array(out)
    return stats.spearmanr(X,Y)[0], np.percentile(out,2.5), np.percentile(out,97.5)

print("\nORDINAL RECIST (co-primary, SAP v1.1)")
for lab,d in [("all pairs",pairs),("tissue-exhausted excluded",u)]:
    r,lo,hi=boot(d,"mts_pct_of_control","rank")
    print(f"   {lab:<28} n={len(d):>2}/{d.study_id.nunique()}  rho {r:+.2f}  cluster-bootstrap 95% CI {lo:+.2f} to {hi:+.2f}")

mw = stats.mannwhitneyu(u[u.benefit].mts_pct_of_control, u[~u.benefit].mts_pct_of_control)[1]
print(f"\nViability by benefit (unflagged): median {u[u.benefit].mts_pct_of_control.median():.0f}% "
      f"(n={u.benefit.sum()}) vs {u[~u.benefit].mts_pct_of_control.median():.0f}% (n={(~u.benefit).sum()})"
      f",  Mann-Whitney p = {mw:.3f}")
