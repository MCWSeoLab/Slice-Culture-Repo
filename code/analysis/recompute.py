"""Recompute all three analyses with patient-level clustering and small-sample inference."""
import pandas as pd, numpy as np, itertools, warnings
from scipy import stats
warnings.filterwarnings('ignore')
rng = np.random.default_rng(20260818)
B = 4000

def spearman(x, y):
    return stats.spearmanr(x, y)[0]

def cluster_bootstrap_rho(df, xcol, ycol, cluster='study_id', B=B):
    """Resample PATIENTS with replacement, not pairs. Pre-indexed for speed."""
    pats = df[cluster].unique()
    idx = [np.flatnonzero((df[cluster] == p).values) for p in pats]
    X = df[xcol].to_numpy(float); Y = df[ycol].to_numpy(float)
    out = np.empty(B); out[:] = np.nan
    picks = rng.integers(0, len(pats), size=(B, len(pats)))
    for b in range(B):
        sel = np.concatenate([idx[k] for k in picks[b]])
        x, y = X[sel], Y[sel]
        if len(np.unique(x)) < 2 or len(np.unique(y)) < 2: continue
        out[b] = stats.spearmanr(x, y)[0]
    ok = out[~np.isnan(out)]
    return (np.percentile(ok, 2.5), np.percentile(ok, 97.5), len(ok))

def one_per_patient(df, xcol, ycol, cluster='study_id'):
    """Exact enumeration over all ways of taking one independent pair per patient."""
    groups = [g.index.tolist() for _, g in df.groupby(cluster)]
    combos = list(itertools.product(*groups))
    rhos, ps = [], []
    for c in combos:
        d = df.loc[list(c)]
        if d[xcol].nunique() < 2 or d[ycol].nunique() < 2: continue
        r, p = stats.spearmanr(d[xcol], d[ycol]); rhos.append(r); ps.append(p)
    return np.array(rhos), np.array(ps), len(combos)

def perm_p(df, xcol, ycol, B=B):
    x = df[xcol].to_numpy(float); y = df[ycol].to_numpy(float)
    obs = stats.spearmanr(x, y)[0]
    rx = stats.rankdata(x); cnt = 0
    for _ in range(B):
        r = stats.spearmanr(rx, rng.permutation(y))[0]
        if abs(r) >= abs(obs) - 1e-12: cnt += 1
    return obs, (cnt + 1) / (B + 1)

def report(name, df, xcol, ycol):
    n, npat = len(df), df.study_id.nunique()
    if n < 4 or df[xcol].nunique() < 2:
        print(f"\n--- {name}: n={n} pairs / {npat} patients — too few to analyse"); return
    obs, pp = perm_p(df, xcol, ycol)
    lo, hi, nb = cluster_bootstrap_rho(df, xcol, ycol)
    rhos, ps, ncomb = one_per_patient(df, xcol, ycol)
    print(f"\n--- {name}")
    print(f"    n = {n} pairs from {npat} patients")
    print(f"    Spearman rho              = {obs:+.3f}")
    print(f"    permutation p (ignores clustering, anticonservative) = {pp:.4f}")
    print(f"    cluster bootstrap 95% CI  = [{lo:+.3f}, {hi:+.3f}]   ({nb}/{B} resamples usable)")
    if len(rhos):
        print(f"    one-pair-per-patient ({ncomb} combinations): rho median {np.median(rhos):+.3f}, "
              f"range [{rhos.min():+.3f}, {rhos.max():+.3f}], median p = {np.median(ps):.3f}")
    print(f"    CI includes 0: {'YES' if lo <= 0 <= hi else 'NO'}")

# ============================================================ A. binary matrix, lab's own calls
print("="*78); print("A. CONFUSION MATRIX ON THE LAB'S OWN Lab Response? CALLS")
print("   (retained only to show what the un-normalized metric yields; not the primary endpoint)")
ann = pd.read_csv('per_sheet_annotations_deid.csv').fillna('')
ann['drug'] = ann.which_matched.str.strip()
ann = ann[(ann.drug != '') & ann.lab_resp2.isin(['Yes','No'])]
ann = ann.sort_values(['study_id','drug','sheet','row'])          # stable order
mas = ann[ann.sheet == 'master'].drop_duplicates(subset=['study_id','drug'])
cat = ann[ann.sheet != 'master'].drop_duplicates(subset=['study_id','drug'])
keys = mas.set_index(['study_id','drug']).index
e = pd.concat([mas, cat[~cat.set_index(['study_id','drug']).index.isin(keys)]])
e = e[e.clinic_resp.isin(['Yes','No'])]                            # filter AFTER dedup
e['sens'] = (e.lab_resp2 == 'Yes').astype(int)
e['resp'] = (e.clinic_resp == 'Yes').astype(int)
TP=((e.sens==1)&(e.resp==1)).sum(); FP=((e.sens==1)&(e.resp==0)).sum()
FN=((e.sens==0)&(e.resp==1)).sum(); TN=((e.sens==0)&(e.resp==0)).sum()
print(f"   {len(e)} pairs / {e.study_id.nunique()} patients   TP={TP} FP={FP} FN={FN} TN={TN}")
print(f"   naive Fisher p = {stats.fisher_exact([[TP,FP],[FN,TN]])[1]:.3f}")
import statsmodels.api as sm
import statsmodels.formula.api as smf
try:
    g = smf.gee("resp ~ sens", groups="study_id", data=e,
                family=sm.families.Binomial(), cov_struct=sm.cov_struct.Exchangeable()).fit()
    b, se = g.params['sens'], g.bse['sens']
    print(f"   GEE (exchangeable, patient cluster): OR = {np.exp(b):.2f} "
          f"(95% CI {np.exp(b-1.96*se):.2f}-{np.exp(b+1.96*se):.2f}), p = {g.pvalues['sens']:.3f}")
    print(f"   -> clustering-corrected p is {'larger' if g.pvalues['sens']>stats.fisher_exact([[TP,FP],[FN,TN]])[1] else 'smaller'} than the naive Fisher p")
except Exception as ex:
    print("   GEE failed:", ex)

# ============================================================ B & C. normalized MTS analyses
pr = pd.read_csv('pairs_recist.csv')
pf = pd.read_csv('pairs_pfs.csv')
qc = pd.read_csv('qc_by_patient.csv').set_index('study_id').qc_verdict.to_dict()
for d in (pr, pf):
    d['qc'] = d.study_id.map(lambda s: qc.get(s, 'UNKNOWN'))

print("\n" + "="*78); print("B. EX VIVO 48h VIABILITY vs RECIST BEST RESPONSE (control-normalized)")
report("all pairs", pr, 'pct_of_control', 'recist_rank')
report("tissue-exhausted excluded", pr[~pr.flag], 'pct_of_control', 'recist_rank')
report("tissue-exhausted excluded AND DMSO-CV QC gate applied",
       pr[(~pr.flag) & (pr.qc != 'FAIL')], 'pct_of_control', 'recist_rank')

print("\n" + "="*78); print("C. EX VIVO 48h VIABILITY vs TIME TO NEXT SYSTEMIC THERAPY")
pf2 = pf[pf.time.notna()]
report("all pairs", pf2, 'pct_of_control', 'time')
report("tissue-exhausted excluded", pf2[~pf2.flag], 'pct_of_control', 'time')
report("tissue-exhausted excluded AND QC gate", pf2[(~pf2.flag)&(pf2.qc!='FAIL')], 'pct_of_control','time')
report("PROSPECTIVE ONLY (drug started at/after biopsy)",
       pf2[(~pf2.flag) & pf2.started_after_biopsy], 'pct_of_control', 'time')
