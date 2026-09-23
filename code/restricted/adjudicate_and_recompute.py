#!/usr/bin/env python3
"""Final concordance analysis: chart abstraction + the three response adjudications.

ADJUDICATIONS from investigator chart review, applied here explicitly so they are auditable:
  1. PT-025 cobimetinib   -> NOT scored. The day +8 PD scan reflects pre-existing
                             progression, not a response to a 6-day course. No pair added.
  2. PT-015 cobimetinib   -> PR, not SD. CT at +84 and +182, both PR, on monotherapy
                             before tovorafenib was added at +245.
  3. PT-014 gem/nab       -> PD, as recorded. The +78 SD is superseded by PD at +139/+166.
"""
import pandas as pd, numpy as np
from scipy import stats
from statsmodels.stats.proportion import proportion_confint as ci

rng=np.random.default_rng(20260818); B=4000
BEN={'CR','PR','SD'}; RANK={'PD':0,'SD':1,'PR':2,'CR':3}; CUT=30.0

p=pd.read_csv('pairs_chart_updated.csv')
p.loc[(p.study_id=='PT-015')&(p.drug=='MEK inhibitor'),'response']='PR'      # adjudication 2
p['benefit']=p.response.isin(BEN); p['rank']=p.response.map(RANK)
p['sens']=p.mts_pct_of_control<(100-CUT)
p.to_csv('pairs_final.csv',index=False)
u=p[~p.flag].copy()

print("FINAL MATCHED PAIRS")
print(u[['study_id','drug','mts_pct_of_control','regimen','day','response','benefit','sens']]
      .to_string(index=False,float_format=lambda v:f"{v:.1f}"))
print(f"\n{len(u)} pairs from {u.study_id.nunique()} patients (tissue-exhausted excluded); "
      f"{len(p)} / {p.study_id.nunique()} including them")

def f(k,n):
    if n==0: return "n/a"
    lo,hi=ci(k,n,method='beta'); return f"{100*k/n:.0f}% (95% CI {100*lo:.0f}-{100*hi:.0f})"

TP=int((u.sens&u.benefit).sum()); FP=int((u.sens&~u.benefit).sum())
FN=int((~u.sens&u.benefit).sum()); TN=int((~u.sens&~u.benefit).sum())
pv=stats.fisher_exact([[TP,FP],[FN,TN]])[1]
print(f"\n2x2 at >{CUT:.0f}% inhibition, binary clinical benefit")
print(f"   TP {TP}   FP {FP}   FN {FN}   TN {TN}")
for nm,k,n in [("sensitivity",TP,TP+FN),("specificity",TN,TN+FP),("PPV",TP,TP+FP),
               ("NPV",TN,TN+FN),("accuracy",TP+TN,len(u))]:
    print(f"   {nm:<12} {f(k,n)}")
print(f"   Fisher exact p = {pv:.3f}")

def boot(df,x,y):
    ids=df.study_id.unique(); idx=[np.flatnonzero((df.study_id==i).values) for i in ids]
    X=df[x].to_numpy(float); Y=df[y].to_numpy(float)
    picks=rng.integers(0,len(ids),size=(B,len(ids))); out=[]
    for b in range(B):
        sel=np.concatenate([idx[k] for k in picks[b]]); xx,yy=X[sel],Y[sel]
        if len(np.unique(xx))<2 or len(np.unique(yy))<2: continue
        out.append(stats.spearmanr(xx,yy)[0])
    out=np.array(out); return stats.spearmanr(X,Y)[0],np.percentile(out,2.5),np.percentile(out,97.5)

print("\nORDINAL RECIST, co-primary")
for lab,d in [("all pairs",p),("tissue-exhausted excluded",u)]:
    r,lo,hi=boot(d,'mts_pct_of_control','rank')
    print(f"   {lab:<28} n={len(d):>2}/{d.study_id.nunique()}  rho {r:+.2f}  cluster-bootstrap 95% CI {lo:+.2f} to {hi:+.2f}")

a=u[u.benefit].mts_pct_of_control; b=u[~u.benefit].mts_pct_of_control
print(f"\nViability by benefit: median {a.median():.0f}% (n={len(a)}) vs {b.median():.0f}% (n={len(b)}), "
      f"Mann-Whitney p = {stats.mannwhitneyu(a,b)[1]:.3f}")

print("\nTHRESHOLD SENSITIVITY")
print(f"{'cutoff':>7} {'TP':>3}{'FP':>3}{'FN':>3}{'TN':>3} {'accuracy':>9} {'Fisher p':>9}")
for c in [20,25,30,32,40,50]:
    s=u.mts_pct_of_control<(100-c)
    t=int((s&u.benefit).sum()); fp=int((s&~u.benefit).sum())
    fn=int((~s&u.benefit).sum()); tn=int((~s&~u.benefit).sum())
    print(f"{c:>6}% {t:>3}{fp:>3}{fn:>3}{tn:>3} {100*(t+tn)/len(u):>8.0f}% "
          f"{stats.fisher_exact([[t,fp],[fn,tn]])[1]:>9.3f}")

pr=u[u.day>=0]
print(f"\nPROSPECTIVE SUBSET (begun at/after biopsy): {len(pr)} pairs, {pr.study_id.nunique()} patients")
print(pr[['study_id','drug','mts_pct_of_control','day','response','sens','benefit']]
      .to_string(index=False,float_format=lambda v:f"{v:.1f}"))
print("   concordant:", int((pr.sens==pr.benefit).sum()), "of", len(pr))
