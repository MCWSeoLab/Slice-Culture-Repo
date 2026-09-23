import pandas as pd, numpy as np
from scipy import stats
from cns_style import apply_style, figure, palette, finalize_axes, save_figure
rng=np.random.default_rng(20260818); B=4000
apply_style("Nature")

pr=pd.read_csv('pairs_recist.csv'); pf=pd.read_csv('pairs_pfs.csv')
qc=pd.read_csv('qc_by_patient.csv').set_index('study_id').qc_verdict.to_dict()
for d in (pr,pf): d['qc']=d.study_id.map(lambda s: qc.get(s,'UNKNOWN'))
pf=pf[pf.time.notna()]

def boot(df,x,y):
    pats=df.study_id.unique()
    idx=[np.flatnonzero((df.study_id==p).values) for p in pats]
    X=df[x].to_numpy(float); Y=df[y].to_numpy(float)
    picks=rng.integers(0,len(pats),size=(B,len(pats))); out=[]
    for b in range(B):
        sel=np.concatenate([idx[k] for k in picks[b]])
        xx,yy=X[sel],Y[sel]
        if len(np.unique(xx))<2 or len(np.unique(yy))<2: continue
        out.append(stats.spearmanr(xx,yy)[0])
    out=np.array(out)
    return stats.spearmanr(X,Y)[0], np.percentile(out,2.5), np.percentile(out,97.5)

rows=[]
for lbl,df,x,y in [
    ("RECIST · all pairs",                         pr,'pct_of_control','recist_rank'),
    ("RECIST · tissue-exhausted excluded",         pr[~pr.flag],'pct_of_control','recist_rank'),
    ("RECIST · + DMSO-CV QC gate",                 pr[(~pr.flag)&(pr.qc!='FAIL')],'pct_of_control','recist_rank'),
    ("TTNT · all pairs",                           pf,'pct_of_control','time'),
    ("TTNT · tissue-exhausted excluded",           pf[~pf.flag],'pct_of_control','time'),
    ("TTNT · + DMSO-CV QC gate",                   pf[(~pf.flag)&(pf.qc!='FAIL')],'pct_of_control','time'),
]:
    r,lo,hi=boot(df,x,y)
    rows.append(dict(label=lbl,rho=r,lo=lo,hi=hi,n=len(df),npat=df.study_id.nunique()))
R=pd.DataFrame(rows)
R.to_csv('forest_rho.csv',index=False)

fig,ax=figure(width_mm=140,height_mm=78)
cols=palette("Nature",n=2)
ypos=np.arange(len(R))[::-1]
for yv,(_,r) in zip(ypos,R.iterrows()):
    c=cols[0] if r.label.startswith('RECIST') else cols[1]
    ax.plot([r.lo,r.hi],[yv,yv],color=c,linewidth=1.1,solid_capstyle='round',zorder=2)
    ax.plot(r.rho,yv,'o',markersize=4.4,markerfacecolor=c,markeredgecolor='black',
            markeredgewidth=0.5,zorder=3)
    ax.text(1.07,yv,f"{r.rho:+.2f}  [{r.lo:+.2f}, {r.hi:+.2f}]",fontsize=5.6,va='center',ha='left')
ax.axvline(0,color='0.35',linewidth=0.5,linestyle='--',zorder=1)
ax.set_yticks(ypos)
ax.set_yticklabels([f"{r.label}\n{int(r.n)} pairs, {int(r.npat)} patients" for _,r in R.iterrows()],
                   fontsize=5.8)
ax.set_xlim(-1.05,1.05); ax.set_ylim(-0.7,len(R)-0.3)
ax.set_xlabel("Spearman $\\rho$ (95% CI, patient-level cluster bootstrap)")
ax.set_title("Ex vivo 48 h viability vs clinical endpoints")
finalize_axes(ax)
fig.subplots_adjust(left=0.345,right=0.745,top=0.89,bottom=0.17)
save_figure(fig,'fig_rho_forest',outdir='figures')
print(R.to_string(index=False,float_format=lambda v:f"{v:.3f}"))
