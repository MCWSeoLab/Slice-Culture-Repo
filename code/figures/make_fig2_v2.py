import itertools
import pandas as pd, numpy as np
from scipy import stats
from matplotlib.lines import Line2D
from cns_style import apply_style, figure, palette, finalize_axes, save_figure, add_panel_labels
apply_style("Nature")
SEED=20260818; B=4000

# Reads the de-identified twin, NOT the PHI workbook. Verified 9 Sep 2026 to reproduce the
# published rho and slice counts exactly (+0.743 / +0.651 / -0.548; 33 / 29 / 30 slices).
m=pd.read_csv('masterlist_linked_deid.csv').rename(columns={
    "mts_pct_of_control":"MTS_%","nuclear_density_pct":"NuclearDensity_%",
    "ki67_pct":"Ki67_%","caspase3_pct":"Caspase3_%",
    "sample_label":"Sample_ID","biopsy_site":"Biopsy_Site"})
for c in ["MTS_%","NuclearDensity_%","Ki67_%","Caspase3_%"]:
    m[c]=pd.to_numeric(m[c],errors='coerce')

SITES=["Pancreas","Peritoneum","Liver"]
cols=dict(zip(SITES,palette("Nature",n=3)))
MARK={"Pancreas":"o","Peritoneum":"^","Liver":"s"}

PANELS=[("NuclearDensity_%","Nuclear density (% of control)","left"),
        ("Ki67_%","Ki-67 (% positive cells)","left"),
        ("Caspase3_%","Cleaved caspase-3 (% positive cells)","right")]

def cluster_boot(d,x,y,cluster="Sample_ID"):
    ids=d[cluster].unique()
    idx=[np.flatnonzero((d[cluster]==i).values) for i in ids]
    X=d[x].to_numpy(float); Y=d[y].to_numpy(float)
    # own generator: a shared rng made each panel's CI depend on the panels before it
    picks=np.random.default_rng(SEED).integers(0,len(ids),size=(B,len(ids))); out=[]
    for b in range(B):
        sel=np.concatenate([idx[k] for k in picks[b]])
        xx,yy=X[sel],Y[sel]
        if len(np.unique(xx))<3 or len(np.unique(yy))<3: continue
        out.append(stats.spearmanr(xx,yy)[0])
    out=np.array(out)
    return np.percentile(out,2.5),np.percentile(out,97.5)

def block_perm(d,x,y,cluster="Sample_ID"):
    """Exact p permuting each biopsy's whole block of y between biopsies of equal size.

    Slices from one biopsy are not independent, so the biopsy - not the slice - is the unit
    that may be permuted. Only size-matched swaps are legal, so enumerate the permutations
    within each size class and take their product; never all G! orderings.
    """
    ids=list(d[cluster].unique())
    blk={i:d.loc[d[cluster]==i,y].to_numpy(float) for i in ids}
    X=d[x].to_numpy(float); obs=abs(stats.spearmanr(X,d[y].to_numpy(float))[0])
    cls={}
    for i in ids: cls.setdefault(len(blk[i]),[]).append(i)
    ks=sorted(cls); cnt=val=0
    for combo in itertools.product(*[itertools.permutations(cls[k]) for k in ks]):
        mp={}
        for k,pm in zip(ks,combo):
            for a,b in zip(cls[k],pm): mp[a]=b
        val+=1
        cnt+= abs(stats.spearmanr(X,np.concatenate([blk[mp[i]] for i in ids]))[0])>=obs-1e-12
    return cnt/val,val
# Cluster = the specimen (`sample_label`), which is what the manuscript counts as a biopsy.
# `study_id` is NOT the right unit here: the de-identification could not resolve one row pair
# and labels two distinct specimens - a primary pancreas and a liver met - `PT-007|PT-008`.
# Pooling them as one cluster is a sensitivity check, not the primary: it gives
# 0.0014 / 0.0312 / 0.0208 against 0.0007 / 0.0153 / 0.0203, i.e. the same conclusion.

fig,axes=figure(width_mm=174,height_mm=68,ncols=3)
for ax,(col,ylab,tpos) in zip(axes,PANELS):
    d=m.dropna(subset=["MTS_%",col]).copy()
    for s in SITES:
        g=d[d.Biopsy_Site==s]
        if not len(g): continue
        ax.plot(g["MTS_%"],g[col],MARK[s],markersize=3.8,markerfacecolor=cols[s],
                markeredgecolor='black',markeredgewidth=0.4,linestyle='none',zorder=3,label=s)
    # least-squares fit across all specimens
    sl,ic,_,_,_=stats.linregress(d["MTS_%"],d[col])
    xs=np.linspace(d["MTS_%"].min(),d["MTS_%"].max(),100)
    ax.plot(xs,sl*xs+ic,color='0.25',linewidth=0.7,zorder=2)
    r,_=stats.spearmanr(d["MTS_%"],d[col])
    lo,hi=cluster_boot(d,"MTS_%",col)
    p,nperm=block_perm(d,"MTS_%",col)
    ptxt=f"exact p = {p:.3f}"
    print(f"  {col:20s} rho={r:+.4f}  exact p={p:.4f} ({nperm} permutations)  "
          f"CI {lo:+.4f} to {hi:+.4f}  {len(d)} slices, {d.Sample_ID.nunique()} clusters")
    tx,ha=(0.03,'left') if tpos=='left' else (0.97,'right')
    ax.text(tx,0.98,f"Spearman $\\rho$ = {r:+.3f}\n{ptxt}\n95% CI {lo:+.2f} to {hi:+.2f}\n"
                    f"{len(d)} slices, {d.Sample_ID.nunique()} biopsies",
            transform=ax.transAxes,fontsize=5.6,va='top',ha=ha)
    ax.set_xlabel("Normalized MTS (% of control)")
    ax.set_ylabel(ylab)
    ax.set_xlim(-8,150)
axes[0].set_ylim(0,205)
axes[1].set_ylim(-0.8,19)
axes[2].set_ylim(-2,54)
handles=[Line2D([],[],marker=MARK[s],linestyle='none',markersize=3.8,markerfacecolor=cols[s],
                markeredgecolor='black',markeredgewidth=0.4,label=s) for s in SITES]
for ax in axes: finalize_axes(ax)
fig.subplots_adjust(bottom=0.30,top=0.95,wspace=0.34,left=0.075,right=0.985)
fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(0.5,0.0),frameon=False,
           fontsize=6,ncol=3,handletextpad=0.5,columnspacing=2.0)
add_panel_labels(fig,axes)
save_figure(fig,'fig2_mts_vs_histology_v2',outdir='figures')
