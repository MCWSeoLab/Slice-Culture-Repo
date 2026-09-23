import pandas as pd, numpy as np
import matplotlib.image as mpimg
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from course_panel import draw_course
from cns_style import apply_style, palette, finalize_axes, save_figure, MM_PER_INCH, plt
apply_style("Nature")

# `stains` lists the H&E tiles available for that case, in the order the MTS bars are drawn:
# (display label, file stem under fig4_images/). Only conditions that were actually
# photographed appear — the two cases were run against different agent panels.
# `cut` truncates the clinical-course panel and sets the top of its linear y axis.
#   PT-024  day 154 is her CA 19-9 nadir (73 U/mL) on gemcitabine/nab-paclitaxel. Draws after
#           it (87, 95, 160 U/mL at days 182-238) are NOT shown. This is a chosen cut, not a
#           clinical event, and the legend says so.
#   PT-015  day 245 is the day tovorafenib was added to cobimetinib. Everything after is held
#           for a later report. Because the cut is that agent's start day, draw_course drops
#           it from the treatment strip too, so it does not appear in the figure at all.
#           Her day -353 draw of 18 002 U/mL sits above the break.
CASES=[dict(sid='PT-024', ihc='Liver bx. [DATE]', title='Case 1 — KRAS G12D liver metastasis',
            sub='gemcitabine/nab-paclitaxel matched, begun 18 days after biopsy', key='Gem/nab-Pac',
            cut=dict(day_max=154, cap=620),
            stains=[('DMSO','PT024_HE_DMSO'), ('FOLFIRINOX','PT024_HE_FOLFIRINOX'),
                    ('Gem/nab-Pac','PT024_HE_GemPac'), ('RMC-7977','PT024_HE_RMC7977'),
                    ('Trametinib','PT024_HE_Trametinib')]),
       dict(sid='PT-015', ihc='Liver Bx.7', title='Case 2 — BRAF V487 liver metastasis',
            sub='MEK inhibitor matched', key='Trametinib',
            cut=dict(day_max=245, cap=300),
            stains=[('DMSO','PT015_HE_DMSO'), ('FOLFIRINOX','PT015_HE_FOLFIRINOX'),
                    ('Trametinib','PT015_HE_Trametinib'), ('Vemurafenib','PT015_HE_Vemurafenib')])]

DRUGCOL=dict(zip(['FOLFIRINOX','Gem/nab-Pac','RMC-7977','Trametinib','Vemurafenib'],
                 palette("Nature",n=5)))
# ordered greyscale so clinical response never collides with the drug palette
RESPCOL={'CR':'0.20','PR':'0.35','SD':'0.62','PD':'0.88','':'1.0'}

def canon(s):
    s=str(s).lower().replace(' ','').replace('-','').replace('/','')
    if 'folfirinox' in s: return 'FOLFIRINOX'
    if 'gem' in s and ('nab' in s or 'pac' in s): return 'Gem/nab-Pac'
    if 'rmc' in s: return 'RMC-7977'
    if 'trametinib' in s: return 'Trametinib'
    if 'vemuraf' in s: return 'Vemurafenib'
    return None
# `px` (mts_normalized_final.csv) and `ml` (Masterlist.xlsx) were read here but never used.
# Dropping them removes this figure's only dependency on a PHI file, so it now builds from
# de-identified inputs alone.
mllink=pd.read_csv('masterlist_linked_deid.csv'); mllink['drug']=mllink.treatment.map(canon)

ihc=pd.read_csv('ihc_tidy.csv')
gg=ihc.groupby(['biopsy','marker','treatment']).value.mean().reset_index()
ref=gg[gg.treatment=='DMSO'].set_index(['biopsy','marker']).value.rename('dmso')
gg=gg.join(ref,on=['biopsy','marker']); gg=gg[(gg.treatment!='DMSO')&gg.dmso.gt(0)].copy()
gg['fold']=gg.value/gg.dmso; gg['drug']=gg.treatment.map(canon)
MARKERS=['Cleaved caspase-3','Ki-67','CD45','EpCAM']

rg=pd.read_csv('regimens_merged.csv').fillna('')
for c in ['day_start_rel_biopsy','day_end_rel_biopsy']: rg[c]=pd.to_numeric(rg[c],errors='coerce')
EXCL=r'whipple|\bbx\b|biopsy|sbrt|chemort|cgy|\bfx\b|resect|NED|LAB'

W=183/MM_PER_INCH; H=236/MM_PER_INCH
fig=plt.figure(figsize=(W,H))
# Row 3 now holds three stacked axes per case (treatment strip / break region / trace),
# so it needs more height than it did as a single log-scaled axes.
gs=GridSpec(4,2,figure=fig,height_ratios=[0.88,1.06,0.88,1.62],hspace=0.60,wspace=0.28,
            left=0.085,right=0.985,top=0.938,bottom=0.135)
letters=iter("ABCDEFGH")
STAIN_X=[(0.048,0.494),(0.552,0.994)]   # x-extent of each case's stain band

for col,C in enumerate(CASES):
    sid=C['sid']
    # ---- MTS by drug
    a=fig.add_subplot(gs[0,col])
    d48=mllink[(mllink.study_id==sid)&(mllink.timepoint_hr==48)&mllink.drug.notna()]
    order=[x for x in ['FOLFIRINOX','Gem/nab-Pac','RMC-7977','Trametinib','Vemurafenib']
           if x in set(d48.drug)]
    for i,dr in enumerate(order):
        v=d48[d48.drug==dr].mts_pct_of_control.mean()
        hatch='//' if dr==C['key'] else None
        a.bar(i,v,width=0.64,color=DRUGCOL[dr],edgecolor='black',linewidth=0.5,zorder=2,hatch=hatch)
        a.text(i,v+3,f"{v:.0f}",ha='center',va='bottom',fontsize=5.6)
        if dr==C['key']:
            a.text(i,v+17,'received',ha='center',va='bottom',fontsize=5.0,style='italic',color='0.30')
            a.annotate('',xy=(i,v+10),xytext=(i,v+16.5),
                       arrowprops=dict(arrowstyle='-|>',color='0.30',lw=0.6,shrinkA=0,shrinkB=0))
    a.axhline(100,color='0.35',lw=0.5,ls='--',zorder=1)
    a.axhline(70,color='0.35',lw=0.5,ls=':',zorder=1)
    a.set_xticks(range(len(order))); a.set_xticklabels(order,rotation=28,ha='right',fontsize=5.6)
    a.set_xlim(-0.62,len(order)-0.38)
    a.set_ylim(0,150); a.set_ylabel('Ex vivo viability at 48 h\n(% of DMSO control)',fontsize=6.2)
    a.set_title(f"{C['title']}\n{C['sub']}",fontsize=7,pad=6)
    if col==0:
        a.text(len(order)-0.4,103,'DMSO control',fontsize=5.2,ha='right',va='bottom',color='0.35')
        a.text(len(order)-0.4,73,'30% inhibition',fontsize=5.2,ha='right',va='bottom',color='0.35')
    finalize_axes(a); a.text(-0.20,1.24,next(letters),transform=a.transAxes,fontsize=9,fontweight='bold',va='top')

    # ---- representative H&E at 48 h, one tile per agent tested
    # Placed by hand rather than through a subgridspec: imshow preserves the 4:3 aspect, so a
    # single row of five tiles across half a page leaves each one tiny. Wrapping to two rows
    # and sizing the tiles explicitly makes the histology legible at print size; reading order
    # is still left to right, top to bottom, matching the bar order in the panel above.
    st=C['stains']
    cell=gs[1,col].get_position(fig)
    BX0,BX1=STAIN_X[col]
    NCOL=3 if len(st)>4 else 2
    NROW=int(np.ceil(len(st)/NCOL))
    GX,GY=0.007,0.012
    imh,imw=mpimg.imread(f"fig4_images/{st[0][1]}.png").shape[:2]
    # one tile width for the whole figure, set by the busiest row, so the two cases are
    # shown at the same size on the page and neither reads as a different magnification
    tw=(BX1-BX0-GX*2)/3
    th=tw*(imh/imw)*(W/H)
    if NROW*th+(NROW-1)*GY > cell.height:                 # height-bound: back off the width
        th=(cell.height-(NROW-1)*GY)/NROW
        tw=th/((imh/imw)*(W/H))
    ytop=cell.y0+cell.height-(NROW*th+(NROW-1)*GY)/2 - (cell.height-(NROW*th+(NROW-1)*GY))/2
    ytop=cell.y0+cell.height - (cell.height-(NROW*th+(NROW-1)*GY))/2
    for k,(lab,stem) in enumerate(st):
        rr,cc=divmod(k,NCOL)
        nthis=min(NCOL,len(st)-rr*NCOL)                   # centre a short final row
        rowx=BX0+((BX1-BX0)-(nthis*tw+(nthis-1)*GX))/2
        b=fig.add_axes([rowx+cc*(tw+GX), ytop-(rr+1)*th-rr*GY, tw, th])
        b.imshow(mpimg.imread(f"fig4_images/{stem}.png"))
        b.set_xticks([]); b.set_yticks([])
        for sp in b.spines.values():
            sp.set_visible(True)                          # the house style hides top/right
            sp.set_linewidth(1.3 if lab==C['key'] else 0.4)
            sp.set_color('black' if lab==C['key'] else '0.45')
        b.set_title(lab,fontsize=6.4,pad=2,
                    fontweight='bold' if lab==C['key'] else 'normal')
    fig.text(BX0-0.052,ytop+0.014,next(letters),fontsize=9,fontweight='bold',va='top')

    # ---- IHC fold change
    a=fig.add_subplot(gs[2,col])
    sub=gg[gg.biopsy==C['ihc']]
    drugs=[x for x in order if x in set(sub.drug)]
    wdt=0.8/max(len(drugs),1)
    for j,dr in enumerate(drugs):
        xs=np.arange(len(MARKERS))+(j-(len(drugs)-1)/2)*wdt
        vals=[sub[(sub.drug==dr)&(sub.marker==m)].fold.mean() for m in MARKERS]
        a.bar(xs,vals,width=wdt*0.9,color=DRUGCOL[dr],edgecolor='black',linewidth=0.4,
              zorder=2,label=dr,hatch='//' if dr==C['key'] else None)
    a.axhline(1.0,color='0.35',lw=0.5,ls='--',zorder=1)
    a.set_xticks(range(len(MARKERS))); a.set_xticklabels(MARKERS,rotation=22,ha='right',fontsize=5.6)
    a.set_ylabel('IHC fold change\nvs DMSO',fontsize=6.2)
    finalize_axes(a); a.text(-0.20,1.16,next(letters),transform=a.transAxes,fontsize=9,fontweight='bold',va='top')

    # ---- clinical course: CA 19-9 on a linear axis, treatment strip above it
    cp = draw_course(fig, gs[3, col], sid, DRUGCOL, show_ylabel=(col == 0), **C['cut'])
    finalize_axes(cp['ax_lo'], tight=False)
    if cp['ax_hi'] is not None:
        finalize_axes(cp['ax_hi'], tight=False)
        cp['ax_hi'].spines['bottom'].set_visible(False)
        cp['ax_hi'].tick_params(axis='x', length=0, labelbottom=False)
    cp['ax_strip'].text(-0.20, 1.02, next(letters), transform=cp['ax_strip'].transAxes,
                        fontsize=9, fontweight='bold', va='top')

h1=[Line2D([],[],marker='s',ls='none',ms=5,mfc=DRUGCOL[d],mec='black',mew=0.4,label=d)
    for d in ['FOLFIRINOX','Gem/nab-Pac','RMC-7977','Trametinib','Vemurafenib']]
h2=[Line2D([],[],marker='s',ls='none',ms=5,mfc='white',mec='black',mew=0.6,label='Hatched / arrow = agent the patient received'),
    Line2D([],[],marker='s',ls='none',ms=5,mfc='0.65',mec='black',mew=0.4,label='Other agent (labelled in panel)'),
    Line2D([],[],marker='o',ls='-',ms=3.0,mfc='white',mec='#1f3b5c',mew=0.8,color='#1f3b5c',lw=0.9,label='Serum CA 19-9 (linear scale)'),
    Line2D([],[],ls='--',color='0.55',lw=0.6,label='CA 19-9 upper limit of normal (37.5 U/mL)'),
    Line2D([],[],ls=(0,(3,2)),color='0.45',lw=0.7,label='End of the period shown'),
    Line2D([],[],ls='none',marker=r'$\mathrm{PD}$',ms=8,color='0.3',label='RECIST assessment at each scan')]
fig.legend(handles=h1+h2,loc='lower center',bbox_to_anchor=(0.5,0.002),frameon=False,fontsize=5.8,
           ncol=5,handletextpad=0.5,columnspacing=1.3)
fig.text(0.085, 0.088,
         "H&E panels: representative slices after 48 h, 10x; the outlined condition is the "
         "agent the patient went on to receive.\n"
         "Course panels are truncated at the dashed vertical rule. Case 1 is shown to the "
         "CA 19-9 nadir at day 154; three later draws, rising to 160 U/mL by day 238, are not "
         "displayed.\n"
         "Case 2 is shown to day 245, when a second targeted agent was added to "
         "cobimetinib; that agent and the period after it are not shown. "
         "A bar with an arrowed end denotes a regimen still running where the panel ends.",
         fontsize=5.2, va="top", linespacing=1.5)
save_figure(fig,'fig4_two_cases',outdir='figures')
