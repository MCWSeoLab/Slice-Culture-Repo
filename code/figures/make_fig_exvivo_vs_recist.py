import pandas as pd, numpy as np
from scipy import stats
from matplotlib.lines import Line2D
from cns_style import apply_style, figure, palette, finalize_axes, save_figure, add_panel_labels

apply_style("Nature")
p = pd.read_csv('pairs_recist.csv')
ORDER   = ['PD','SD','PR','CR']
present = [r for r in ORDER if (p.response == r).any()]
xpos    = {r: i for i, r in enumerate(present)}

CLASSES = sorted(p.drug_c.unique())
cols    = dict(zip(CLASSES, palette("Nature", n=len(CLASSES))))
CLASSMATCH_TRI = 'class-level (trametinib ex vivo / cobimetinib clinical)'
NOISE_BOUND = 68.3          # 100 - 1.96*sigma; sigma = 16.2 pp (control-normalized MTS > 100%)
YMAX = 158

def dodge(df):
    """Deterministic horizontal offsets within each response category."""
    out = []
    for resp, g in df.groupby('response'):
        g = g.sort_values('pct_of_control')
        k = len(g)
        offs = np.zeros(1) if k == 1 else np.linspace(-0.20, 0.20, k)
        for off, (_, r) in zip(offs, g.iterrows()):
            out.append((xpos[resp] + off, r))
    return out

def draw(ax, df, grey_flagged, label_lines=False):
    for x, r in dodge(df):
        flagged = bool(r.flag)
        face = '0.78' if (grey_flagged and flagged) else cols[r.drug_c]
        edge = '0.45' if (grey_flagged and flagged) else 'black'
        mk   = '^' if r.class_match == CLASSMATCH_TRI else 'o'
        ax.plot(x, r.pct_of_control, mk, markersize=4.2, markerfacecolor=face,
                markeredgecolor=edge, markeredgewidth=0.5, linestyle='none', zorder=3)
    ax.axhline(100, color='0.35', linewidth=0.5, linestyle='--', zorder=1)
    ax.axhline(NOISE_BOUND, color='0.35', linewidth=0.5, linestyle=':', zorder=1)
    ax.set_xticks(range(len(present))); ax.set_xticklabels(present)
    ax.set_xlim(-0.5, len(present) - 0.5)
    ax.set_ylim(-8, YMAX)
    ax.set_xlabel('Clinical best response (RECIST)')
    ax.set_ylabel('Ex vivo viability at 48 h\n(% of DMSO control)')
    # label the shared reference lines once, on the right-hand panel only
    if label_lines:
        xr = len(present) - 0.5
        ax.text(xr + 0.08, 100, 'DMSO control, 100%', fontsize=5.5, ha='left',
                va='center', color='0.35', clip_on=False)
        ax.text(xr + 0.08, NOISE_BOUND, f'Assay noise bound, {NOISE_BOUND}%', fontsize=5.5,
                ha='left', va='center', color='0.35', clip_on=False)

fig, axes = figure(width_mm=183, height_mm=76, ncols=2)
axA, axB = axes

u = p[~p.flag]
draw(axA, u, grey_flagged=False, label_lines=False)
rho, pv = stats.spearmanr(u.pct_of_control, u.recist_rank)
axA.set_title('Tissue-exhausted specimens excluded')
axA.text(0.02, 0.03, f"Spearman $\\rho$ = {rho:.2f}, p = {pv:.2f}\n"
                     f"n = {len(u)} drug–patient pairs, {u.study_id.nunique()} patients",
         transform=axA.transAxes, fontsize=6, va='bottom', ha='left')

draw(axB, p, grey_flagged=True, label_lines=True)
rho2, pv2 = stats.spearmanr(p.pct_of_control, p.recist_rank)
axB.set_title('All pairs; tissue-exhausted specimens in grey')
axB.text(0.02, 0.03, f"Spearman $\\rho$ = {rho2:.2f}, p = {pv2:.2f}\n"
                     f"n = {len(p)} drug–patient pairs, {p.study_id.nunique()} patients",
         transform=axB.transAxes, fontsize=6, va='bottom', ha='left')

handles = [Line2D([], [], marker='o', linestyle='none', markersize=4.2,
                  markerfacecolor=cols[c], markeredgecolor='black',
                  markeredgewidth=0.5, label=c) for c in CLASSES]
handles += [
    Line2D([], [], marker='o', linestyle='none', markersize=4.2, markerfacecolor='white',
           markeredgecolor='black', markeredgewidth=0.5, label='Same agent ex vivo and clinically'),
    Line2D([], [], marker='^', linestyle='none', markersize=4.2, markerfacecolor='white',
           markeredgecolor='black', markeredgewidth=0.5,
           label='Class-level match: trametinib ex vivo, cobimetinib clinically'),
    Line2D([], [], marker='o', linestyle='none', markersize=4.2, markerfacecolor='0.78',
           markeredgecolor='0.45', markeredgewidth=0.5, label='Tissue-exhausted specimen'),
]
for ax in (axA, axB): finalize_axes(ax)
fig.subplots_adjust(bottom=0.26, top=0.90, wspace=0.34, left=0.085, right=0.815)
fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 0.0),
           frameon=False, fontsize=6, ncol=2, handletextpad=0.4,
           labelspacing=0.4, columnspacing=1.4)
add_panel_labels(fig, axes)
save_figure(fig, 'fig_exvivo_vs_recist', outdir='figures')
