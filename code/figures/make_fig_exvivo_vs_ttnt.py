import pandas as pd, numpy as np
from scipy import stats
from matplotlib.lines import Line2D
from cns_style import apply_style, figure, palette, finalize_axes, save_figure, add_panel_labels

apply_style("Nature")
m = pd.read_csv('pairs_pfs.csv')
m = m[m.time.notna()].copy()
m['stop'] = m.start + m.time

CLASSES = sorted(m.drug_c.unique())
cols = dict(zip(CLASSES, palette("Nature", n=len(CLASSES))))
NOISE_BOUND = 68.3

fig, axes = figure(width_mm=183, height_mm=92, ncols=2)
axA, axB = axes

# ---------------- Panel A: swimmer, ordered by ex vivo viability
d = m.sort_values('pct_of_control', ascending=False).reset_index(drop=True)
for i, r in d.iterrows():
    grey = bool(r.flag)
    c = '0.78' if grey else cols[r.drug_c]
    axA.plot([r.start, r.stop], [i, i], linewidth=2.6, color=c,
             solid_capstyle='butt', zorder=2, alpha=1.0)
    if r.event == 1:
        axA.plot(r.stop, i, marker='|', markersize=5, color='black',
                 markeredgewidth=0.8, zorder=3)
    else:
        axA.plot(r.stop, i, marker='>', markersize=3.4, color='black', zorder=3)
XLO, XHI = -430, 330
offscale = d[(d.start < XLO)]
for i, r in offscale.iterrows():
    axA.annotate('', xy=(XLO, i), xytext=(XLO + 55, i),
                 arrowprops=dict(arrowstyle='-|>', color='0.55', linewidth=0.8,
                                 shrinkA=0, shrinkB=0))
    axA.text(XLO + 62, i, f'starts day {r.start:.0f}, off scale', fontsize=5.2,
             ha='left', va='center', color='0.45')
axA.set_xlim(XLO, XHI)
axA.axvline(0, color='black', linewidth=0.6, linestyle='-', zorder=1)
axA.text(4, len(d) - 0.35, 'Biopsy / ex vivo assay, day 0', fontsize=5.5,
         ha='left', va='center', color='black')
axA.set_yticks(range(len(d)))
axA.set_yticklabels([f"{r.study_id}  {r.drug_c}  ({r.pct_of_control:.0f}%)"
                     for _, r in d.iterrows()], fontsize=5.5)
axA.set_ylim(-0.7, len(d) - 0.1)
axA.invert_yaxis()
axA.set_xlabel('Days relative to biopsy')
axA.set_title('Matched systemic therapy courses, ordered by ex vivo viability')

# ---------------- Panel B: ex vivo viability vs time to next therapy
u = m[~m.flag]
for _, r in m.iterrows():
    grey = bool(r.flag)
    face = '0.78' if grey else (cols[r.drug_c] if r.started_after_biopsy else 'white')
    edge = '0.45' if grey else 'black'
    axB.plot(r.pct_of_control, r.time, 'o', markersize=4.4, markerfacecolor=face,
             markeredgecolor=edge, markeredgewidth=0.6, linestyle='none', zorder=3)
axB.axvline(100, color='0.35', linewidth=0.5, linestyle='--', zorder=1)
axB.axvline(NOISE_BOUND, color='0.35', linewidth=0.5, linestyle=':', zorder=1)
axB.set_xlabel('Ex vivo viability at 48 h (% of DMSO control)')
axB.set_ylabel('Time to next systemic therapy (days)')
axB.set_title('Time to next therapy by ex vivo viability')
rho, p = stats.spearmanr(u.pct_of_control, u.time)
rho_a, p_a = stats.spearmanr(m.pct_of_control, m.time)
axB.text(0.97, 0.97,
         f"Excluding tissue-exhausted\nSpearman $\\rho$ = {rho:.2f}, p = {p:.3f}, n = {len(u)}\n"
         f"All pairs\nSpearman $\\rho$ = {rho_a:.2f}, p = {p_a:.3f}, n = {len(m)}",
         transform=axB.transAxes, fontsize=5.8, va='top', ha='right')
axB.set_xlim(-8, 152); axB.set_ylim(0, 310)
axB.text(101, 6, 'DMSO control, 100%', fontsize=5.5, rotation=90, ha='left', va='bottom', color='0.35')
axB.text(NOISE_BOUND + 1, 6, f'Noise bound, {NOISE_BOUND}%', fontsize=5.5, rotation=90,
         ha='left', va='bottom', color='0.35')

handles = [Line2D([], [], color=cols[c], linewidth=2.6, label=c) for c in CLASSES]
handles += [
    Line2D([], [], color='0.78', linewidth=2.6, label='Tissue-exhausted specimen'),
    Line2D([], [], marker='|', color='black', linestyle='none', markersize=5,
           markeredgewidth=0.8, label='Next systemic line started'),
    Line2D([], [], marker='>', color='black', linestyle='none', markersize=3.4,
           label='No subsequent line recorded (censored)'),
    Line2D([], [], marker='o', linestyle='none', markersize=4.4, markerfacecolor='white',
           markeredgecolor='black', markeredgewidth=0.6,
           label='Course started before the biopsy (panel B, open)'),
]
for ax in (axA, axB): finalize_axes(ax)
fig.subplots_adjust(bottom=0.30, top=0.90, wspace=0.55, left=0.20, right=0.965)
fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 0.0), frameon=False,
           fontsize=6, ncol=2, handletextpad=0.5, labelspacing=0.4, columnspacing=1.4)
add_panel_labels(fig, axes, x=-0.02)
save_figure(fig, 'fig_exvivo_vs_ttnt', outdir='figures')
