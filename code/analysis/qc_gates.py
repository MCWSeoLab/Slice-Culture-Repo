"""Per-experiment QC gating for the ex vivo MTS data."""
import pandas as pd, numpy as np

CV_LIMIT = 0.20          # exclude an experiment whose DMSO replicate CV exceeds this
MIN_REPS = 1             # conditions with fewer wells than this are flagged (see report)

raw = pd.read_csv('prism_raw_long_deid.csv').fillna('')
mts = raw[raw.file.str.contains('MTS', case=False) | raw.sheet.str.contains('MTS', case=False)].copy()
CTRL = ['dmso','control','ctrl','vehicle','untreated']
mts['is_ctrl'] = mts.column.str.strip().str.lower().isin(CTRL)

# ---- DMSO replicate CV per experiment (file x sheet x timepoint)
c = mts[mts.is_ctrl].groupby(['file','sheet','row_label','study_id'],dropna=False).value.agg(['count','mean','std']).reset_index()
c['dmso_cv'] = c['std']/c['mean']
c['qc'] = np.where(c['count'] < 2, 'UNKNOWN (single control well)',
           np.where(c.dmso_cv > CV_LIMIT, 'FAIL (CV > 20%)', 'PASS'))
c = c.rename(columns={'count':'n_control_wells','row_label':'timepoint'})
c.to_csv('qc_dmso_cv.csv', index=False)
print("=== QC 1: DMSO replicate CV per experiment ===")
print(c[['study_id','sheet','timepoint','n_control_wells','mean','dmso_cv','qc']]
      .to_string(index=False, float_format=lambda x: f"{x:.3f}"))

# ---- negative absorbance
neg = mts[mts.value < 0]
print(f"\n=== QC 2: negative absorbance readings: {len(neg)} of {len(mts)} ({len(neg)/len(mts)*100:.1f}%) ===")
print("   rule: floored at 0 before normalization (a negative reading is background over-subtraction,")
print("         not negative viability); flagged in the per-condition table.")

# ---- replication depth per treated condition
n = pd.read_csv('prism_mts_normalized_deid.csv')
n = n[n.file.str.contains('MTS', case=False)]
rep = n.groupby(['file','sheet','timepoint_label','drug']).size().rename('n_wells').reset_index()
print(f"\n=== QC 3: replication depth, treated conditions ===")
print(rep.n_wells.value_counts().sort_index().rename('conditions').to_frame().to_string())
print(f"   single-well conditions: {(rep.n_wells==1).sum()} of {len(rep)} ({(rep.n_wells==1).mean()*100:.0f}%)")

# ---- which patients pass / fail / unknown
pat = c.groupby('study_id').qc.agg(lambda s: 'FAIL' if (s=='FAIL (CV > 20%)').any()
                                   else ('PASS' if (s=='PASS').any() else 'UNKNOWN'))
pat = pat[pat.index.astype(str).str.startswith('PT-')]
print("\n=== QC verdict by patient (only where replicate-level control data exists) ===")
print(pat.to_frame('qc_verdict').to_string())
pat.to_frame('qc_verdict').to_csv('qc_by_patient.csv')
