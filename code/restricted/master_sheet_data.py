#!/usr/bin/env python3
"""Assemble the per-biopsy record set for the master data sheet.
Every value is DERIVED from a source file in deid/ (or, for the identifier block,
from PHI_CROSSWALK locally). Nothing is retyped by hand."""
import io, os, re, math
from pathlib import Path
import pandas as pd

ROOT = Path(os.environ['HOME']) / 'mnt' / 'Organotypic Slice Culture'
DEID = ROOT / 'deid'
CWFILE = ROOT / 'PHI_CROSSWALK' / 'PHI_CROSSWALK_DO_NOT_SHARE.csv'

def rd(name, **kw):
    return pd.read_csv(DEID / name, keep_default_na=False, **kw)

def rd_commented(name):
    txt = [l for l in open(DEID / name) if not l.startswith('#')]
    return pd.read_csv(io.StringIO(''.join(txt)), keep_default_na=False)

SPINE = [
    ('PDAC Bx 1',       'Primary PDAC', 'PT-001', 'PDAC 1'),
    ('Liver Bx 1',      'Liver',        'PT-002', 'L.Bx 1'),
    ('Liver Bx 2',      'Liver',        'PT-003', 'L.Bx 2'),
    ('PDAC Bx 2',       'Primary PDAC', 'PT-004', 'PDAC 2'),
    ('Liver Bx 3',      'Liver',        'PT-005', 'L.Bx 3'),
    ('PDAC Bx 3',       'Primary PDAC', 'PT-006', 'PDAC 3'),
    ('Liver Bx 4',      'Liver',        'PT-008', 'L.Bx 4'),
    ('PDAC Bx 4',       'Primary PDAC', 'PT-007', 'PDAC 4'),
    ('Peritoneum Bx 1', 'Peritoneum',   'PT-009', 'Met Peritonium Bx-1'),
    ('PDAC Bx 5',       'Primary PDAC', 'PT-011', 'PDAC 5'),
    ('Liver Bx 5',      'Liver',        'PT-012', 'L.Bx 5'),
    ('Peritoneum Bx 2', 'Peritoneum',   'PT-013', 'Met Peritonium Bx-2'),
    ('Liver Bx 6',      'Liver',        'PT-014', 'L.Bx 6'),
    ('Liver Bx 7',      'Liver',        'PT-015', 'L.Bx 7'),
    ('PDAC Bx 6',       'Primary PDAC', 'PT-017', 'PDAC-6'),
    ('PDAC Bx 7',       'Primary PDAC', 'PT-019', 'PDAC 7'),
    ('Liver Bx 8',      'Liver',        'PT-023', 'L.Bx 9'),
    ('Liver Bx 9',      'Liver',        'PT-024', 'L.Bx 10'),
    ('Liver Bx 10',     'Liver',        'PT-025', 'L.Bx 11'),
]
ID_INFERRED = {'Liver Bx 4', 'PDAC Bx 4'}
DRUGS = ['FOLFIRINOX', 'Gem/nab-Pac', 'RMC-7977', 'RMC-6236', 'Trametinib',
         'Vemurafenib', 'Staurosporine']
# RMC-6236 (daraxonrasib) added 1 Sep 2026 - the agent going forward. No historical data:
# no cohort biopsy was treated with it, so its block is empty by design, not by a lookup miss.
DRUG_LABEL = {'RMC-7977': 'RMC-7977 (pan-RAS, preclinical)',
              'RMC-6236': 'RMC-6236 (daraxonrasib, pan-RAS)'}

IHC_BIOPSY = {'Liver Bx.3 IC50': 'Liver Bx 3', 'Liver Bx.4': 'Liver Bx 4',
              'Liver Bx.5': 'Liver Bx 5', 'Liver Bx.6': 'Liver Bx 6',
              'Liver Bx.7': 'Liver Bx 7', 'Liver bx. [DATE]': 'Liver Bx 9',
              'Liver Bx. [DATE]': 'Liver Bx 10', 'PDAC 4': 'PDAC Bx 4',
              'PDAC 5': 'PDAC Bx 5', 'Umbilicus': 'Peritoneum Bx 2'}
ML_SAMPLE = {'PDAC 1': 'PDAC Bx 1', 'PDAC 2 IC50': 'PDAC Bx 2', 'Liver 3': 'Liver Bx 3',
             'PDAC 4': 'PDAC Bx 4', 'Liver 4 resection 1': 'Liver Bx 4', 'PDAC 5': 'PDAC Bx 5',
             'Liver 5': 'Liver Bx 5', 'Umbilicus  1': 'Peritoneum Bx 2',
             'Liver 6 resection 2': 'Liver Bx 6', 'Liver 7': 'Liver Bx 7',
             'Liver 8': 'Liver Bx 8', 'Liver 9': 'Liver Bx 9', 'Liver 10': 'Liver Bx 10'}
ML_DRUG = {'FOLFIRINOX': 'FOLFIRINOX', 'Gem/nab': 'Gem/nab-Pac', 'RMC7977': 'RMC-7977',
           'Trametinib': 'Trametinib', 'Vemurafenib': 'Vemurafenib',
           'STS': 'Staurosporine', 'Day 7 STS': 'Staurosporine'}
PAIR_DRUG = {'FOLFIRINOX': 'FOLFIRINOX', 'Gem/nab-Pac': 'Gem/nab-Pac',
             'MEK inhibitor': 'Trametinib'}
ERK_SAMPLE = {'KRAS G12D': 'Liver Bx 6', 'KRAS G12R': 'Liver Bx 10'}

t1 = rd_commented('table1_source_deid.csv')
t2 = rd_commented('table2_source_deid.csv')
t3 = rd_commented('table3_assaylog_source_deid.csv')
ml = rd('masterlist_linked_deid.csv')
wells = rd('mts_treated_pct_deid.csv')
dmso = rd('mts_dmso_cv_wells_deid.csv')
ihc = rd('ihc_tidy.csv')
tc = rd('tumor_content_deid.csv')
erk = rd('erk_tidy.csv')
pairs = rd('pairs_final.csv')
regs = rd('regimens_merged.csv')
ther = rd('chart_therapy_deid.csv')
mark = rd('chart_markers_deid.csv')
fup = rd('chart_followup_deid.csv')
sf2nd = rd('suppfig2_peritoneal_nuclear_density_deid.csv')
sf1nd = rd('suppfig1_pdac_nuclear_density_deid.csv')

ml['sample'] = ml.sample_label.map(ML_SAMPLE)
ml['arm'] = ml.treatment.map(ML_DRUG)
ihc['sample'] = ihc.biopsy.map(IHC_BIOPSY)
tc['sample'] = tc.biopsy.map(IHC_BIOPSY)
pairs['arm'] = pairs.drug.map(PAIR_DRUG)
pairs['sample'] = pairs.sample_label.map(ML_SAMPLE)
t1i = t1.set_index('sample'); t2i = t2.set_index('sample'); t3i = t3.set_index('sample')
cw = pd.read_csv(CWFILE, dtype=str, keep_default_na=False).set_index('study_id')

def num(x):
    if x is None or x == '' or (isinstance(x, float) and math.isnan(x)):
        return None
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(v) else v

def rnd(x, n=2):
    v = num(x)
    return None if v is None else round(v, n)

KRAS_RE = re.compile(r'\b((?:K|N)RAS)\s+(G12[DVRC]|Q61[HRLK]|G13[DC])')

def kras_of(genomic):
    m = KRAS_RE.search(genomic)
    if m:
        return f'{m.group(1)} {m.group(2)}'
    if 'KRAS wild-type' in genomic:
        return 'KRAS wild-type'
    return 'Not detected / not reported'

DOSE_AS = re.compile(r'Doses as (Liver Bx \d+)')

def build():
    recs = {}
    for i, (samp, grp, sid, lab) in enumerate(SPINE, start=1):
        r2 = t2i.loc[samp]; r1 = t1i.loc[samp]; r3 = t3i.loc[samp]
        r = {'n': i, 'sample': samp, 'group': grp, 'study_id': sid, 'lab_id': lab,
             'study_day': int(r2.study_day)}
        assert int(r1.age_years) == int(r2.age_years) and r1.sex == r2.sex, samp
        c = cw.loc[sid]
        r['mrn'] = c.MRN; r['tb'] = c.TB_number
        r['initials'] = c.Patient_Initials; r['bx_date'] = c.Biopsy_Date
        r['id_note'] = ('study_id INFERRED: PT-007 and PT-008 share this biopsy date. '
                        'Assigned from the PT-008 chart note naming liver metastases. CONFIRM.'
                        if samp in ID_INFERRED else '')
        r['age'] = int(r1.age_years); r['sex'] = r1.sex
        r['stage'] = r1.stage; r['grade'] = r1.grade
        r['approach'] = r1.approach; r['prior_sys'] = r1.prior_systemic
        r['histology'] = r2.histologic_type
        r['genomic'] = r2.genomic; r['kras'] = kras_of(r2.genomic)
        r['met_site'] = r2.met_site
        r['ca199'] = rnd(r1.ca19_9, 1)
        r['ca199_censored'] = ('Yes - reported as >100000, not a measurement'
                               if (num(r1.ca19_9) or 0) >= 100000 else '')
        r['ca125'] = rnd(r1.ca125, 1); r['cea'] = rnd(r1.cea, 1)
        r['assays'] = r3.assay; r['assay_details'] = r3.details
        r['tissue_type'] = 'Primary' if grp == 'Primary PDAC' else 'Metastatic'
        r['bx_site'] = {'Liver': 'Liver', 'Peritoneum': 'Peritoneum',
                        'Primary PDAC': 'Pancreas'}[grp]
        r['exhausted'] = 'Yes' if 'exhaust' in r3.details.lower() else ''
        recs[samp] = r

    w = wells.copy()
    w['pct'] = pd.to_numeric(w.pct_of_control, errors='coerce')
    g = w.groupby(['sample', 'drug', 'timepoint'])
    agg = g['pct'].agg(['mean', 'count']).reset_index()
    conc = g['concentration'].agg(lambda s: '; '.join(sorted({x for x in s if x.strip()}))).reset_index().set_index(['sample', 'drug', 'timepoint'])
    for _, row in agg.iterrows():
        r = recs.get(row['sample'])
        if r is None:
            continue
        tp = row['timepoint']
        r['mts_%s_%s' % (row.drug, tp)] = round(row['mean'], 2)
        if tp == '48h' or ('mts_n_%s' % row.drug) not in r:
            r['mts_n_%s' % row.drug] = int(row['count'])
        c = conc.loc[(row['sample'], row['drug'], tp), 'concentration']
        if c and not r.get('dose_%s' % row.drug):
            r['dose_%s' % row.drug] = c
    for samp, r in recs.items():
        m = DOSE_AS.search(r['assay_details'])
        if m:
            src = recs.get(m.group(1))
            for d in DRUGS:
                if not r.get('dose_%s' % d) and src and src.get('dose_%s' % d):
                    r['dose_%s' % d] = '%s (assay log: "%s")' % (src['dose_%s' % d], m.group(0))

    for _, row in ml.iterrows():
        r = recs.get(row['sample'])
        if r is None or not row['arm']:
            continue
        a = row['arm']
        r['mlmts_%s' % a] = rnd(row.mts_pct_of_control, 2)
        r['nd_%s' % a] = rnd(row.nuclear_density_pct, 2)
        if row.notes:
            r['ml_notes'] = '; '.join([x for x in [r.get('ml_notes', ''), '%s: %s' % (a, row.notes)] if x])
        tp = num(row.timepoint_hr)
        if tp and int(tp) not in (24, 48):
            r['mts_other_%s' % a] = '%s%% of control at %d h' % (rnd(row.mts_pct_of_control, 2), int(tp))

    ihc['v'] = pd.to_numeric(ihc.value, errors='coerce')
    MK = {'Cleaved caspase-3': 'cc3', 'Ki-67': 'ki67', 'CD45': 'cd45', 'EpCAM': 'epcam'}
    for (samp, mk, tr), sub in ihc.groupby(['sample', 'marker', 'treatment']):
        r = recs.get(samp)
        if r is None or mk not in MK:
            continue
        r['ihc_%s_%s' % (MK[mk], tr)] = round(sub['v'].mean(), 2)
        r['ihc_n_%s' % tr] = int(sub['v'].count())

    d = dmso.copy()
    for _, row in d.iterrows():
        r = recs.get(row.get('sample', ''))
        if r is None:
            continue
        tp = row.get('timepoint', '')
        if tp in ('24h', '48h', '72h'):
            n = num(row.get('n_control_wells'))
            if n:
                r['dmso_n_%s' % tp] = int(n)
            cv = num(row.get('dmso_cv'))
            if cv is not None:
                r['dmso_cv_%s' % tp] = round(cv * 100, 1) if cv <= 1.5 else round(cv, 1)

    for _, row in tc.iterrows():
        r = recs.get(row['sample'])
        if r is None:
            continue
        r['tc_mean'] = rnd(row.epcam_pct_mean, 2)
        r['tc_sd'] = rnd(row.epcam_pct_sd, 2)
        n = num(row.n_slices)
        r['tc_n'] = int(n) if n else None

    for df, samp, grpname in ((sf2nd, 'Peritoneum Bx 2', 'DMSO'), (sf1nd, 'PDAC Bx 2', 'Day 7')):
        sub = df[df.group == grpname]
        v = pd.to_numeric(sub.value, errors='coerce').dropna()
        if len(v):
            recs[samp]['nd_abs_dmso'] = round(v.mean(), 1)
            recs[samp]['nd_abs_note'] = 'mean of %d slices, group "%s"' % (len(v), grpname)

    erk['v'] = pd.to_numeric(erk.value, errors='coerce')
    for (al, tr, mkr), sub in erk.groupby(['allele', 'treatment', 'marker']):
        samp = ERK_SAMPLE.get(al)
        if not samp:
            continue
        k = 'perk' if mkr == 'p-ERK' else 'terk'
        recs[samp]['erk_%s_%s' % (k, 'dmso' if tr == 'DMSO' else 'rmc')] = round(sub['v'].mean(), 2)

    for _, row in pairs.iterrows():
        r = recs.get(row['sample'])
        if r is None or not row['arm']:
            continue
        a = row['arm']
        r['clinreg_%s' % a] = row.regimen
        r['clinresp_%s' % a] = row.response
        r['n_pairs'] = r.get('n_pairs', 0) + 1
        r['in_concordance'] = 'No - tissue exhausted' if str(row.flag) == 'True' else 'Yes'

    def lines_for(sid):
        t = ther[ther.study_id == sid]
        if len(t):
            out = []
            for _, x in t.iterrows():
                ds = num(x.day_start)
                if ds is None or 'resection' in x.regimen.lower():
                    continue
                out.append((ds, x.regimen, num(x.day_stop), x.best_response, num(x.day_progression)))
            return sorted(out), 'chart review'
        g = regs[(regs.study_id == sid) & (regs.src_sheet == 'master')]
        out = []
        for _, x in g.iterrows():
            ds = num(x.day_start_rel_biopsy)
            if ds is None or re.search(r'\bBX\b|BIOPSY', x.regimen, re.I):
                continue
            out.append((ds, x.regimen, num(x.day_end_rel_biopsy), x.response, num(x.day_progression)))
        return sorted(out), 'lab workbook'

    for samp, r in recs.items():
        lines, src = lines_for(r['study_id'])
        r['tx_source'] = src if lines else ''
        pre = [l for l in lines if l[0] <= 0]
        post = [l for l in lines if l[0] > 0]
        r['n_prior_lines'] = len(pre) if lines else None
        if pre:
            d0, rg, d1, rp, dp = pre[-1]
            r['pre_reg'], r['pre_start'], r['pre_stop'], r['pre_resp'] = rg, d0, d1, rp
        if post:
            d0, rg, d1, rp, dp = post[0]
            r['post_reg'], r['post_start'], r['post_stop'] = rg, d0, d1
            r['post_resp'], r['post_prog'] = rp, dp
        f = fup[fup.study_id == r['study_id']]
        if len(f):
            x = f.iloc[0]
            r['vital'] = x.vital_status.strip()
            r['day_death'] = num(x.day_death)
            r['cod'] = x.cause_of_death
        m = mark[(mark.study_id == r['study_id']) & (mark.test == 'CA 19-9')].copy()
        if len(m):
            m['d'] = pd.to_numeric(m.day, errors='coerce')
            m['val'] = pd.to_numeric(m.result, errors='coerce')
            m = m.dropna(subset=['d', 'val']).sort_values('d')
            r['ca199_n'] = len(m)
            pre_m = m[m.d <= 0]
            if len(pre_m):
                r['ca199_pre'] = round(pre_m.iloc[-1].val, 1)
                r['ca199_pre_day'] = int(pre_m.iloc[-1].d)
            post_m = m[m.d > 0]
            if len(post_m):
                nad = post_m.loc[post_m.val.idxmin()]
                r['ca199_nadir'] = round(nad.val, 1)
                r['ca199_nadir_day'] = int(nad.d)
                r['ca199_last'] = round(post_m.iloc[-1].val, 1)
                r['ca199_last_day'] = int(post_m.iloc[-1].d)
    # ---- assay-setup fields derivable from what is present
    tp_order = ['24h', '48h', '72h']
    tps = {}
    for (samp, tp), _ in wells.groupby(['sample', 'timepoint']):
        tps.setdefault(samp, set()).add(tp)
    for samp, r in recs.items():
        r['in_cohort'] = 'Yes'
        got = sorted(tps.get(samp, ()), key=lambda t: tp_order.index(t) if t in tp_order else 9)
        if got:
            r['timepoints'] = '; '.join(t.replace('h', ' h') for t in got)
        if any(k.startswith('ihc_') and not k.startswith('ihc_n_') for k in r):
            r['has_ihc'] = 'Yes'
        if any(k.startswith('nd_') for k in r):
            r['has_nd'] = 'Yes'
        r['has_erk'] = 'Yes' if any(k.startswith('erk_') for k in r) else ''
        for d in DRUGS:
            keys = ['mts_%s_24h', 'mts_%s_48h', 'mts_%s_72h', 'mlmts_%s', 'nd_%s',
                    'ihc_cc3_%s', 'ihc_ki67_%s', 'ihc_cd45_%s', 'ihc_epcam_%s', 'dose_%s']
            if any(r.get(k % d) not in (None, '') for k in keys):
                r['tested_%s' % d] = 'Yes'
        blob = (r['assays'] + ' ' + r['assay_details'])
        if 'H&E' in blob:
            r['has_he'] = 'Yes'
        m = re.search(r'Viability, (\d+) days', blob)
        if m:
            r['days_culture'] = int(m.group(1))
    return [recs[s] for s, *_ in SPINE]

if __name__ == '__main__':
    rows = build()
    for r in rows:
        print(r['n'], r['sample'], r['study_id'], 'day', r['study_day'], '|', len(r), 'fields')
