#!/usr/bin/env python3
"""Extract well-level MTS readings from the two lab assay workbooks.

Replaces an earlier version that sniffed headers and got two sheets wrong. Every cohort sheet
now has an EXPLICIT column spec, written after reading that sheet by eye. Sniffing was the
defect: it inferred timepoints from nearby header text, mispaired controls with treated wells
on L.Bx 5 and L.Bx 6, and double-counted the layout blocks on the right of L.Bx 10.

Normalisation is done here rather than read from the sheets, because the sheets disagree about
what "Normalised" means — per well on some, per drug-group on others. Everywhere:

    normalised(well)  = raw(well) - mean(blank wells at that timepoint)
    pct_of_control    = 100 * normalised(well) / mean(normalised DMSO wells at that timepoint)

which is the SAP 4.1 definition. Outputs are de-identified: manuscript sample labels and
timepoints only, never a calendar date or a lab TB number.

Four cohort sheets are excluded BY DESIGN, not by failure: they are day-series feasibility
runs with no treated-versus-DMSO comparison, so they cannot contribute to either derivation.
Liver Bx 1 has no sheet in either workbook. Lab sheet "L.Bx 8" is a separate biopsy that was
not cultured and is not in the cohort; it is why the lab numbering runs one ahead from
Liver Bx 8 onward.
"""
import sys

import numpy as np
import openpyxl
import pandas as pd

LIVER = ("/root/.claude/uploads/5ca89d3c-111f-578e-bab9-226ffe80bccc/"
         "272a5f6c-Liver_Bx__Diaphgram_Bx__MTS_assay.xlsx")
PDAC = ("/root/.claude/uploads/5ca89d3c-111f-578e-bab9-226ffe80bccc/"
        "8afd034f-PDAC_Bx__Peritoneum_Bx__MTS_assay.xlsx")

def _sheet(book, name):
    """Return the worksheet whose title is `name`, or the single title beginning
    with it. Some sheet titles in the source workbooks carry a trailing collection
    date; matching on the stable prefix keeps that date out of this file."""
    if name in book.sheetnames:
        return book[name]
    matches = [t for t in book.sheetnames if t.startswith(name)]
    if len(matches) != 1:
        raise KeyError(f"{name!r} matched {len(matches)} sheets; expected exactly one")
    return book[matches[0]]

# ---------------------------------------------------------------------------- specs
# kind "raw"    : per-well RAW absorbance at the listed columns; blanks present on the sheet
# kind "norm"   : per-well values ALREADY blank-corrected; no blank subtraction to do
# kind "matrix" : rows x drug-columns, two stacked timepoint blocks, already blank-corrected
# kind "none"   : day-series feasibility run, no drug panel — excluded by design
SPEC = {
    "Liver Bx 2":      dict(book=LIVER, sheet="L.Bx 2", kind="none",
                            why="day-series feasibility run, no treated-vs-DMSO comparison"),
    "Liver Bx 3":      dict(book=LIVER, sheet="L.Bx 3 STS IC50", kind="none",
                            why="staurosporine IC50 dilution series, not a drug panel"),
    "Liver Bx 4":      dict(book=LIVER, sheet="L.Bx-4 Folfirinox", kind="raw",
                            drug=1, conc=2, well=3, tp={"24h": 4, "48h": 7}),
    "Liver Bx 5":      dict(book=LIVER, sheet="L.Bx 5 All Inhibitors", kind="raw",
                            drug=1, conc=2, well=3, tp={"24h": 4}),
    "Liver Bx 6":      dict(book=LIVER, sheet="L.Bx 6", kind="raw",
                            drug=1, conc=2, well=3, tp={"24h": 4, "48h": 7}),
    "Liver Bx 7":      dict(book=LIVER, sheet="L.Bx 7", kind="raw",
                            drug=1, conc=2, well=3, tp={"24h": 4, "48h": 6}),
    "Liver Bx 8":      dict(book=LIVER, sheet="L.Bx 9", kind="raw",
                            drug=1, conc=2, well=3, tp={"24h": 4, "48h": 6}),
    "Liver Bx 9":      dict(book=LIVER, sheet="L.Bx 10", kind="raw",
                            drug=1, conc=2, well=3, tp={"24h": 4, "48h": 6}),
    "Liver Bx 10":     dict(book=LIVER, sheet="L.Bx 11", kind="matrix",
                            drug_cols={5: "FOLFIRINOX", 6: "Gem/nab-Pac", 7: "RMC-7977",
                                       8: "Trametinib", 9: "DMSO"},
                            row_col=4, block_col=3),
    "Peritoneum Bx 1": dict(book=PDAC, sheet="Met Peritonium Bx-1", kind="none",
                            why="day-series feasibility run, no treated-vs-DMSO comparison"),
    "Peritoneum Bx 2": dict(book=PDAC, sheet="Met Peritonium Bx-2", kind="raw",
                            drug=1, conc=2, well=3, tp={"24h": 4, "48h": 7}),
    "PDAC Bx 1":       dict(book=PDAC, sheet="PDAC 1", kind="none",
                            why="day-series feasibility run, no treated-vs-DMSO comparison"),
    "PDAC Bx 2":       dict(book=PDAC, sheet="PDAC 2", kind="none",
                            why="day-series feasibility run, no treated-vs-DMSO comparison"),
    "PDAC Bx 3":       dict(book=PDAC, sheet="PDAC 3 Folfirnox 1", kind="raw",
                            drug=1, conc=2, well=3, tp={"24h": 4, "48h": 7}),
    "PDAC Bx 4":       dict(book=PDAC, sheet="PDAC 4 Folfirnox 2", kind="raw",
                            drug=1, conc=2, well=3, tp={"24h": 4, "48h": 7}),
    "PDAC Bx 5":       dict(book=PDAC, sheet="PDAC 5 All Inhibitors", kind="raw",
                            drug=1, conc=2, well=3, tp={"24h": 4, "48h": 7}),
    "PDAC Bx 6":       dict(book=PDAC, sheet="PDAC-6", kind="raw",
                            drug=1, conc=2, well=3, tp={"24h": 4, "48h": 7}),
    "PDAC Bx 7":       dict(book=PDAC, sheet="PDAC 7", kind="norm",
                            drug=1, conc=2, well=None,
                            tp={"24h": 3, "48h": 4, "72h": 5}),
}
ORDER = ([f"Liver Bx {i}" for i in range(2, 11)] + ["Peritoneum Bx 1", "Peritoneum Bx 2"]
         + [f"PDAC Bx {i}" for i in range(1, 8)])

books = {LIVER: openpyxl.load_workbook(LIVER, data_only=True),
         PDAC: openpyxl.load_workbook(PDAC, data_only=True)}


def canon(s):
    t = "".join(ch for ch in str(s).lower() if ch.isalnum())
    if not t: return None
    if "blank" in t: return "BLANK"
    if "dmso" in t or t in ("control", "untreated"): return "DMSO"
    if "folfirinox" in t or "folfirnox" in t: return "FOLFIRINOX"
    if "gem" in t: return "Gem/nab-Pac"
    if "rmc" in t or "rasinhib" in t: return "RMC-7977"
    if "trametinib" in t or t == "mek" or "mekinh" in t: return "Trametinib"
    if "vemuraf" in t or "rafinh" in t: return "Vemurafenib"
    if "staurosporine" in t or t == "sts": return "Staurosporine"
    if "leuco" in t or t == "5fu" or "sn38" in t or "irinotecan" in t or "oxalipl" in t:
        return "FOLFIRINOX component"
    return None


def num(v):
    return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None


rows, notes = [], []
for sample in ORDER:
    sp = SPEC[sample]
    if sp["kind"] == "none":
        notes.append((sample, sp["sheet"], sp["why"]))
        continue
    grid = [[c.value for c in r] for r in _sheet(books[sp["book"]], sp["sheet"]).iter_rows()]

    if sp["kind"] == "matrix":
        block = None
        for r in grid:
            lab = str(r[sp["block_col"]] or "")
            if "24 Hour" in lab: block = "24h"
            elif "48 Hour" in lab: block = "48h"
            if block is None: continue
            rid = str(r[sp["row_col"]] or "").strip()
            if rid not in ("A", "B", "C", "D"): continue
            for col, drug in sp["drug_cols"].items():
                v = num(r[col]) if col < len(r) else None
                if v is not None:
                    rows.append(dict(sample=sample, timepoint=block, drug=drug,
                                     concentration="", well=rid, normalised=v))
        continue

    # "raw" and "norm" share the row walk; only blank handling differs
    recs = []
    for r in grid:
        d = sp["drug"]
        drug_raw = r[d] if d < len(r) else None
        well = ""
        if sp.get("well") is not None and sp["well"] < len(r) and r[sp["well"]] is not None:
            well = str(r[sp["well"]]).strip()
        drug = canon(drug_raw)
        # the blank can be labelled in the drug column or in the well column
        if drug is None and well.lower() == "blank":
            drug = "BLANK"
        if drug is None or drug == "FOLFIRINOX component":
            continue
        conc = ""
        if sp.get("conc") is not None and sp["conc"] < len(r) and r[sp["conc"]] is not None:
            conc = str(r[sp["conc"]]).strip()
        for tp, col in sp["tp"].items():
            v = num(r[col]) if col < len(r) else None
            if v is not None:
                recs.append(dict(sample=sample, timepoint=tp, drug=drug, concentration=conc,
                                 well=well, value=v))

    df = pd.DataFrame(recs)
    if df.empty:
        notes.append((sample, sp["sheet"], "SPEC PARSED NOTHING — check the column indices"))
        continue
    for tp, g in df.groupby("timepoint"):
        blanks = g.loc[g.drug == "BLANK", "value"]
        if sp["kind"] == "raw":
            if blanks.empty:
                notes.append((sample, sp["sheet"], f"no blank well at {tp} — skipped"))
                continue
            off = float(blanks.mean())
        else:
            off = 0.0
        for _, r in g[g.drug != "BLANK"].iterrows():
            rows.append(dict(sample=sample, timepoint=tp, drug=r.drug,
                             concentration=r.concentration, well=r.well,
                             normalised=r.value - off))

w = pd.DataFrame(rows)

print("=" * 78)
print("COVERAGE")
print("=" * 78)
print(f"cohort sheets with a drug panel : {sum(1 for s in ORDER if SPEC[s]['kind'] != 'none')}")
print(f"parsed                          : {w['sample'].nunique()}")
print(f"well readings                   : {len(w)}")
print("\nexcluded by design (no treated-vs-DMSO comparison exists on these sheets):")
for s, sh, why in notes:
    print(f"    {s:<16s} {sh:<22s} {why}")
print("    Liver Bx 1       —                      no sheet in either workbook")

print()
print("=" * 78)
print("AUDIT — DMSO wells per experiment, for eyeballing against the sheets")
print("=" * 78)
for (s, tp), g in w[w.drug == "DMSO"].groupby(["sample", "timepoint"]):
    vals = ", ".join(f"{v:.4f}" for v in g.normalised)
    print(f"    {s:<16s} {tp:<5s} n={len(g)}  {vals}")

# ------------------------------------------------------------------ derivation B
cv = []
for (s, tp), g in w[w.drug == "DMSO"].groupby(["sample", "timepoint"]):
    v = g.normalised.dropna()
    cv.append(dict(sample=s, timepoint=tp, n_control_wells=len(v), mean=v.mean(),
                   sd=v.std(ddof=1) if len(v) > 1 else np.nan,
                   dmso_cv=(v.std(ddof=1) / v.mean()) if len(v) > 1 and v.mean() else np.nan))
cv = pd.DataFrame(cv)
cv["qc"] = np.where(cv.dmso_cv.isna(), "UNKNOWN (single control well)",
                    np.where(cv.dmso_cv > 0.20, "FAIL (CV > 20%)", "PASS"))

print()
print("=" * 78)
print("B. DMSO REPLICATE CV")
print("=" * 78)
print(cv.sort_values("dmso_cv").to_string(index=False, float_format=lambda x: f"{x:.4f}"))
p = cv[cv.qc == "PASS"]
print(f"\nexperiments with >= 2 DMSO wells : {int((cv.n_control_wells >= 2).sum())}")
print(f"passing the 20% QC gate          : {len(p)}  "
      f"(from {p['sample'].nunique()} biopsies)")
print(f"largest control-well count       : {int(cv.n_control_wells.max())}")
if len(p):
    MED = float(p.dmso_cv.median())
    SD = 100 * MED * np.sqrt(2)
    print(f"\nmedian CV                        : {MED*100:.2f}%")
    print(f"  -> SD of a single-well ratio   : {SD:.2f} pp")
    print(f"  -> 2 SD                        : {2*SD:.1f}% inhibition")
    print(f"  published value                : 10.8% -> 15.3 pp -> 30%")

# ------------------------------------------------------------------ derivation A
ctrl = (w[w.drug == "DMSO"].groupby(["sample", "timepoint"]).normalised.mean()
        .rename("control_mean"))
t = w[w.drug != "DMSO"].join(ctrl, on=["sample", "timepoint"])
t = t[t.control_mean.notna() & (t.control_mean > 0)].copy()
t["pct_of_control"] = 100 * t.normalised / t.control_mean
fail = set(map(tuple, cv.loc[cv.qc.str.startswith("FAIL"), ["sample", "timepoint"]].values))
t["qc_fail"] = [(a, b) in fail for a, b in zip(t["sample"], t.timepoint)]

rng = np.random.default_rng(20260818)
def hn(d):
    d = np.asarray(d, float)
    return float(np.sqrt(np.mean(d ** 2))) if len(d) else float("nan")

print()
print("=" * 78)
print("A. HALF-NORMAL FIT TO TREATED WELLS ABOVE 100% OF CONTROL")
print("=" * 78)
for lab, sub in [("all treated wells", t), ("QC-passing experiments only", t[~t.qc_fail])]:
    d = (sub.loc[sub.pct_of_control > 100, "pct_of_control"] - 100).to_numpy(float)
    s = hn(d)
    if len(d) > 1:
        bs = [hn(d[rng.integers(0, len(d), len(d))]) for _ in range(4000)]
        lo, hi = np.percentile(bs, 2.5), np.percentile(bs, 97.5)
    else:
        lo = hi = float("nan")
    print(f"\n  {lab}: {len(sub)} wells, {len(d)} above 100%")
    if len(d):
        print(f"    sigma {s:.2f} pp (bootstrap {lo:.1f}-{hi:.1f})  ->  2 sigma = {2*s:.1f}%")
        print(f"    excesses: {np.round(np.sort(d), 1)}")
print("\n  published value: sigma = 16.2 pp -> 32%")

w.to_csv("mts_wells_deid.csv", index=False)
cv.to_csv("mts_dmso_cv_wells_deid.csv", index=False)
t.to_csv("mts_treated_pct_deid.csv", index=False)
print("\nwrote mts_wells_deid.csv, mts_dmso_cv_wells_deid.csv, mts_treated_pct_deid.csv")
