#!/usr/bin/env python3
"""Generate one clinical-course CSV per case patient from the chart abstraction.

Combines therapy courses, CA 19-9 / ctDNA results and imaging RECIST calls into the tidy
format make_clinical_course.py consumes. All days are relative to that patient's biopsy.
"""
import pandas as pd, numpy as np, re

TH = pd.read_csv("chart_therapy_deid.csv")
MK = pd.read_csv("chart_markers_deid.csv")
IM = pd.read_csv("chart_imaging_deid.csv")

# ULN 35 U/mL is the value the lab's own PT015_Clinical-Course figure draws; it is not in the
# abstraction, so it is carried here as an assumption and labelled as such in the figure.
ULN = 35.0

TRACK = [(r"folfirinox|folfiri|5-?fu|irinotecan|capecitabin|gem|cisplatin|nal-?iri", "Systemic chemotherapy"),
         (r"cobimetinib|trametinib",                                                 "MEK inhibitor"),
         (r"tovorafenib|vemurafenib|daraxonrasib|rmc",                               "RAF / RAS inhibitor"),
         (r"vaccine|leidos|trial|pembro",                                            "Trial / immunotherapy"),
         (r"sbrt|radiat|chemort",                                                    "Radiation")]

def track_of(reg):
    s = str(reg).lower()
    for pat, name in TRACK:
        if re.search(pat, s): return name
    return "Other"

def parse_result(v):
    """returns (value, censor) where censor is '', 'above' or 'below'"""
    s = str(v).strip().replace(",", "")
    c = "above" if s.startswith(">") else ("below" if s.startswith("<") else "")
    try:    return float(s.lstrip("<>")), c
    except Exception: return np.nan, c

CASES = {"PT-024": "KRAS G12D metastatic PDAC",
         "PT-015": "BRAF V487 metastatic PDAC",
         "PT-014": "KRAS G12D metastatic PDAC",
         "PT-025": "KRAS G12R metastatic PDAC"}

for sid, label in CASES.items():
    rows = [dict(kind="event", label="Index biopsy (slice culture)", day=0, etype="biopsy")]

    for _, r in TH[TH.study_id == sid].iterrows():
        reg = str(r.regimen)
        if re.search(r"resect|whipple", reg, re.I):
            rows.append(dict(kind="event", label=reg, day=r.day_start, etype="surgery")); continue
        if pd.isna(r.day_start): continue
        rows.append(dict(kind="therapy", track=track_of(reg), label=reg,
                         day=r.day_start, day_end=r.day_stop,
                         annotate=(r.best_response if pd.notna(r.best_response) else ""),
                         ongoing=1 if pd.isna(r.day_stop) else ""))

    mk = MK[MK.study_id == sid].sort_values("day")
    for _, r in mk.iterrows():
        v, cen = parse_result(r.result)
        if not np.isfinite(v): continue
        rows.append(dict(kind="marker", series=("CA 19-9 (U/mL)" if "19-9" in str(r.test)
                                                else str(r.test)),
                         day=r.day, value=v, censored=cen))
    if len(mk):
        rows.append(dict(kind="hline", series="CA 19-9 (U/mL)", value=ULN, label="ULN 35"))

    for _, r in IM[IM.study_id == sid].iterrows():
        if pd.isna(r.recist): continue
        rows.append(dict(kind="scan", day=r.day, label=str(r.recist).upper()))

    df = pd.DataFrame(rows)
    for c in ["kind","series","track","label","day","day_end","value","annotate","etype",
              "ongoing","censored"]:
        if c not in df: df[c] = ""
    out = f"course_{sid.replace('-','')}.csv"
    df[["kind","series","track","label","day","day_end","value","annotate","etype",
        "ongoing","censored"]].to_csv(out, index=False)
    n = lambda k: int((df.kind == k).sum())
    print(f"{out:<22} therapy {n('therapy')}  markers {n('marker')}  scans {n('scan')}  events {n('event')}")
