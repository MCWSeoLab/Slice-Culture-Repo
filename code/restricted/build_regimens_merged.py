#!/usr/bin/env python3
"""Single regimen table for the figures: chart abstraction where it exists, workbook elsewhere.

Emits regimens_merged.csv in the workbook's column schema so the figure scripts need only
change which file they read.
"""
import pandas as pd, numpy as np

ch = pd.read_csv("chart_therapy_deid.csv")
chart_pts = set(ch.study_id)

out = pd.DataFrame({
    "study_id": ch.study_id,
    "phase": np.where(ch.day_start >= 0, "post", "pre"),
    "src_sheet": "chart",
    "regimen": ch.regimen,
    "day_start_rel_biopsy": ch.day_start,
    "day_end_rel_biopsy": ch.day_stop,
    "cycles_completed": ch.cycles,
    "response": ch.best_response,
    "day_progression": ch.day_progression,
    "reason_stopped": ch.reason_stopped,
})

old = pd.read_csv("regimens_deid.csv")
# keep workbook rows for non-chart patients, PLUS the surgical / radiation events for chart
# patients - the abstraction sheet only captured systemic lines, so Whipple and SBRT live
# nowhere else and would otherwise vanish from the timelines
EVENT = r"whipple|sbrt|resect|chemort"
is_event = old.regimen.astype(str).str.contains(EVENT, case=False, na=False)
already  = set(zip(ch.study_id, ch.regimen.astype(str).str.lower()))
dup = pd.Series([(s, str(r).lower()) in already
                 for s, r in zip(old.study_id, old.regimen)], index=old.index)
old = old[(~old.study_id.isin(chart_pts)) | (is_event & ~dup)].copy()
old["day_progression"] = np.nan
old["reason_stopped"] = np.nan
keep = ["study_id","phase","src_sheet","regimen","day_start_rel_biopsy","day_end_rel_biopsy",
        "cycles_completed","response","day_progression","reason_stopped"]

m = pd.concat([out[keep], old[keep]], ignore_index=True)
m.to_csv("regimens_merged.csv", index=False)
print(f"regimens_merged.csv: {len(out)} chart rows for {len(chart_pts)} patients "
      f"+ {len(old)} workbook rows for {old.study_id.nunique()} patients")
for s in sorted(chart_pts):
    d = out[out.study_id == s]
    print(f"   {s}: {len(d)} rows, "
          f"{d.day_end_rel_biopsy.notna().sum()} stop dates, {d.response.notna().sum()} responses")
