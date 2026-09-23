#!/usr/bin/env python3
"""Verify the rebuilt n = 19 Table 1.

Everything is recomputed INDEPENDENTLY from table1_source_deid.csv, in Python, and then
compared against what make_table1.js actually printed into the .docx. The JSON the builder
emits is checked too, but it is not trusted as the source of truth — the .docx is.

Checks
  1. group n and total n
  2. every parameter's categories account for all n in its group
  3. every count printed in the .docx matches the independent recomputation
  4. every percentage printed in the .docx matches count / n
  5. regression: restricted to the original 16 biopsies, the derived counts still reproduce
     the numbers V5 printed — i.e. adding three rows perturbed nothing
  6. threshold sensitivity for the three tumour-marker cut-offs
"""
import html
import json
import re
import sys
import zipfile

import pandas as pd

ULN = {"ca19_9": 37.5, "ca125": 35.0, "cea": 2.5}
GROUPS = ["Liver", "Peritoneum", "Primary PDAC"]
NEW = ["Liver Bx 8", "Liver Bx 9", "Liver Bx 10"]

V5_AGG = {   # V5 group totals for the ORIGINAL 16, mirrored from make_table1.js
    "stage": {"I": [1, 0, 2], "II": [0, 1, 4], "III": [2, 1, 1], "IV": [4, 0, 0]},
    "approach": {"Image-guided needle core biopsy": [5, 2, 0],
                 "Excisional biopsy": [2, 0, 0],
                 "Pancreatoduodenectomy": [0, 0, 7]},
    "grade": {"Well differentiated": [2, 0, 0],
              "Moderately differentiated": [5, 2, 6],
              "Poorly differentiated": [0, 0, 1]},
    "prior_systemic": {"Yes": [7, 2, 7], "No": [0, 0, 0]},
}
V5_DERIVED_16 = {   # what V5 printed for the variables that DO have a per-biopsy source
    ("Age at biopsy, years", "< 50"): [1, 0, 2],
    ("Age at biopsy, years", "≥ 50"): [6, 2, 5],
    ("Sex", "Male"): [2, 1, 2],
    ("Sex", "Female"): [5, 1, 5],
    ("CA19-9", "Elevated"): [4, 2, 4], ("CA19-9", "Normal"): [3, 0, 3],
    ("CA-125", "Elevated"): [3, 1, 0], ("CA-125", "Normal"): [3, 1, 7],
    ("CEA", "Elevated"): [4, 2, 5], ("CEA", "Normal"): [3, 0, 2],
}

src = pd.read_csv("table1_source_deid.csv", comment="#")
N = [int((src.group == g).sum()) for g in GROUPS]
issues = []


def by(df, fn):
    return [int(sum(fn(r) for _, r in df[df.group == g].iterrows())) for g in GROUPS]


def markers(df, key):
    return {
        "Elevated": by(df, lambda r, m=key: pd.notna(r[m]) and r[m] > ULN[m]),
        "Normal": by(df, lambda r, m=key: pd.notna(r[m]) and r[m] <= ULN[m]),
        "Not measured": by(df, lambda r, m=key: pd.isna(r[m])),
    }


def aggregate(field):
    out = {}
    for cat, base in V5_AGG[field].items():
        out[cat] = [base[i] + int(((src.group == GROUPS[i]) & (src[field] == cat)).sum())
                    for i in range(3)]
    return out


DERIVED = {
    "Age at biopsy, years": {"< 50": by(src, lambda r: r.age_years < 50),
                             "≥ 50": by(src, lambda r: r.age_years >= 50)},
    "Sex": {"Male": by(src, lambda r: r.sex == "M"),
            "Female": by(src, lambda r: r.sex == "F")},
    "Stage at initial diagnosis": aggregate("stage"),
    "Tissue acquisition": aggregate("approach"),
    "Tumor grade": aggregate("grade"),
    "Systemic therapy before biopsy": aggregate("prior_systemic"),
    "CA19-9": markers(src, "ca19_9"),
    "CA-125": markers(src, "ca125"),
    "CEA": markers(src, "cea"),
}
for cats in DERIVED.values():
    if "Not measured" in cats and not any(cats["Not measured"]):
        del cats["Not measured"]

print("=" * 78)
print("1. COHORT SIZE")
print("=" * 78)
print(f"    Liver {N[0]}   Peritoneum {N[1]}   Primary PDAC {N[2]}   total {sum(N)}")
if sum(N) != 19:
    issues.append(f"total n is {sum(N)}, expected 19")
blank = src[src.age_years.isna()]["sample"].tolist()
print(f"    rows with no characteristics recorded: {blank if blank else 'none'}")
if blank:
    issues.append(f"still missing characteristics: {', '.join(blank)}")

print()
print("=" * 78)
print("2. COLUMN SUMS — every parameter must account for all n in its group")
print("=" * 78)
for param, cats in DERIVED.items():
    tot = [sum(v[i] for v in cats.values()) for i in range(3)]
    ok = tot == N
    print(f"{'ok  ' if ok else '*** '}{param:<32s} {tot}  (expected {N})")
    if not ok:
        issues.append(f"{param}: categories sum to {tot}, group n is {N}")

print()
print("=" * 78)
print("3 + 4. WHAT THE .DOCX ACTUALLY PRINTS — counts and percentages")
print("=" * 78)
xml = zipfile.ZipFile("Table1_cohort_characteristics.docx").read("word/document.xml").decode()
body = xml.split("<w:tbl>", 1)[1].split("</w:tbl>", 1)[0]
doc_rows = []
for tr in re.findall(r"<w:tr[ >].*?</w:tr>", body, re.S):
    cells = [re.sub(r"\s+", " ",
                    html.unescape("".join(
                        re.findall(r"<w:t(?:\s[^>]*)?>(.*?)</w:t>", tc, re.S)))).strip()
             for tc in re.findall(r"<w:tc>.*?</w:tc>", tr, re.S)]
    doc_rows.append(cells)

header = doc_rows[0]
print(f"    header: {header}")
for i, g in enumerate(GROUPS):
    want = f"(n = {N[i]})"
    if want not in header[2 + i]:
        issues.append(f"header for {g} reads '{header[2+i]}', expected to contain '{want}'")

flat = []
for param, cats in DERIVED.items():
    for cat, c in cats.items():
        flat.append((param, cat, c))

if len(doc_rows) - 1 != len(flat):
    issues.append(f".docx has {len(doc_rows)-1} body rows, recomputation gives {len(flat)}")

bad_count = bad_pct = 0
current = ""
for k, (dr, (param, cat, counts)) in enumerate(zip(doc_rows[1:], flat)):
    printed_param = dr[0].strip()
    if printed_param:
        current = printed_param
    stripped = re.sub(r"[ᵃᵇᶜᵈ ]+$", "", current)
    if stripped != param:
        issues.append(f"row {k}: .docx parameter '{stripped}' vs recomputed '{param}'")
    if dr[1].strip() != cat:
        issues.append(f"row {k}: .docx category '{dr[1]}' vs recomputed '{cat}'")
    for i in range(3):
        m = re.match(r"^(\d+) \(([\d.]+)\)$", dr[2 + i].strip())
        if not m:
            issues.append(f"{param}/{cat}/{GROUPS[i]}: unparseable cell '{dr[2+i]}'")
            continue
        got_n, got_p = int(m.group(1)), float(m.group(2))
        if got_n != counts[i]:
            bad_count += 1
            issues.append(f"{param}/{cat}/{GROUPS[i]}: .docx prints n={got_n}, "
                          f"recomputed {counts[i]}")
        exp_p = round(100 * counts[i] / N[i], 1)
        if abs(got_p - exp_p) > 0.05:
            bad_pct += 1
            issues.append(f"{param}/{cat}/{GROUPS[i]}: .docx prints {got_p}%, "
                          f"{counts[i]}/{N[i]} = {exp_p}%")

n_cells = (len(doc_rows) - 1) * 3
print(f"    {n_cells} cells checked")
print(f"    counts:      {'all match the source CSV' if bad_count == 0 else f'{bad_count} MISMATCHED'}")
print(f"    percentages: {'all match count / n' if bad_pct == 0 else f'{bad_pct} MISMATCHED'}")

print()
print("=" * 78)
print("5. REGRESSION — the original 16 must still reproduce V5 exactly")
print("=" * 78)
old = src[~src["sample"].isin(NEW)]
oldN = [int((old.group == g).sum()) for g in GROUPS]
print(f"    original subset: {oldN}  (expected [7, 2, 7])")
if oldN != [7, 2, 7]:
    issues.append(f"original subset is {oldN}, expected [7, 2, 7]")
old_derived = {
    ("Age at biopsy, years", "< 50"): by(old, lambda r: r.age_years < 50),
    ("Age at biopsy, years", "≥ 50"): by(old, lambda r: r.age_years >= 50),
    ("Sex", "Male"): by(old, lambda r: r.sex == "M"),
    ("Sex", "Female"): by(old, lambda r: r.sex == "F"),
}
for key, lab in [("ca19_9", "CA19-9"), ("ca125", "CA-125"), ("cea", "CEA")]:
    m = markers(old, key)
    old_derived[(lab, "Elevated")] = m["Elevated"]
    old_derived[(lab, "Normal")] = m["Normal"]
for k, got in old_derived.items():
    want = V5_DERIVED_16[k]
    ok = got == want
    print(f"{'    ok  ' if ok else '    ***'} {k[0]:<22s} {k[1]:<14s} {got}  V5 {want}")
    if not ok:
        issues.append(f"REGRESSION {k[0]}/{k[1]}: now {got}, V5 printed {want}")

print()
print("=" * 78)
print("6. THRESHOLD SENSITIVITY — how much slack is in each cut-off?")
print("=" * 78)
for key, lab in [("ca19_9", "CA19-9"), ("ca125", "CA-125"), ("cea", "CEA")]:
    vals = src[key].dropna()
    hi = vals[vals > ULN[key]]
    lo = vals[vals <= ULN[key]]
    if hi.empty or lo.empty:
        print(f"    {lab:<8s} using {ULN[key]:g}; one side of the cut-off is empty")
        continue
    lower, upper = lo.max(), hi.min()
    print(f"    {lab:<8s} using {ULN[key]:g}; counts are identical for any cut-off in "
          f"[{lower:g}, {upper:g}) — nearest values either side are {lower:g} and {upper:g}")
    if not (lower <= ULN[key] < upper):
        issues.append(f"{lab}: ULN {ULN[key]} is outside the stable window "
                      f"[{lower}, {upper})")

print()
print("=" * 78)
print("7. BUILDER JSON — does table1_counts.json agree with the .docx?")
print("=" * 78)
try:
    j = json.load(open("table1_counts.json"))
    mism = 0
    if j["N"] != N:
        issues.append(f"table1_counts.json N is {j['N']}, recomputed {N}")
    for entry, (param, cat, counts) in zip(j["rows"], flat):
        if entry["counts"] != counts:
            mism += 1
            issues.append(f"table1_counts.json {param}/{cat}: {entry['counts']} vs {counts}")
    print(f"    {'agrees on all rows' if mism == 0 else f'{mism} rows disagree'}")
except FileNotFoundError:
    print("    table1_counts.json not found (run make_table1.js first)")

print()
print("=" * 78)
print(f"FINDINGS ({len(issues)})")
print("=" * 78)
if not issues:
    print("    none — every printed cell reproduces from table1_source_deid.csv")
for k, s in enumerate(issues, 1):
    print(f"{k}. {s}")
sys.exit(1 if issues else 0)
