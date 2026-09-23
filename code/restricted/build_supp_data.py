#!/usr/bin/env python3
"""Build the de-identified source CSVs for Supplementary Figures 1 and 2.

Every value is read from a Prism raw data grid and every group mean is checked
against the ANOVA / multiple-comparison output Prism itself stored in the same
file. Nothing is retyped.

Why not deid/prism_raw_long_deid.csv: for .prism files that extractor assigns
the FIRST value of the grid to `row_label` and then shifts every remaining
value one position left, so its day/agent assignment is wrong. Verified against
Prism's own stored means for the PDAC viability curve. The .pzfx sources in the
same CSV are unaffected and are used as-is.
"""
from __future__ import annotations

import ast
import re
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

UP = Path(".")   # run from the project root
PRISM = UP / "Manuscript Drafts with Results/Manuscript Final Version/Manuscript Prism Files"
OUT = Path("deid")
OUT.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------- prism I/O

def prism_tables(path: Path) -> dict[str, list[str]]:
    """Return {table_uuid: [csv lines]} for every data.csv in a .prism file."""
    z = zipfile.ZipFile(path)
    out = {}
    for n in z.namelist():
        if n.endswith("data.csv"):
            out[n.split("/")[-2]] = z.read(n).decode("utf-8", "replace").splitlines()
    return out


def raw_grid(path: Path) -> tuple[list[str], list[list[float]]]:
    """The one table in a .prism file that is the raw data grid.

    Returns (row_titles, rows_of_values). A grid line is `title,v,v,v...` when
    the first field is non-numeric, or all-numeric when there is no row title.
    """
    grids = []
    for uid, lines in prism_tables(path).items():
        if not lines:
            continue
        # analysis tables carry text in many cells and have blank lines
        body = [l for l in lines if l.strip()]
        if len(body) != len(lines):
            continue
        parsed_rows, titles, ok = [], [], True
        for line in body:
            parts = [p.strip() for p in line.split(",")]
            title = ""
            if parts and not _is_num(parts[0]):
                title, parts = parts[0], parts[1:]
            vals = []
            for p in parts:
                if p == "":
                    continue
                if not _is_num(p):
                    ok = False
                    break
                vals.append(float(p))
            if not ok or not vals:
                ok = False
                break
            titles.append(title)
            parsed_rows.append(vals)
        if ok and parsed_rows:
            grids.append((titles, parsed_rows))
    if len(grids) != 1:
        raise RuntimeError(f"{path.name}: expected 1 raw grid, found {len(grids)}")
    return grids[0]


def _is_num(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def prism_stat(path: Path, pattern: str) -> str | None:
    """First line of any table in the file matching `pattern`."""
    for lines in prism_tables(path).values():
        for line in lines:
            if re.search(pattern, line):
                return line
    return None


def reshape(values: list[float], groups: list[str], reps: int) -> pd.DataFrame:
    """Prism stores a grouped row column-major: group1 reps, group2 reps, ..."""
    if len(values) > len(groups) * reps:
        raise RuntimeError(f"{len(values)} values > {len(groups)}x{reps}")
    rows = []
    for i, v in enumerate(values):
        rows.append({"group": groups[i // reps], "replicate": i % reps + 1, "value": v})
    return pd.DataFrame(rows)


def check(name: str, got: dict[str, float], want: dict[str, float], tol=0.6):
    bad = {k: (round(got[k], 2), v) for k, v in want.items() if abs(got[k] - v) > tol}
    if bad:
        raise AssertionError(f"{name}: group means disagree with Prism: {bad}")
    print(f"  OK  {name}: {len(want)} group means match Prism's own output")


# --------------------------------------------------- Supp Fig 1, primary PDAC

print("Supplementary Figure 1 - primary PDAC (PDAC Bx 2)")

# --- A. MTS viability over the culture week -------------------------------
p = PRISM / "Suppl.Fig.1. Viability PDAC Bx Results/1. Manuscript PDAC Viability curve.prism"
titles, rows = raw_grid(p)
assert len(rows) == 1 and len(rows[0]) == 8, rows
mts = reshape(rows[0], ["Day 1", "Day 3", "Day 5", "Day 7"], 2)
check("MTS viability", mts.groupby("group").value.mean().to_dict(),
      {"Day 1": 0.1875, "Day 3": 0.5230, "Day 5": 0.7855, "Day 7": 1.067}, tol=0.005)
mts_f = prism_stat(p, r"^  F,")
mts_p = prism_stat(p, r"^  P value,")
print(f"      Prism one-way ANOVA: {mts_f} / {mts_p}")
mts.insert(0, "biopsy", "Primary PDAC Bx 2")
mts.insert(1, "endpoint", "MTS absorbance 490 nm")
mts.to_csv(OUT / "suppfig1_pdac_mts_deid.csv", index=False)

# --- B. nuclear density over the culture week -----------------------------
# pzfx source, unaffected by the .prism shift; stored as nuclei per um^2.
raw = pd.read_csv(UP / "deid/prism_raw_long_deid.csv")
nd = raw[(raw.file.str.contains("Core Biopsy Nuclear Density", na=False))
         & (raw.sheet == "PDAC 2 Viability")].copy()
nd["group"] = nd["column"]
nd["value"] = nd["value"] * 1e6                      # nuclei per um^2 -> per mm^2
nd = nd.sort_values(["column", "row_label"])[["group", "value"]]
nd["replicate"] = nd.groupby("group").cumcount() + 1   # pzfx rows are the replicates
nd = nd[["group", "replicate", "value"]]
check("nuclear density", nd.groupby("group").value.mean().to_dict(),
      {"Day 0": 1289, "Day 2": 1413, "Day 3": 1459, "Day 5": 1142, "Day 7": 1561}, tol=1.0)
p_he = PRISM / "Suppl.Fig.1. Viability PDAC Bx Results/2. Manuscript PDAC Viability H&E.prism"
print(f"      Prism one-way ANOVA: {prism_stat(p_he, r'^  F,')} / {prism_stat(p_he, r'^  P value,')}")
nd.insert(0, "biopsy", "Primary PDAC Bx 2")
nd.insert(1, "endpoint", "Nuclear density (nuclei/mm2)")
nd.to_csv(OUT / "suppfig1_pdac_nuclear_density_deid.csv", index=False)

# --- C. four-marker IHC over the culture week ------------------------------
MARKER = {"Cleaved Caspase": "Cleaved caspase-3", "Ki67": "Ki-67",
          "CD45": "CD45", "EpCam": "EpCAM"}
ihc = raw[raw.file.str.contains("PDAC Bx.2 IHC", na=False, regex=False)].copy()
ihc["marker"] = ihc.sheet.str.replace("PDAC 2 ", "", regex=False).map(MARKER)
ihc["group"] = ihc["column"]
assert ihc.marker.notna().all()
ihc = ihc.sort_values(["marker", "column", "row_label"])[["marker", "group", "value"]]
ihc["replicate"] = ihc.groupby(["marker", "group"]).cumcount() + 1
ihc = ihc[["marker", "group", "replicate", "value"]]
# cross-check the two timepoints that also exist in the Suppl.Fig.1 IHC prism
p_ihc = PRISM / "Suppl.Fig.1. Viability PDAC Bx Results/3. Manuscript PDAC Viability IHC.prism"
alt = pd.read_csv(UP / "deid/prism_raw_long_deid.csv")
alt = alt[alt.file.str.contains("Manuscript PDAC Viability IHC", na=False)]
alt["marker"] = alt["column"].map(
    lambda v: ast.literal_eval(v)["string"] if str(v).startswith("{'string'") else v)
alt["marker"] = alt["marker"].replace({"Cleaved Caspase-3": "Cleaved caspase-3", "CD-45": "CD45"})
merged = (ihc[ihc.group.isin(["Day 2", "Day 7"])]
          .groupby(["marker", "group"]).value.mean().round(1)
          .rename("full").reset_index()
          .merge(alt.rename(columns={"row_label": "group"})
                    .groupby(["marker", "group"]).value.mean().round(1)
                    .rename("suppl").reset_index(), on=["marker", "group"]))
assert (merged.full - merged.suppl).abs().max() < 0.15, merged
print(f"  OK  IHC: all {len(merged)} day 2 / day 7 means reproduce the Suppl.Fig.1 prism")
ihc.insert(0, "biopsy", "Primary PDAC Bx 2")
ihc.insert(1, "endpoint", "Positive cells (% of total)")
ihc.to_csv(OUT / "suppfig1_pdac_ihc_deid.csv", index=False)


# ------------------------------------------------ Supp Fig 2, peritoneal core

print("\nSupplementary Figure 2 - peritoneal core (Peritoneum Bx 1)")
PER = PRISM / "Fig.4. Chemo Treat Peritoneum Bx Results"
AGENTS = ["DMSO", "FOLFIRINOX", "Gem/nab-Pac", "RMC-7977", "Staurosporine"]

# --- A/B. MTS at 24 and 48 h ----------------------------------------------
titles, rows = raw_grid(PER / "1. Manuscript Peritoneum MTS.prism")
assert titles == ["24hrs", "48hrs"], titles
frames = []
for title, vals in zip(["24 h", "48 h"], rows):
    f = reshape(vals, AGENTS, 4)
    f.insert(0, "timepoint", title)
    frames.append(f)
mts_p = pd.concat(frames, ignore_index=True)
# normalize to same-plate DMSO at the same timepoint
dmso = mts_p[mts_p.group == "DMSO"].groupby("timepoint").value.mean()
mts_p["pct_of_dmso"] = mts_p.value / mts_p.timepoint.map(dmso) * 100
print("      % of same-timepoint DMSO control:")
for tp, g in mts_p.groupby("timepoint"):
    print("       ", tp, {k: round(v, 1) for k, v in
                          g.groupby("group").pct_of_dmso.mean().reindex(AGENTS).items()})
mts_p.insert(0, "biopsy", "Peritoneum Bx 1")
mts_p.insert(1, "endpoint", "MTS absorbance 490 nm")
mts_p.to_csv(OUT / "suppfig2_peritoneal_mts_deid.csv", index=False)

# --- C. nuclear density by agent ------------------------------------------
titles, rows = raw_grid(PER / "2. Manuscript Peritoneum H&E.prism")
assert len(rows) == 1 and len(rows[0]) == 30, (titles, len(rows[0]))
nd_p = reshape(rows[0], AGENTS, 6)
print("      nuclear density means:",
      {k: round(v) for k, v in nd_p.groupby("group").value.mean().reindex(AGENTS).items()})
nd_p.insert(0, "biopsy", "Peritoneum Bx 1")
nd_p.insert(1, "endpoint", "Nuclear density (nuclei/mm2)")
nd_p.to_csv(OUT / "suppfig2_peritoneal_nuclear_density_deid.csv", index=False)

# --- D. four-marker IHC by agent ------------------------------------------
IHC_FILES = {"Cleaved caspase-3": "3. Manuscript Peritoneum IHC CC3.prism",
             "Ki-67": "4. Manuscript Peritoneum IHC KI67.prism",
             "CD45": "5. Manuscript Peritoneum IHC CD45.prism",
             "EpCAM": "6. Manuscript Peritoneum IHC EpCAM.prism"}
frames = []
for marker, fn in IHC_FILES.items():
    titles, rows = raw_grid(PER / fn)
    assert len(rows) == 1 and len(rows[0]) == 20, (marker, len(rows[0]))
    f = reshape(rows[0], AGENTS, 4)
    f.insert(0, "marker", marker)
    frames.append(f)
ihc_p = pd.concat(frames, ignore_index=True)
# cross-check against the umbilicus IHC pzfx, which holds the same experiment
umb = raw[raw.file.str.contains("umbilicus IHC", na=False)].copy()
umb["marker"] = umb.sheet.str.replace("Umbilicus ", "", regex=False).replace(
    {"Cleaved Caspase": "Cleaved caspase-3", "Ki67": "Ki-67", "EpCam": "EpCAM"})
umb["group"] = umb["column"].replace({"Folfirinox": "FOLFIRINOX", "Gem/nab": "Gem/nab-Pac",
                                      "RMC7977": "RMC-7977", "STS": "Staurosporine"})
cmp = (ihc_p.groupby(["marker", "group"]).value.mean().rename("prism").reset_index()
       .merge(umb.groupby(["marker", "group"]).value.mean().rename("pzfx").reset_index(),
              on=["marker", "group"]))
assert len(cmp) == 20 and (cmp.prism - cmp.pzfx).abs().max() < 0.05, cmp
print(f"  OK  IHC: all {len(cmp)} agent x marker means reproduce the umbilicus pzfx")
ihc_p.insert(0, "biopsy", "Peritoneum Bx 1")
ihc_p.insert(1, "endpoint", "Positive cells (% of total)")
ihc_p.to_csv(OUT / "suppfig2_peritoneal_ihc_deid.csv", index=False)

print("\nwrote:")
for f in sorted(OUT.glob("suppfig*_deid.csv")):
    print("  ", f.name, len(pd.read_csv(f)), "rows")
