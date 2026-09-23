#!/usr/bin/env python3
"""Derive the released data set from the restricted de-identified tables.

The restricted tables live outside this repository. They are already free of the
18 HIPAA Safe Harbor identifiers, but they carry more clinical metadata than the
release is meant to expose. This script applies the release policy:

  clinical variables released : demographics, progression-free interval / time to
                                next therapy, best overall response, genomics
  clinical variables withheld : serial tumor markers, scan-level RECIST target
                                measurements, therapy logs, follow-up and vital
                                status, free-text notes and correlation columns

Laboratory variables (MTS absorbance, Prism exports, IHC scoring, nuclear
density, tumor content, QC) are released unchanged; they carry no clinical
metadata.

Usage:
    python build_release_dataset.py --source <restricted_deid_dir> --out <repo>/data
"""
from __future__ import annotations

import argparse
import io
import re
import sys
from pathlib import Path

import pandas as pd

# ---------------------------------------------------------------- release policy

COPY_AS_IS = [
    "mts_wells_deid.csv",
    "mts_treated_pct_deid.csv",
    "mts_dmso_cv_wells_deid.csv",
    "mts_normalized_final.csv",
    "prism_raw_long_deid.csv",
    "prism_mts_normalized_deid.csv",
    "prism_index_deid.csv",
    "ihc_tidy.csv",
    "erk_tidy.csv",
    "tumor_content_deid.csv",
    "tumour_content_deid.csv",
    "fig1_timecourse_deid.csv",
    "suppfig1_pdac_mts_deid.csv",
    "suppfig1_pdac_ihc_deid.csv",
    "suppfig1_pdac_nuclear_density_deid.csv",
    "suppfig2_peritoneal_mts_deid.csv",
    "suppfig2_peritoneal_ihc_deid.csv",
    "suppfig2_peritoneal_nuclear_density_deid.csv",
    "qc_dmso_cv.csv",
    "qc_by_patient.csv",
    "sensitivity_grid.csv",
    "forest_rho.csv",
    "pairs_final.csv",
    "pairs_recist.csv",
    "pairs_pfs.csv",
]

DROP_COLUMNS = {
    "table1_source_deid.csv": ["ca19_9", "ca125", "cea"],
    "table2_source_deid.csv": ["markers"],
    "table3_assaylog_source_deid.csv": [],
    "masterlist_linked_deid.csv": ["notes"],
    # Retained for the response-concordance analysis: the worksheet's own laboratory
    # call and the recorded clinical response. Free-text status and the derived
    # correlation columns are removed.
    "per_sheet_annotations_deid.csv": ["status", "pre_match", "post_match", "correlation"],
}

WITHHELD = {
    "chart_markers_deid.csv": "serial tumor marker results",
    "chart_imaging_deid.csv": "scan-level RECIST target measurements",
    "chart_therapy_deid.csv": "line-by-line therapy log with free-text reasons",
    "chart_followup_deid.csv": "follow-up, vital status and cause of death",
    "chart_erk_methods.csv": "assay requisition detail",
    "chart_ctdna_deid.csv": "serial circulating tumor DNA results",
    "regimens_deid.csv": "full regimen log",
    "regimens_merged.csv": "full regimen log",
    "patients_deid.csv": "per-patient free-text notes and correlation columns",
    "clinic_correlation_deid.csv": "free-text clinical correlation worksheet",
    "biopsy_exvivo_deid.csv": "extended sequencing report text beyond the reported driver",
    "exvivo_master_specimens_deid.csv": "program-wide cohort outside this report",
    "exvivo_master_drugtests_deid.csv": "program-wide cohort outside this report",
    "accrual_studyday_deid.csv": "program-wide accrual timeline outside this report",
}

COMMENTED = {
    "table1_source_deid.csv",
    "table2_source_deid.csv",
    "table3_assaylog_source_deid.csv",
}

GENOMIC_SOURCE = "table2_source_deid.csv"
GENOMIC_OUT = "genomics_deid.csv"

# Alterations are recorded in the worksheet as a gene symbol optionally followed by
# a protein-level variant or the word "wild-type". Entries are separated by commas
# or semicolons. The gene symbol is the first whitespace-delimited token; whatever
# follows is the variant, recorded verbatim. Nothing is inferred or added.
_SPLIT = re.compile(r"[,;]")


def strip_comments(path):
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines(True)
    return "".join(l for l in lines if not l.lstrip().startswith("#"))


def read_table(path, commented):
    text = strip_comments(path) if commented else path.read_text(encoding="utf-8", errors="replace")
    return pd.read_csv(io.StringIO(text), dtype=str, keep_default_na=False)


def build_genomics(src):
    """Tidy the per-biopsy `genomic` column of the Table 2 source worksheet.

    One row per biopsy per reported gene. This is the driver panel recorded in the
    study worksheet at the time of the report; no other sequencing output is used.
    """
    t2 = read_table(src, commented=True)
    if "genomic" not in t2.columns:
        raise SystemExit(f"{src.name}: no 'genomic' column")
    rows = []
    for _, r in t2.iterrows():
        cell = (r.get("genomic") or "").strip()
        if not cell or cell.upper() in {"N/A", "NA", "NONE", "-"}:
            continue
        for token in (t.strip() for t in _SPLIT.split(cell)):
            if not token:
                continue
            gene, _, variant = token.partition(" ")
            rows.append({"sample": r.get("sample", ""),
                         "gene": gene.strip(),
                         "variant": variant.strip()})
    out = pd.DataFrame(rows, columns=["sample", "gene", "variant"])
    return out.drop_duplicates().sort_values(["sample", "gene"]).reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True, type=Path, help="restricted de-identified directory")
    ap.add_argument("--out", required=True, type=Path, help="repository data directory")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    written, missing = [], []

    for name in COPY_AS_IS:
        src = a.source / name
        if not src.exists():
            missing.append(name)
            continue
        (a.out / name).write_text(
            strip_comments(src) if name in COMMENTED else src.read_text(encoding="utf-8", errors="replace"),
            encoding="utf-8",
        )
        written.append((name, "verbatim", ""))

    for name, drop in DROP_COLUMNS.items():
        src = a.source / name
        if not src.exists():
            missing.append(name)
            continue
        df = read_table(src, commented=name in COMMENTED)
        present = [c for c in drop if c in df.columns]
        df.drop(columns=present).to_csv(a.out / name, index=False)
        written.append((name, "columns removed", ", ".join(present) or "none"))

    gsrc = a.source / GENOMIC_SOURCE
    if gsrc.exists():
        g = build_genomics(gsrc)
        g.to_csv(a.out / GENOMIC_OUT, index=False)
        written.append((GENOMIC_OUT, "derived", f"{len(g)} gene calls from {GENOMIC_SOURCE}"))

    print(f"{'file':46s} {'action':16s} detail")
    for n, act, det in sorted(written):
        print(f"{n:46s} {act:16s} {det}")
    print(f"\n{len(written)} files written to {a.out}")

    print("\nwithheld by release policy:")
    for n, why in sorted(WITHHELD.items()):
        print(f"  {n:46s} {why}")
    if missing:
        print("\nnot found in source (skipped):")
        for n in missing:
            print(f"  {n}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
