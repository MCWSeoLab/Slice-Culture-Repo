#!/usr/bin/env python3
"""Compare the regenerated results against the values reported in the manuscript.

Reads the stage logs written by run_all.py and checks each reported statistic.
Numeric checks use an absolute tolerance chosen to be tighter than the precision
at which the value is printed in the report, so a genuine change in the analysis
fails the check rather than rounding to agreement.

    python run_all.py && python code/checks/verify_results.py
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

# stage log, description, regex capturing one number, reported value, tolerance
CHECKS = [
    ("fig2_mts_vs_histology", "Figure 2A nuclear density rho",
     r"NuclearDensity_%\s+rho=([+-][\d.]+)", 0.743, 0.001),
    ("fig2_mts_vs_histology", "Figure 2A exact p",
     r"NuclearDensity_%.*?exact p=([\d.]+)", 0.001, 0.001),
    ("fig2_mts_vs_histology", "Figure 2A slices",
     r"NuclearDensity_%.*?(\d+) slices", 33, 0),
    ("fig2_mts_vs_histology", "Figure 2A biopsies",
     r"NuclearDensity_%.*?(\d+) clusters", 13, 0),
    ("fig2_mts_vs_histology", "Figure 2B Ki-67 rho",
     r"Ki67_%\s+rho=([+-][\d.]+)", 0.651, 0.001),
    ("fig2_mts_vs_histology", "Figure 2B exact p",
     r"Ki67_%.*?exact p=([\d.]+)", 0.015, 0.001),
    ("fig2_mts_vs_histology", "Figure 2C cleaved caspase-3 rho",
     r"Caspase3_%\s+rho=([+-][\d.]+)", -0.548, 0.001),
    ("fig2_mts_vs_histology", "Figure 2C exact p",
     r"Caspase3_%.*?exact p=([\d.]+)", 0.020, 0.001),

    ("fig3_concordance", "Figure 3A exact p",
     r"panel A: exact p=([\d.]+)", 0.20, 0.005),
    ("fig3_concordance", "Figure 3B rho",
     r"panel B: rho=([+-][\d.]+)", -0.88, 0.005),
    ("fig3_concordance", "Figure 3B CI lower bound",
     r"panel B:.*?CI ([+-][\d.]+)", -1.00, 0.005),
    ("fig3_concordance", "Figure 3B CI upper bound",
     r"panel B:.*?CI [+-][\d.]+ to ([+-][\d.]+)", -0.51, 0.005),
    ("fig3_concordance", "Figure 3B exact p",
     r"panel B:.*?exact p=([\d.]+)", 0.06, 0.005),
    ("fig3_concordance", "Figure 3B courses",
     r"panel B:.*?n=(\d+) courses", 9, 0),
    ("fig3_concordance", "Figure 3B patients",
     r"panel B:.*?n=\d+ courses (\d+) patients", 6, 0),

    ("fig3_confusion", "2x2 true positives", r"TP (\d+)\s+FP", 4, 0),
    ("fig3_confusion", "2x2 false positives", r"TP \d+\s+FP (\d+)", 1, 0),
    ("fig3_confusion", "2x2 false negatives", r"FN (\d+)\s+TN", 1, 0),
    ("fig3_confusion", "2x2 true negatives", r"FN \d+\s+TN (\d+)", 4, 0),
    ("fig3_confusion", "2x2 sensitivity (%)", r"Sensitivity\s+\d+/\d+ =\s+(\d+)%", 80, 0),
    ("fig3_confusion", "2x2 accuracy (%)", r"Accuracy\s+\d+/\d+ =\s+(\d+)%", 80, 0),
    ("fig3_confusion", "2x2 Fisher exact p", r"Fisher exact p = ([\d.]+)", 0.206, 0.001),

    ("threshold_figure", "Threshold, DMSO CV route (% inhibition)",
     r"CV route\s+median [\d.]+%\s+-> SD [\d.]+ pp -> ([\d.]+)%", 31.5, 0.2),
    ("threshold_figure", "Threshold, half-normal route (% inhibition)",
     r"half-normal route sigma [\d.]+ pp.*?-> ([\d.]+)%", 34.1, 0.2),
    ("threshold_figure", "Half-normal sigma (percentage points)",
     r"half-normal route sigma ([\d.]+) pp", 17.0, 0.1),

    ("qc_gates", "Negative absorbance readings (%)",
     r"negative absorbance readings: \d+ of \d+ \(([\d.]+)%\)", 2.1, 0.05),
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--logs", type=Path,
                    default=Path(__file__).resolve().parents[2] / "results" / "logs")
    a = ap.parse_args()

    if not a.logs.is_dir():
        print(f"no logs at {a.logs}; run run_all.py first")
        return 2

    cache: dict[str, str] = {}
    width = max(len(d) for _, d, _, _, _ in CHECKS)
    passed = failed = missing = 0

    for stage, desc, rx, expected, tol in CHECKS:
        if stage not in cache:
            f = a.logs / f"{stage}.log"
            cache[stage] = f.read_text(encoding="utf-8", errors="replace") if f.exists() else ""
        m = re.search(rx, cache[stage], re.S)
        if not m:
            print(f"  MISSING {desc:{width}s}  not found in {stage}.log")
            missing += 1
            continue
        got = float(m.group(1))
        if abs(got - expected) <= tol:
            print(f"  ok      {desc:{width}s}  {got:g}")
            passed += 1
        else:
            print(f"  DIFFERS {desc:{width}s}  got {got:g}, report says {expected:g}")
            failed += 1

    print(f"\n{passed} reproduced, {failed} differ, {missing} not found")
    return 0 if failed == 0 and missing == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
