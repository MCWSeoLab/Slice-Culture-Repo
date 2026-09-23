#!/usr/bin/env python3
"""Reproduce the analysis and figures of this report from the released data.

Every stage runs in a scratch working directory (results/work) populated with a
copy of data/, so the released tables are never modified. Console output for each
stage is written to results/logs/, regenerated figures to results/figures/, and
regenerated tables to results/tables/.

    python run_all.py              # run everything
    python run_all.py --list       # show the stages without running them
    python run_all.py --only fig3  # run the stages whose name contains "fig3"
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CODE, DATA, RESULTS = ROOT / "code", ROOT / "data", ROOT / "results"
WORK = RESULTS / "work"

# name, script, environment overrides
STAGES = [
    ("qc_gates",            "analysis/qc_gates.py", {}),
    ("threshold_figure",    "figures/make_threshold_figure.py", {}),
    ("clustered_p",         "analysis/clustered_p.py", {}),
    ("recompute",           "analysis/recompute.py", {}),
    ("sensitivity",         "analysis/sensitivity.py", {}),
    ("cellularity",         "analysis/adjust_for_cellularity.py", {}),
    ("fig2_mts_vs_histology", "figures/make_fig2_v2.py", {}),
    ("fig3_concordance",    "figures/make_fig3_v5.py",
     {"FIG3_OUTCOME": "response", "FIG3_LAYOUT": "row",
      "FIG3_STATS": "computed", "FIG3_TICKS": "short"}),
    ("fig3_confusion",      "figures/make_confusion_response.py", {}),
    ("figS_tumor_content",  "figures/make_figS_cellularity.py", {}),
    ("figS_replicate_agreement", "figures/make_fig_reproducibility.py", {}),
    ("fig_flow_diagram",    "figures/make_fig_flow.py", {}),
    ("sap_forest",          "figures/make_fig_forest.py", {}),
    ("sap_exvivo_vs_recist", "figures/make_fig_exvivo_vs_recist.py", {}),
    ("sap_exvivo_vs_ttnt",  "figures/make_fig_exvivo_vs_ttnt.py", {}),
]


def prepare_work() -> None:
    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir(parents=True)
    for csv in sorted(DATA.glob("*.csv")):
        shutil.copy2(csv, WORK / csv.name)
    # Some scripts address the data directory by name from the project root.
    link = WORK / "deid"
    try:
        link.symlink_to(".", target_is_directory=True)
    except (OSError, NotImplementedError):
        (WORK / "deid").mkdir(exist_ok=True)
        for csv in sorted(DATA.glob("*.csv")):
            shutil.copy2(csv, WORK / "deid" / csv.name)


def run(name: str, script: str, overrides: dict) -> tuple[str, float, str]:
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(CODE / "lib")] + ([env["PYTHONPATH"]] if env.get("PYTHONPATH") else []))
    env.setdefault("MPLBACKEND", "Agg")
    env.update(overrides)
    log = RESULTS / "logs" / f"{name}.log"
    t0 = time.time()
    proc = subprocess.run([sys.executable, str(CODE / script)], cwd=WORK, env=env,
                          capture_output=True, text=True)
    log.write_text(proc.stdout + proc.stderr, encoding="utf-8")
    elapsed = time.time() - t0
    if proc.returncode == 0:
        return "ok", elapsed, ""
    tail = (proc.stderr or proc.stdout).strip().splitlines()
    return "FAILED", elapsed, tail[-1] if tail else f"exit {proc.returncode}"


def collect() -> None:
    for sub, pattern in (("figures", "figures/*"), ("tables", "*.csv")):
        dest = RESULTS / sub
        dest.mkdir(parents=True, exist_ok=True)
        for f in sorted(WORK.glob(pattern)):
            if f.is_file():
                shutil.copy2(f, dest / f.name)
    # Tables that were copied in as inputs are not outputs.
    for f in DATA.glob("*.csv"):
        stale = RESULTS / "tables" / f.name
        if stale.exists() and stale.read_bytes() == f.read_bytes():
            stale.unlink()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--only", default=None, help="substring filter on the stage name")
    a = ap.parse_args()

    stages = [s for s in STAGES if not a.only or a.only in s[0]]
    if a.list:
        for name, script, env in stages:
            extra = " ".join(f"{k}={v}" for k, v in env.items())
            print(f"  {name:28s} {script}{'  ' + extra if extra else ''}")
        return 0
    if not stages:
        print(f"no stage matches {a.only!r}")
        return 1

    (RESULTS / "logs").mkdir(parents=True, exist_ok=True)
    prepare_work()

    failures = 0
    print(f"{'stage':28s} {'status':8s} {'sec':>6s}")
    for name, script, env in stages:
        status, elapsed, why = run(name, script, env)
        print(f"{name:28s} {status:8s} {elapsed:6.1f}" + (f"   {why}" if why else ""))
        failures += status != "ok"
    collect()

    print(f"\nlogs     {RESULTS / 'logs'}")
    print(f"figures  {RESULTS / 'figures'}")
    print(f"tables   {RESULTS / 'tables'}")
    if failures:
        print(f"\n{failures} stage(s) failed")
        return 1
    print("\nall stages completed; run code/checks/verify_results.py to compare against the report")
    return 0


if __name__ == "__main__":
    sys.exit(main())
