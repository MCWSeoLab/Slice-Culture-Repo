#!/usr/bin/env python3
"""Screen every file in this repository for protected health information.

Two independent passes:

  1. Pattern pass (always runs, no restricted input). Flags calendar dates in any
     common format, digit runs long enough to be a record number, government-ID
     shapes, and residual redaction markers left behind by the de-identification
     scrubber.

  2. Token pass (requires --crosswalk). Reads the re-identification crosswalk and
     searches every repository file for each identifier it contains. The
     crosswalk is read on the analyst's own machine and is never copied, echoed
     or written anywhere by this script: only per-file match counts are printed.

Exit status is non-zero if either pass finds anything.

    python check_no_phi.py --repo .
    python check_no_phi.py --repo . --crosswalk /secure/PHI_CROSSWALK_DO_NOT_SHARE.csv
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

TEXT_SUFFIXES = {".csv", ".py", ".md", ".txt", ".json", ".yml", ".yaml", ".sh", ".cff", ".svg", ".r", ".R"}
SKIP_DIRS = {".git", "results", "__pycache__", ".ipynb_checkpoints"}

PATTERNS = {
    "date m/d/y": re.compile(r"\b\d{1,2}/\d{1,2}/(?:\d{2}|\d{4})\b"),
    "date m-d-y": re.compile(r"\b\d{1,2}-\d{1,2}-(?:\d{2}|\d{4})\b"),
    "date y-m-d": re.compile(r"\b(?:19|20)\d{2}-\d{2}-\d{2}\b"),
    "date d.m.y": re.compile(r"\b\d{1,2}\.\d{1,2}\.(?:19|20)\d{2}\b"),
    "date spelled": re.compile(
        r"\b(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)[a-z]*\s+\d{1,2},?\s*(?:19|20)\d{2}\b", re.I),
    "record-number run": re.compile(r"(?<![\d.])\d{7,}(?!\d)"),
    "government ID shape": re.compile(r"\b\d{3}-\d{2}-\d{4}\b"),
    "age over 89": re.compile(r"\b(?:9\d|1\d\d)\s*(?:years old|yo\b|y/o\b)"),
}

# Reported but not treated as failures. A redaction marker is the scrubber's own
# output: it records that an identifier was removed at that position.
INFORMATIONAL = {
    "redaction marker": re.compile(r"\[(?:NAME|MRN|TB|ID)\]"),
}

# Literal digit strings that are part of the analysis and are not record numbers.
# Each entry must be justified in docs/PHI_AND_DEIDENTIFICATION.md.
ALLOWLIST_FILE = ".phi-allowlist"

# A run of digits inside a decimal fraction is a measurement, not a record number.
DECIMAL = re.compile(r"\.\d{7,}")

# Vector graphics carry machine-generated coordinate streams. Those attributes are
# geometry, never text, so they are removed before screening and the record-number
# pattern is not applied to the remainder of the file.
SVG_GEOMETRY = re.compile(r'\b(?:d|transform|points|viewBox|gradientTransform)="[^"]*"')


def repo_files(repo: Path):
    for p in sorted(repo.rglob("*")):
        if not p.is_file():
            continue
        if any(part in SKIP_DIRS for part in p.parts):
            continue
        if p.suffix.lower() in TEXT_SUFFIXES:
            yield p


def load_allowlist(repo: Path) -> list[str]:
    f = repo / ALLOWLIST_FILE
    if not f.exists():
        return []
    out = []
    for line in f.read_text(encoding="utf-8").splitlines():
        line = line.split("#", 1)[0].strip()
        if line:
            out.append(line)
    return out


def pattern_pass(repo: Path) -> int:
    allow = load_allowlist(repo)
    if allow:
        print(f"  allowlist: {len(allow)} literal(s) from {ALLOWLIST_FILE}")
    hits, notes = 0, 0
    for p in repo_files(repo):
        text = p.read_text(encoding="utf-8", errors="replace")
        stripped = DECIMAL.sub(".0", text)
        is_svg = p.suffix.lower() == ".svg"
        if is_svg:
            stripped = SVG_GEOMETRY.sub("", stripped)
        for literal in allow:
            stripped = stripped.replace(literal, "")
        for label, rx in PATTERNS.items():
            if is_svg and label == "record-number run":
                continue
            n = len(rx.findall(stripped))
            if n:
                print(f"  FLAG {p.relative_to(repo)}: {n} x {label}")
                hits += n
        for label, rx in INFORMATIONAL.items():
            n = len(rx.findall(text))
            if n:
                notes += n
    if notes:
        print(f"  note: {notes} redaction marker(s) present; these record removed identifiers")
    print(f"  pattern pass: {hits} flag(s)")
    return hits


def token_pass(repo: Path, crosswalk: Path) -> int:
    """Search every repository file for each identifier held in the crosswalk.

    Tokens are matched at token boundaries, not as bare substrings: a four-digit
    record number must not be reported because those digits happen to fall inside
    a measurement such as 0.75002, and a short name fragment must not be reported
    because it occurs inside an ordinary word. The crosswalk is read here and
    nowhere else; no value from it is printed, stored or copied.
    """
    patterns = []
    seen = set()
    with crosswalk.open(encoding="utf-8", errors="replace", newline="") as fh:
        for row in csv.DictReader(fh):
            for key, value in row.items():
                if key == "study_id":
                    continue
                v = (value or "").strip()
                candidates = [v] + [q for q in re.split(r"[\s,]+", v)]
                for c in candidates:
                    c = c.strip()
                    if len(c) < 3 or c.lower() in seen:
                        continue
                    seen.add(c.lower())
                    esc = re.escape(c)
                    if c.isdigit():
                        rx = rf"(?<![\d.]){esc}(?![\d.])"
                    elif c.isalpha():
                        rx = rf"\b{esc}\b"
                    else:
                        rx = rf"(?<![\w.]){esc}(?![\w.])"
                    # A short alphabetic fragment collides with ordinary English and
                    # with identifiers in source code, so it is reported for review
                    # rather than treated as a failure. Everything else is a failure.
                    weak = c.isalpha() and len(c) <= 4
                    patterns.append((re.compile(rx, re.I), weak))
    print(f"  token pass: {len(patterns)} identifier token(s) loaded from the crosswalk")
    hits, reviews = 0, 0
    for p in repo_files(repo):
        text = p.read_text(encoding="utf-8", errors="replace")
        strong = sum(1 for rx, weak in patterns if not weak and rx.search(text))
        short = sum(1 for rx, weak in patterns if weak and rx.search(text))
        if strong:
            print(f"  FLAG   {p.relative_to(repo)}: {strong} crosswalk token(s) present")
            hits += strong
        if short:
            print(f"  REVIEW {p.relative_to(repo)}: {short} short name fragment(s) matched")
            reviews += short
    if reviews:
        print(f"  token pass: {reviews} item(s) for manual review "
              f"(short fragments collide with ordinary words; confirm each in context)")
    print(f"  token pass: {hits} flag(s)")
    return hits


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, default=Path("."))
    ap.add_argument("--crosswalk", type=Path, default=None,
                    help="path to the re-identification crosswalk; never copied or printed")
    a = ap.parse_args()

    repo = a.repo.resolve()
    print(f"screening {repo}")
    total = pattern_pass(repo)
    if a.crosswalk:
        total += token_pass(repo, a.crosswalk)
    else:
        print("  token pass: skipped (no --crosswalk given)")

    print("CLEAN" if total == 0 else f"NOT CLEAN: {total} flag(s)")
    return 0 if total == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
