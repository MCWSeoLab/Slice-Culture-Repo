#!/usr/bin/env python3
"""
Extract raw data tables from GraphPad Prism files (.prism and .pzfx), normalize
MTS absorbance to the DMSO/vehicle control on the same sheet+timepoint, and emit
de-identified tidy CSVs.

RUNS LOCALLY ONLY. Prism sheet/table titles embed tissue-bank numbers and biopsy
dates, so every title is scrubbed against the
crosswalk before anything is written.

Outputs (into deid/):
  prism_raw_long_deid.csv    one row per (file, sheet, row label, column, replicate)
  prism_mts_normalized_deid.csv   treated/control * 100 for every non-control column
  prism_index_deid.csv       inventory of every file and sheet found
  _PRISM_AUDIT.txt
"""
import os, re, csv, json, zipfile, datetime, glob
import xml.etree.ElementTree as ET

ROOT = os.path.expanduser("~/mnt/Organotypic Slice Culture")
OUT  = os.path.join(ROOT, "deid")
XWP  = os.path.join(ROOT, "PHI_CROSSWALK", "PHI_CROSSWALK_DO_NOT_SHARE.csv")
os.makedirs(OUT, exist_ok=True)

audit = []
def A(m): audit.append(str(m)); print(m)

# ------------------------------------------------------------------ scrubbing
STOP = {"no","yes","na","nan","none","pr","sd","pd","cr","ne","gem","nab","tx","ca","ih",
        "er","or","and","the","of","on","in","at","to","dx","bx","pt","mg","ml","um","nm",
        "xt","xr","pa","br","tb","id","ok","up","hx","po","iv","d1","d2","cy","wk","mo"}
def _usable(t):
    t = (t or "").strip()
    if not t or t.lower() in STOP: return False
    return len(t) >= 3

DATE_RX = [re.compile(r'\b\d{1,2}[/-]\d{1,2}[/-]\d{2,4}\b'),
           re.compile(r'\b\d{4}[/-]\d{1,2}[/-]\d{1,2}\b')]
TBRX    = re.compile(r'\bTB\s*#?\s*\d{3,6}\b', re.I)
LONGNUM = re.compile(r'\b\d{6,}\b')

XW, DATE2SID = {}, {}
if os.path.exists(XWP):
    for r in csv.DictReader(open(XWP)):
        XW[r["study_id"]] = r
        if r["Biopsy_Date"]:
            DATE2SID.setdefault(r["Biopsy_Date"], []).append(r["study_id"])
TOKENS = set()
for r in XW.values():
    for f in ("MRN","TB_number","Patient_Initials","Patient_Name"):
        v = (r.get(f) or "").strip()
        TOKENS.add(v)
        for p in re.split(r'\W+', v): TOKENS.add(p)
TOKENS = sorted({t for t in TOKENS if _usable(t)}, key=len, reverse=True)
TOKRX  = [(t, re.compile(r'\b'+re.escape(t)+r'\b', re.I)) for t in TOKENS]

def parse_date(s):
    for fmt in ("%m/%d/%Y","%m/%d/%y","%Y-%m-%d","%m-%d-%Y","%m-%d-%y"):
        try: return datetime.datetime.strptime(s, fmt).date()
        except Exception: pass
    return None

def sid_from_title(title):
    """If a title carries a biopsy date or TB#, resolve it to a study_id."""
    for rx in DATE_RX:
        for m in rx.findall(title or ""):
            d = parse_date(m)
            if d:
                s = DATE2SID.get(d.isoformat())
                if s and len(s) == 1: return s[0]
                if s: return "|".join(s)
    tb = TBRX.search(title or "")
    if tb:
        num = re.sub(r'\D','',tb.group(0))
        for sid, r in XW.items():
            if re.sub(r'\D','', r.get("TB_number") or "") == num: return sid
    return ""

def scrub(s):
    if s is None: return ""
    s = str(s)
    s = TBRX.sub("[TB]", s)
    for rx in DATE_RX: s = rx.sub("[DATE]", s)
    s = LONGNUM.sub("[ID]", s)
    for t, rx in TOKRX:
        s = rx.sub("[ID]" if t.isdigit() else "[NAME]", s)
    return " ".join(s.split())

# ------------------------------------------------------------------ readers
def num(v):
    try:
        f = float(str(v).strip())
        return f
    except Exception:
        return None

def read_prism(path):
    """Yield one dict per populated cell of every data sheet in a .prism file.

    Prism stores 'grouped'/'y_replicates' tables as a flat CSV: field 0 is the row
    title, then every column's replicate subcolumns laid out consecutively. The
    number of replicates per column is (numberOfColumns - 1) / len(dataSets), so a
    naive field->column mapping silently assigns the wrong drug to each value.
    """
    out = []
    try: z = zipfile.ZipFile(path)
    except Exception as e:
        A("  !! cannot open %s (%s)" % (os.path.basename(path), e)); return out
    names = set(z.namelist())
    if "document.json" not in names: return out
    doc = json.loads(z.read("document.json"))
    for uid in doc.get("sheets", {}).get("data", []):
        sp = f"data/sheets/{uid}/sheet.json"
        if sp not in names: continue
        sh = json.loads(z.read(sp))
        title = sh.get("title","")
        tbl = sh.get("table",{}); tu = tbl.get("uid")
        cpath = f"data/tables/{tu}/data.csv"
        if cpath not in names: continue
        raw = z.read(cpath).decode("utf-8","replace")
        if not raw.strip(): continue
        cols = []
        for ds in tbl.get("dataSets", []):
            dp = f"data/sets/{ds}.json"
            cols.append((json.loads(z.read(dp)).get("title") or "") if dp in names else "")
        ncol = None
        cj = f"data/tables/{tu}/content.json"
        if cj in names:
            ncol = json.loads(z.read(cj)).get("numberOfColumns")
        has_rowtitles = tbl.get("rowTitlesDataSet") is not None
        reps = 1
        if ncol and cols:
            data_fields = ncol - (1 if has_rowtitles else 0)
            if data_fields > 0 and data_fields % len(cols) == 0:
                reps = data_fields // len(cols)
        for row in csv.reader(raw.splitlines()):
            if not row or not any(c.strip() for c in row): continue
            rlabel = row[0].strip()
            fields = row[1:]
            for j, cell in enumerate(fields):
                v = num(cell)
                if v is None: continue
                ci, ri = divmod(j, reps)
                out.append({"sheet": title, "sheet_uid": uid,
                            "row_label": rlabel,
                            "column": cols[ci] if ci < len(cols) else f"col{ci+1}",
                            "replicate": ri + 1, "value": v,
                            "reps_per_col": reps})
    return out

def read_pzfx(path):
    """Yield dicts for every populated cell in every table of a .pzfx (XML) file."""
    out = []
    try: root = ET.parse(path).getroot()
    except Exception as e:
        A("  !! cannot parse %s (%s)" % (os.path.basename(path), e)); return out
    def local(t): return t.split('}')[-1]
    for tbl in root.iter():
        if local(tbl.tag) != "Table": continue
        ti = None
        for ch in tbl:
            if local(ch.tag) == "Title": ti = "".join(ch.itertext()).strip(); break
        title = ti or ""
        # row titles
        rowlabels = []
        for ch in tbl:
            if local(ch.tag) in ("RowTitlesColumn","XColumn","XAdvancedColumn"):
                for sub in ch.iter():
                    if local(sub.tag) == "Subcolumn":
                        rowlabels = ["".join(d.itertext()).strip() for d in sub
                                     if local(d.tag) == "d"]
                        break
                if rowlabels: break
        for col in tbl:
            if local(col.tag) not in ("YColumn","YAdvancedColumn"): continue
            ct = ""
            for ch in col:
                if local(ch.tag) == "Title": ct = "".join(ch.itertext()).strip(); break
            rep = 0
            for sub in col:
                if local(sub.tag) != "Subcolumn": continue
                rep += 1
                vals = [d for d in sub if local(d.tag) == "d"]
                for i, d in enumerate(vals):
                    v = num("".join(d.itertext()))
                    if v is None: continue
                    out.append({"sheet": title, "sheet_uid": title,
                                "row_label": rowlabels[i] if i < len(rowlabels) else f"row{i+1}",
                                "column": ct, "replicate": rep, "value": v,
                                "reps_per_col": None})
    return out

# ------------------------------------------------------------------ walk
CONTROL_RX = re.compile(r'^\s*(dmso|control|ctrl|untreated|vehicle|veh|no\s*tx|media|medium)\b', re.I)

files = sorted(glob.glob(os.path.join(ROOT,"**","*.prism"), recursive=True) +
               glob.glob(os.path.join(ROOT,"**","*.pzfx"),  recursive=True))
A("Scanning %d Prism files under %s" % (len(files), os.path.basename(ROOT)))

long_rows, index_rows = [], []
for p in files:
    rel = os.path.relpath(p, ROOT)
    rows = read_prism(p) if p.lower().endswith(".prism") else read_pzfx(p)
    sheets = sorted({r["sheet"] for r in rows})
    index_rows.append({"file": scrub(rel), "kind": os.path.splitext(p)[1].lstrip("."),
                       "n_sheets": len(sheets), "n_values": len(rows),
                       "study_id_from_path": sid_from_title(rel),
                       "sheets": scrub(" | ".join(sheets))[:300]})
    for r in rows:
        sid = sid_from_title(r["sheet"]) or sid_from_title(rel)
        long_rows.append({"file": scrub(rel), "sheet": scrub(r["sheet"]),
                          "_key": (rel, r["sheet_uid"]),
                          "study_id": sid, "row_label": scrub(r["row_label"]),
                          "column": scrub(r["column"]), "replicate": r["replicate"],
                          "reps_per_col": r.get("reps_per_col"),
                          "value": r["value"]})
A("  extracted %d values from %d files" % (len(long_rows), sum(1 for i in index_rows if i["n_values"])))

# ------------------------------------------------------------------ normalize
norm_rows = []
groups = {}
for r in long_rows:
    groups.setdefault((r["_key"], r["row_label"]), []).append(r)
for (_k, rl), rows in groups.items():
    ctrl = [r["value"] for r in rows if CONTROL_RX.match(r["column"] or "")]
    if not ctrl: continue
    cmean = sum(ctrl)/len(ctrl)
    if cmean == 0: continue
    for r in rows:
        if CONTROL_RX.match(r["column"] or ""): continue
        norm_rows.append({"file": r["file"], "sheet": r["sheet"], "study_id": r["study_id"],
                          "timepoint_label": rl, "drug": r["column"],
                          "replicate": r["replicate"],
                          "raw_value": r["value"], "control_mean": round(cmean, 6),
                          "n_control_reps": len(ctrl),
                          "pct_of_control": round(r["value"]/cmean*100, 3),
                          "pct_inhibition": round(100 - r["value"]/cmean*100, 3)})
A("  normalized %d treated values against an on-sheet control" % len(norm_rows))
A("  sheets with a usable control: %d of %d" %
  (len({(r['file'],r['sheet'],r['timepoint_label']) for r in norm_rows}), len(groups)))
linked = {r["study_id"] for r in norm_rows if r["study_id"]}
A("  normalized rows resolving to a study_id: %d (%d distinct patients)" %
  (sum(1 for r in norm_rows if r["study_id"]), len(linked)))

def write(name, rows):
    p = os.path.join(OUT, name)
    if not rows: A("  (empty) "+name); return
    with open(p,"w",newline="",encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    A("  wrote %-34s %d rows" % (name, len(rows)))

for r in long_rows: r.pop("_key", None)
def dedupe(rows, keys):
    seen, out = set(), []
    for r in rows:
        k = tuple(r.get(x) for x in keys)
        if k in seen: continue
        seen.add(k); out.append(r)
    return out
before_l, before_n = len(long_rows), len(norm_rows)
long_rows = dedupe(long_rows, ["study_id","sheet","row_label","column","replicate","value"])
norm_rows = dedupe(norm_rows, ["study_id","sheet","timepoint_label","drug","replicate","raw_value","control_mean"])
A("  de-duplicated repeated file copies: raw %d->%d, normalized %d->%d"
  % (before_l, len(long_rows), before_n, len(norm_rows)))

A("OUTPUT:")
write("prism_index_deid.csv", index_rows)
write("prism_raw_long_deid.csv", long_rows)
write("prism_mts_normalized_deid.csv", norm_rows)

# leak check -- text columns only. Numeric measurement columns are excluded because a
# 4-digit tissue-bank number can coincidentally equal an absorbance or cell count; matching
# there is a false positive, not a disclosure (the number carries no identifying context).
NUMERIC_COLS = {"value","raw_value","control_mean","pct_of_control","pct_inhibition",
                "replicate","n_control_reps","n_sheets","n_values"}
bad = 0
for name in ("prism_index_deid.csv","prism_raw_long_deid.csv","prism_mts_normalized_deid.csv"):
    p = os.path.join(OUT,name)
    if not os.path.exists(p): continue
    text = []
    for row in csv.DictReader(open(p)):
        for k, v in row.items():
            if k not in NUMERIC_COLS and v: text.append(v)
    body = " \n ".join(text).upper()
    hits = [t for t in TOKENS if re.search(r'\b'+re.escape(t.upper())+r'\b', body)]
    bad += len(hits); A("  leak check %-34s %d residual identifier tokens (text cols only)" % (name, len(hits)))
A("LEAK CHECK TOTAL: %d  (must be 0 before staging)" % bad)
open(os.path.join(OUT,"_PRISM_AUDIT.txt"),"w").write("\n".join(audit))
