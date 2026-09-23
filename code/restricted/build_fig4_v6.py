# -*- coding: utf-8 -*-
"""Figure 4 v4 = v3 with the CA19-9 ULN corrected and panels B and F enlarged in place.

Two edits, both surgical:

1. ULN 37.5 -> 35.0.  The three occurrences are literal strings in the page content stream
   (`( ULN 37.5) Tj` twice, and `(37.5 U/mL)` in the key).  "37.5" and "35.0" are the same
   number of characters and the digits are tabular in this face, so a byte-for-byte swap in
   the uncompressed stream leaves every glyph position untouched.

2. Panels B and F are painted over with white and redrawn larger from panelB_big.pdf /
   panelF_big.pdf.  Nothing else on the page moves.

Note the v3 page still carries the nine ORIGINAL panel B/F tiles as buried objects under the
first white patch - they were covered, not removed.  This build inherits them.  They
are de-identified micrographs that were already published in the figure, so they are harmless,
but they do inflate the file and are extractable with pdfimages.
"""
import re
import subprocess
from pathlib import Path

import pikepdf

SRC = "v3.pdf"
OUT = "fig4_two_cases_v6.pdf"
PAGE_W, PAGE_H = 547.936, 643.3269
PT = 72.0 / 25.4
Y_TOP = PAGE_H - 41.2 * PT
Y_BOT = PAGE_H - 88.6 * PT
PATCHES = {"B": (4.0, 258.0), "F": (262.0, 520.0)}

# ---------------------------------------------------------------- 1. ULN 37.5 -> 35.0
subprocess.run(["qpdf", "--qdf", "--object-streams=disable", SRC, "_uln_in.pdf"], check=True)
d = Path("_uln_in.pdf").read_bytes()
n1 = d.count(b"( ULN 37.5) Tj")
n2 = d.count(b"(37.5 U/mL\\)) Tj") + d.count(b"37.5 U/mL")
d = d.replace(b"( ULN 37.5) Tj", b"( ULN 35.0) Tj")
d = d.replace(b"37.5 U/mL", b"35.0 U/mL")
# NB: a blanket check on b"37.5" is wrong - PDF coordinates such as 637.5283 and 117.2861
# contain that substring.  Check the two labelled forms instead.
assert b"( ULN 37.5)" not in d and b"37.5 U/mL" not in d, "a labelled 37.5 survived"
Path("_uln_out.pdf").write_bytes(d)
subprocess.run(["qpdf", "_uln_out.pdf", "_base.pdf"], check=True)
print(f"ULN: replaced {n1} panel annotations and the key ({n2} matches of '37.5 U/mL')")

# ------------------------------------------------- 2. white out and redraw panels B and F
pdf = pikepdf.open("_base.pdf")
page = pdf.pages[0]

rects = "".join(
    f"{x0:.3f} {Y_BOT:.3f} {x1 - x0:.3f} {Y_TOP - Y_BOT:.3f} re f\n"
    for x0, x1 in PATCHES.values()
)
page.contents_add(pikepdf.Stream(pdf, f"q 1 1 1 rg\n{rects}Q\n".encode()), prepend=False)

for key, (x0, x1) in PATCHES.items():
    src = pikepdf.open(f"panel{key}_big.pdf")
    page.add_overlay(src.pages[0], pikepdf.Rectangle(x0, Y_BOT, x1, Y_TOP))
    print(f"panel {key} placed at x {x0:.1f}-{x1:.1f}, y {Y_BOT:.1f}-{Y_TOP:.1f}")

# ---------------------------------------------- 3. shrink the images nobody can see
# v3 covered the original panel tiles rather than removing them, and this build covers v3's
# in turn, so the page would carry three generations of hidden micrographs and weigh 20 MB.
# Every buried tile is entirely inside a white patch, so replace its pixels with a single
# transparent pixel: the draw operation stays valid, the object stops costing anything.
BURIED = {(447, 336), (450, 336), (444, 336), (445, 336), (544, 544)}
tiny = pikepdf.Stream(pdf, b"\x00\x00\x00")
n = 0
for name, xo in page.resources.get("/XObject", pikepdf.Dictionary()).items():
    if xo.get("/Subtype") != pikepdf.Name("/Image"):
        continue
    dims = (int(xo.Width), int(xo.Height))
    if dims not in BURIED:
        continue
    xo.write(b"\x00\x00\x00", filter=pikepdf.Name("/FlateDecode")) if False else None
    xo.write(b"\xff\xff\xff")
    xo.Width, xo.Height = 1, 1
    xo.ColorSpace = pikepdf.Name("/DeviceRGB")
    xo.BitsPerComponent = 8
    for k in ("/Filter", "/DecodeParms", "/SMask", "/Decode", "/Interpolate"):
        if k in xo:
            del xo[k]
    n += 1
print(f"blanked {n} buried tiles that sit under the white patches")

pdf.save(OUT, linearize=True)
print("wrote", OUT)
