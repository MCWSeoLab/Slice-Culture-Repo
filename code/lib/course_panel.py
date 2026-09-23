#!/usr/bin/env python3
"""Shared clinical-course panel: CA 19-9 on a linear axis with treatment bands above it.

Used by the Figure 4 case columns. Everything comes from the de-identified chart extracts;
days are relative to that patient's biopsy.

    draw_course(fig, gs[3, col], "PT-024", DRUGCOL, day_max=154, cap=600)

Layout, top to bottom, as three stacked axes sharing one x range:

    treatment strip   thin coloured bands, one row per concurrent course, best response after
    break region      OPTIONAL — only drawn when a value exceeds `cap`; carries those points
                      above a pair of break marks so one outlier cannot flatten the rest
    CA 19-9 trace     LINEAR y from 0 to `cap`, open circles, dashed ULN line,
                      RECIST letters along the foot

Two things this module is deliberate about:

*   **The y axis is linear.** It used to be log, which made a fall from 574 to 73 look the
    same size as a fall from 57 to 7. Linear is honest about magnitude, at the cost of
    needing a break when one draw is orders of magnitude above the rest.
*   **Truncation is visible.** `day_max` cuts the panel at a stated day. Any treatment
    course still running at that point is drawn with an open, arrowed right end rather than
    a stop tick, so a truncated course cannot be misread as a course that ended. The caller
    is responsible for saying in the legend where the cut is and why.
"""
import re

import numpy as np
import matplotlib.patheffects as pe
import pandas as pd
from matplotlib.patches import Rectangle

PITCH, BARH = 0.30, 0.16         # treatment-strip row pitch / bar height, axes fraction
ULN = 37.5                       # U/mL — the laboratory's stated upper limit of normal
TRACE = "#1f3b5c"

_MK = pd.read_csv("chart_markers_deid.csv")
_TH = pd.read_csv("chart_therapy_deid.csv")
_IM = pd.read_csv("chart_imaging_deid.csv")


def _canon(s):
    s = str(s).lower()
    if re.search(r"folfirinox|folforinox", s): return "FOLFIRINOX"
    if re.search(r"gem", s) and re.search(r"nab|pac|abrax", s): return "Gem/nab-Pac"
    if re.search(r"rmc", s): return "RMC-7977"
    if re.search(r"trametinib|cobimetinib", s): return "Trametinib"
    if re.search(r"vemurafenib", s): return "Vemurafenib"
    return None


def _short(s):
    """A compact label for the treatment strip. This is the drug the patient actually got,
    not the assay class it maps to — colour carries the class, text carries the agent."""
    t = re.sub(r"\s*\([^)]*\)", "", str(s)).strip()          # drop (neoadjuvant), (adjuvant)
    t = re.sub(r"\s*\btrial\b", "", t, flags=re.I)
    t = re.sub(r"\s*/\s*", "/", re.sub(r"\s+", " ", t)).strip()
    t = re.sub(r"(?i)gem/nab-?pac\w*", "Gem/nab-Pac", t)
    return t if len(t) <= 17 else t[:16] + "…"


_OTHER = ["0.52", "0.70", "0.38", "0.62"]


def _val(v):
    """returns (value, censor) — censor is '', 'above' or 'below'"""
    s = str(v).strip().replace(",", "")
    c = "above" if s.startswith(">") else ("below" if s.startswith("<") else "")
    try:    return float(s.lstrip("<>")), c
    except Exception: return np.nan, c


def _pack(th, span):
    """Assign each course to a strip row so concurrent courses do not overlap."""
    rows, ends = [], []
    for _, r in th.iterrows():
        st = r.day_start
        en = r.day_stop if pd.notna(r.day_stop) else st + max(span * 0.02, 8)
        for k, e in enumerate(ends):
            if st >= e:
                rows.append(k); ends[k] = en; break
        else:
            rows.append(len(ends)); ends.append(en)
    return rows, max(len(ends), 1)


def draw_course(fig, spec, sid, drugcol, day_max=None, cap=None, show_ylabel=True, fs=5.4):
    """Draw one patient's course into the gridspec cell `spec`.

    day_max  truncate the panel here; courses still running are drawn open-ended, and
             courses starting at or after the cut are omitted entirely
    cap      top of the linear y axis; anything above goes in a broken-out region

    Returns dict(ax_strip=, ax_hi=, ax_lo=, n_above=, xlim=).
    """
    th = _TH[_TH.study_id == sid].copy()
    th = th[~th.regimen.astype(str).str.contains(r"resect|whipple", case=False, na=False)]
    th = th.dropna(subset=["day_start"]).sort_values("day_start")
    mk = _MK[(_MK.study_id == sid) & _MK.test.astype(str).str.contains("19-9")].copy()
    mk[["v", "cen"]] = mk.result.apply(lambda x: pd.Series(_val(x)))
    mk = mk.dropna(subset=["v", "day"]).sort_values("day")
    im = _IM[(_IM.study_id == sid) & _IM.recist.notna()].sort_values("day")

    # --- truncation ------------------------------------------------------------------
    # Courses are kept if they START at or before the cut, so a regimen running through the
    # cut still appears; its bar is clipped and flagged as open-ended below.
    if day_max is not None:
        mk = mk[mk.day <= day_max]
        im = im[im.day <= day_max]
        # Strictly before the cut: a course starting exactly at `day_max` has zero duration
        # inside the window, so it would contribute a zero-width bar and a label and nothing
        # else. This is what keeps tovorafenib out of Figure 4H, whose cut IS its start day.
        th = th[th.day_start < day_max]

    # --- x range ---------------------------------------------------------------------
    ends_for_range = th.day_stop.clip(upper=day_max) if day_max is not None else th.day_stop
    days = pd.concat([th.day_start, ends_for_range, mk.day, im.day]).dropna()
    lo, hi = float(days.min()), float(days.max())
    if day_max is not None:
        hi = min(hi, day_max) if hi > day_max else hi
        hi = max(hi, day_max)
    span0 = max(hi - lo, 1.0)
    pad = max(span0 * 0.05, 10)
    XLO, XHI = lo - pad, hi + pad * 1.2

    # --- how tall must the strip be, and is a break needed? ---------------------------
    rows, nrow = _pack(th, span0)
    if cap is None:
        cap = float(mk.v.max()) * 1.15 if len(mk) else 1.0
    above = mk[mk.v > cap]
    n_above = len(above)

    strip_h = 0.30 + 0.26 * nrow          # relative height for the treatment strip
    hi_h = 0.30 if n_above else 0.0
    ratios = [strip_h] + ([hi_h] if n_above else []) + [1.55]
    sub = spec.subgridspec(len(ratios), 1, height_ratios=ratios, hspace=0.10)

    k = 0
    ax_strip = fig.add_subplot(sub[k]); k += 1
    ax_hi = fig.add_subplot(sub[k]) if n_above else None
    if n_above: k += 1
    ax_lo = fig.add_subplot(sub[k])

    for ax in filter(None, (ax_strip, ax_hi, ax_lo)):
        ax.set_xlim(XLO, XHI)

    # --- treatment strip --------------------------------------------------------------
    # Fit the y range to the rows actually in use. With a fixed 0-1 range a single-row
    # strip left most of its axes empty, opening a large gap above the trace.
    ax_strip.set_ylim(1.0 - PITCH * nrow - 0.10, 1.10)
    ax_strip.axis("off")
    oth = {}
    for (_, r), row in zip(th.iterrows(), rows):
        st = r.day_start
        en_true = r.day_stop if pd.notna(r.day_stop) else st + max(span0 * 0.02, 8)
        open_end = pd.isna(r.day_stop) or (day_max is not None and en_true > day_max)
        en = min(en_true, day_max) if day_max is not None else en_true
        y0 = 1.0 - PITCH * (row + 1)
        name, canon = _short(r.regimen), _canon(r.regimen)
        resp = str(r.best_response).strip() if pd.notna(r.best_response) else ""
        col = drugcol.get(canon, "0.6") if canon else \
            oth.setdefault(name, _OTHER[len(oth) % len(_OTHER)])
        ax_strip.add_patch(Rectangle((st, y0), max(en - st, span0 * 0.004), BARH,
                                     facecolor=col, edgecolor="black", lw=0.35, zorder=3))
        if open_end:
            ax_strip.annotate("", xy=(en + span0 * 0.030, y0 + BARH / 2),
                              xytext=(en, y0 + BARH / 2), zorder=4,
                              arrowprops=dict(arrowstyle="-|>", color=col, lw=0.8,
                                              shrinkA=0, shrinkB=0))
        else:
            ax_strip.plot([en, en], [y0 - 0.03, y0 + BARH + 0.03], "-", color="black",
                          lw=0.7, zorder=4)
        span = XHI - XLO
        if (en - st) > span * 0.13 and name:
            ax_strip.text((st + en) / 2, y0 + BARH / 2, name, ha="center", va="center",
                          fontsize=fs - 1.4, color="white", fontweight="bold", zorder=6,
                          path_effects=[pe.withStroke(linewidth=1.1, foreground="0.15")])
            if resp and not open_end:
                ax_strip.text(en, y0 + BARH + 0.05, resp, ha="right", va="bottom",
                              fontsize=fs - 0.8, fontweight="bold", zorder=6)
        else:
            txt = f"{name} {resp}".strip() if not open_end else name
            cx = (st + en) / 2
            if cx > XLO + span * 0.85:   tx, ha = XHI, "right"
            elif cx < XLO + span * 0.15: tx, ha = XLO, "left"
            else:                        tx, ha = cx, "center"
            ax_strip.text(tx, y0 + BARH + 0.05, txt, ha=ha, va="bottom",
                          fontsize=fs - 1.4, color="0.25", zorder=6)

    # --- CA 19-9, linear, optionally broken -------------------------------------------
    def _trace(ax):
        ax.plot(mk.day, mk.v, "-", lw=0.9, color=TRACE, zorder=3)
        plain = mk[mk.cen == ""]
        ax.plot(plain.day, plain.v, "o", ms=2.6, mfc="white", mec=TRACE, mew=0.8,
                ls="none", zorder=4)
        for cen, m in (("above", "^"), ("below", "v")):
            q = mk[mk.cen == cen]
            if len(q):
                ax.plot(q.day, q.v, m, ms=3.6, mfc="none", mec=TRACE, mew=0.9,
                        ls="none", zorder=5)

    if len(mk):
        _trace(ax_lo)
        ax_lo.set_ylim(0, cap)
        ax_lo.axhline(ULN, color="0.55", lw=0.5, ls="--", zorder=1)
        ax_lo.text(XLO, ULN, " ULN 37.5", fontsize=4.4, color="0.45", va="bottom", ha="left")
        ax_lo.tick_params(axis="y", colors=TRACE, labelsize=fs - 0.6)
        if show_ylabel:
            ax_lo.set_ylabel("CA 19-9 (U/mL)", fontsize=fs + 0.6, color=TRACE)

        if n_above:
            _trace(ax_hi)
            top = float(above.v.max())
            ax_hi.set_ylim(top * 0.72, top * 1.18)
            ax_hi.set_yticks([top])
            ax_hi.set_yticklabels([f"{top:,.0f}"])
            ax_hi.tick_params(axis="y", colors=TRACE, labelsize=fs - 0.6)
            ax_hi.tick_params(axis="x", length=0, labelbottom=False)
            ax_hi.spines["bottom"].set_visible(False)
            # Break marks on the left edge only. The house style hides the right spine, so a
            # mark there reads as a stray glyph floating beside the panel.
            kw = dict(marker=[(-1, -0.6), (1, 0.6)], markersize=4, linestyle="none",
                      color="0.35", mec="0.35", mew=0.8, clip_on=False)
            ax_hi.plot([0], [0], transform=ax_hi.transAxes, **kw)
            ax_lo.plot([0], [1], transform=ax_lo.transAxes, **kw)
            ax_lo.spines["top"].set_visible(False)
    else:
        ax_lo.set_ylim(0, 1); ax_lo.set_yticks([])

    # --- RECIST row along the foot of the trace ---------------------------------------
    if len(im):
        span, prev, lvl = XHI - XLO, None, 0
        for _, r in im.iterrows():
            lvl = 0 if (prev is None or (r.day - prev) > span * 0.07) else (lvl + 1) % 3
            ax_lo.text(r.day, 0.015 + 0.052 * lvl, str(r.recist).upper(),
                       transform=ax_lo.get_xaxis_transform(), ha="center", va="bottom",
                       fontsize=fs - 0.9, fontweight="bold", color="0.3")
            prev = r.day

    for ax in filter(None, (ax_strip, ax_hi, ax_lo)):
        ax.axvline(0, color="black", lw=0.7, zorder=2)
    if day_max is not None:
        for ax in filter(None, (ax_strip, ax_hi, ax_lo)):
            ax.axvline(day_max, color="0.45", lw=0.7, ls=(0, (3, 2)), zorder=2)
    ax_lo.set_xlabel("Days relative to biopsy", fontsize=fs + 0.8)

    return dict(ax_strip=ax_strip, ax_hi=ax_hi, ax_lo=ax_lo, n_above=n_above,
                xlim=(XLO, XHI), cap=cap)
