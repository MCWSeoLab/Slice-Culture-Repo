"""
cns_style — publication figure styling for matplotlib/seaborn.

Styling follows the CNS Plots conventions (https://github.com/faridrashidi/cnsplots):
small type, hairline spines, no grid, no top/right spine, journal palettes,
Illustrator-editable vector text.

Two hard rules this module enforces:
  1. Font is Arial (or a metric-identical fallback that is *labelled* Arial in
     vector output, so the file opens with real Arial on a machine that has it).
  2. Every figure is written as PNG *and* PDF *and* SVG by ``save_figure``.

Typical use
-----------
    from cns_style import apply_style, figure, save_figure

    apply_style()                       # or apply_style(palette="Nature")
    fig, ax = figure(width_mm=89, height_mm=70)
    ax.bar(labels, values, color=palette()[0])
    ax.set_ylabel("Tumor volume (mm$^3$)")
    save_figure(fig, "tumor_volume", outdir="figures")
"""

from __future__ import annotations

import re
import shutil
import subprocess
import warnings
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
from cycler import cycler
from matplotlib import font_manager

__all__ = [
    "apply_style",
    "figure",
    "save_figure",
    "palette",
    "PALETTES",
    "SEQUENTIAL",
    "finalize_axes",
    "resolve_font",
    "JOURNAL_WIDTHS_MM",
]

# --------------------------------------------------------------------------
# Palettes (hex values from cnsplots)
# --------------------------------------------------------------------------

PALETTES: dict[str, list[str]] = {
    # Nature Publishing Group — the default; colour-blind tolerable, prints well
    "Nature": ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
               "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85"],
    "Cell": ["#C84C3A", "#2F7E8F", "#E1A22E", "#4E5A8A", "#5F9862",
             "#D07A6A", "#8B6FA8", "#7B8C9E", "#B85F7A", "#6B6B6B"],
    "Science": ["#3B4992", "#EE0000", "#008B45", "#631879", "#008280",
                "#BB0021", "#5F559B", "#A20056", "#808180", "#1B1919"],
    "Ecotyper1": ["#D6372E", "#5189BB", "#70B460", "#985EA8", "#F08F35",
                  "#FADD4B", "#A3A3A3", "#B7D3E5", "#E6D8C2"],
    "Ecotyper2": ["#EB7D5B", "#FED23F", "#B5D33D", "#6CA2EA", "#442288"],
    "Ecotyper3": ["#D13570", "#569AB4", "#70AC58", "#74509D", "#ED7E30",
                  "#F5C945", "#9C5732", "#E787E5"],
    # Okabe–Ito: safe for all common forms of colour vision deficiency
    "OkabeIto": ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                 "#D55E00", "#CC79A7", "#000000"],
    "Grey": ["#3A3A3A", "#7A7A7A", "#ADADAD", "#D4D4D4"],
}

SEQUENTIAL: dict[str, str] = {
    "sequential": "viridis",      # continuous, perceptually uniform
    "sequential_warm": "magma",
    "diverging": "RdBu_r",        # centre on a meaningful zero
    "heatmap": "RdBu_r",
}

# Journal live-area widths, millimetres
JOURNAL_WIDTHS_MM = {
    "nature_single": 89, "nature_double": 183,
    "cell_single": 85, "cell_double": 174,
    "science_single": 55, "science_double": 120, "science_triple": 180,
    "generic_single": 85, "generic_double": 180,
}

MM_PER_INCH = 25.4

# Fonts that are metrically identical to Arial; used only when Arial is absent.
_ARIAL_EQUIVALENTS = ("Arial", "Liberation Sans", "Nimbus Sans", "Helvetica",
                      "Arimo", "Albany AMT")

_FONT_STATE = {"resolved": None, "is_true_arial": False}


def resolve_font(verbose: bool = True) -> str:
    """Return the font family matplotlib will actually render with.

    Prefers real Arial. Falls back to a metric-compatible substitute so layout
    is unchanged; vector output is still *labelled* Arial by ``save_figure``.
    """
    if _FONT_STATE["resolved"]:
        return _FONT_STATE["resolved"]

    # Pick up any font files dropped next to the working directory or in the
    # usual per-user locations (e.g. an Arial.ttf the user supplied).
    for extra in (Path.cwd() / "fonts", Path.home() / ".fonts",
                  Path.home() / "Library" / "Fonts", Path("/Library/Fonts"),
                  Path("/mnt/user-data/uploads")):
        if extra.is_dir():
            for f in list(extra.glob("*.ttf")) + list(extra.glob("*.otf")):
                if "arial" in f.name.lower():
                    try:
                        font_manager.fontManager.addfont(str(f))
                    except Exception:  # pragma: no cover - font file unreadable
                        pass

    available = {f.name for f in font_manager.fontManager.ttflist}
    for name in _ARIAL_EQUIVALENTS:
        if name in available:
            _FONT_STATE["resolved"] = name
            _FONT_STATE["is_true_arial"] = name == "Arial"
            if name != "Arial" and verbose:
                warnings.warn(
                    f"Arial is not installed here; rendering with '{name}', which is "
                    "metrically identical. PDF/SVG text is tagged 'Arial' so the "
                    "files use real Arial on a machine that has it.",
                    stacklevel=2,
                )
            return name

    _FONT_STATE["resolved"] = "DejaVu Sans"
    if verbose:
        warnings.warn(
            "No Arial-compatible font found; falling back to DejaVu Sans. "
            "Metrics will differ slightly. Install fonts-liberation or supply "
            "Arial.ttf in ./fonts/ for exact output.",
            stacklevel=2,
        )
    return "DejaVu Sans"


def palette(name: str = "Nature", n: int | None = None) -> list[str]:
    """Return a qualitative palette as a list of hex strings."""
    if name not in PALETTES:
        raise KeyError(f"Unknown palette {name!r}. Options: {', '.join(PALETTES)}")
    colors = PALETTES[name]
    if n is None:
        return list(colors)
    if n <= len(colors):
        return list(colors[:n])
    reps = -(-n // len(colors))
    return (colors * reps)[:n]


# --------------------------------------------------------------------------
# Style
# --------------------------------------------------------------------------

def apply_style(
    palette_name: str = "Nature",
    base_fontsize: float = 8,
    small_fontsize: float = 7,
    linewidth: float = 0.5,
    colormap: str = "viridis",
    bold_titles: bool = True,
) -> None:
    """Install the CNS Plots rcParams globally. Call once before plotting."""
    font = resolve_font()
    mpl.rcParams.update({
        # ---- type -------------------------------------------------------
        "font.family": "sans-serif",
        "font.sans-serif": [font, "Arial", "Liberation Sans", "Helvetica",
                            "DejaVu Sans"],
        "font.size": base_fontsize,
        # math text uses the same sans face, so superscripts (mm^3) match
        "mathtext.fontset": "custom",
        "mathtext.default": "regular",
        "mathtext.rm": font,
        "mathtext.it": f"{font}:italic",
        "mathtext.bf": f"{font}:bold",
        "mathtext.sf": font,
        "mathtext.cal": font,
        "mathtext.tt": "monospace",
        # ---- vector output: keep text editable in Illustrator -----------
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.01,
        "savefig.dpi": 600,
        # ---- axes -------------------------------------------------------
        "axes.titlesize": base_fontsize,
        "axes.titleweight": "bold" if bold_titles else "normal",
        "axes.titlelocation": "center",
        "axes.titlepad": 4,
        "axes.labelsize": base_fontsize,
        "axes.labelpad": 2,
        "axes.labelcolor": "black",
        "axes.edgecolor": "black",
        "axes.linewidth": linewidth,
        "axes.grid": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.xmargin": 0.05,
        "axes.ymargin": 0.05,
        "axes.prop_cycle": cycler(color=PALETTES[palette_name]),
        "image.cmap": colormap,
        # ---- ticks ------------------------------------------------------
        "xtick.bottom": True, "xtick.top": False, "xtick.color": "black",
        "xtick.direction": "out", "xtick.labelsize": small_fontsize,
        "xtick.major.size": 2, "xtick.major.width": 0.6, "xtick.major.pad": 1,
        "xtick.minor.size": 1, "xtick.minor.width": 0.5,
        "ytick.left": True, "ytick.right": False, "ytick.color": "black",
        "ytick.direction": "out", "ytick.labelsize": small_fontsize,
        "ytick.major.size": 2, "ytick.major.width": 0.6, "ytick.major.pad": 1,
        "ytick.minor.size": 1, "ytick.minor.width": 0.5,
        # ---- legend -----------------------------------------------------
        "legend.fontsize": small_fontsize,
        "legend.title_fontsize": base_fontsize,
        "legend.frameon": False,
        "legend.markerscale": 0.8,
        "legend.handlelength": 0.9,
        "legend.handleheight": 0.7,
        "legend.handletextpad": 0.35,
        "legend.borderpad": 0.2,
        "legend.labelspacing": 0.3,
        # ---- lines / markers -------------------------------------------
        "lines.linewidth": 1.0,
        "lines.markersize": 3,
        "patch.linewidth": 0.5,
        "boxplot.flierprops.markersize": 2,
        "errorbar.capsize": 1.5,
        # ---- figure -----------------------------------------------------
        "figure.dpi": 150,
        "figure.facecolor": "white",
        "figure.titlesize": base_fontsize,
        "figure.titleweight": "bold" if bold_titles else "normal",
    })


def figure(
    width_mm: float | None = None,
    height_mm: float | None = None,
    width_px: float | None = None,
    height_px: float | None = None,
    nrows: int = 1,
    ncols: int = 1,
    **subplot_kw,
):
    """Create a correctly sized figure.

    Size may be given in millimetres (journal-facing) or in CSS pixels at 72 dpi
    (the cnsplots convention). Defaults to a Nature single column, 89 x 65 mm.
    """
    if width_px is not None:
        width_mm = width_px / 72 * MM_PER_INCH
    if height_px is not None:
        height_mm = height_px / 72 * MM_PER_INCH
    width_mm = 89 if width_mm is None else width_mm
    height_mm = 65 if height_mm is None else height_mm
    figsize = (width_mm / MM_PER_INCH, height_mm / MM_PER_INCH)
    fig, ax = plt.subplots(nrows, ncols, figsize=figsize, **subplot_kw)
    return fig, ax


def finalize_axes(ax, despine: bool = True, tight: bool = True) -> None:
    """Apply the last-mile CNS conventions to one axes."""
    axes = ax if isinstance(ax, Iterable) else [ax]
    for a in getattr(axes, "flat", axes):
        if despine:
            a.spines["top"].set_visible(False)
            a.spines["right"].set_visible(False)
        a.tick_params(which="both", direction="out")
    if tight:
        plt.tight_layout(pad=0.3)


def add_panel_labels(fig, axes: Sequence, labels: Sequence[str] | None = None,
                     fontsize: float = 9, x: float = -0.02, y: float = 1.06) -> None:
    """Add bold A, B, C... panel labels in journal position."""
    letters = labels or [chr(ord("A") + i) for i in range(len(axes))]
    for a, letter in zip(axes, letters):
        a.text(x, y, letter, transform=a.transAxes, fontsize=fontsize,
               fontweight="bold", va="bottom", ha="right")


# --------------------------------------------------------------------------
# Export — always PNG + PDF + SVG
# --------------------------------------------------------------------------

_FONT_FAMILY_RE = re.compile(r"font-family:\s*[^;}\n<]+")


def _force_arial_in_svg(path: Path) -> None:
    """Rewrite font-family declarations to Arial in an SVG written with
    svg.fonttype='none'. Text stays live and picks up real Arial downstream."""
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return
    text = _FONT_FAMILY_RE.sub("font-family: Arial", text)
    text = re.sub(r'font-family="[^"]*"', 'font-family="Arial"', text)
    path.write_text(text, encoding="utf-8")


def save_figure(
    fig=None,
    stem: str = "figure",
    outdir: str | Path = ".",
    formats: Sequence[str] = ("png", "pdf", "svg"),
    dpi: int = 600,
    transparent_vector: bool = True,
    png_facecolor: str = "white",
    verbose: bool = True,
) -> dict[str, Path]:
    """Save one figure to PNG, PDF and SVG. Returns {format: path}.

    PNG is written on an opaque background (safe for slides and Word); PDF and
    SVG are written transparent with live, Arial-tagged text for Illustrator.
    """
    fig = fig or plt.gcf()
    outdir = Path(outdir).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)
    stem = Path(stem).stem
    written: dict[str, Path] = {}

    for fmt in formats:
        fmt = fmt.lower().lstrip(".")
        path = outdir / f"{stem}.{fmt}"
        if fmt == "png":
            fig.savefig(path, format="png", dpi=dpi, bbox_inches="tight",
                        pad_inches=0.01, transparent=False,
                        facecolor=png_facecolor)
        else:
            fig.savefig(path, format=fmt, bbox_inches="tight", pad_inches=0.01,
                        transparent=transparent_vector)
        if fmt == "svg":
            _force_arial_in_svg(path)
        written[fmt] = path

    if verbose:
        for fmt, path in written.items():
            print(f"  {fmt.upper():<4} {path}")
    return written


def check_font_embedding(pdf_path: str | Path) -> list[str]:
    """List fonts embedded in a PDF (needs pdffonts from poppler-utils).
    Useful when a journal rejects a figure for font problems."""
    if not shutil.which("pdffonts"):
        return ["pdffonts not installed (apt-get install -y poppler-utils)"]
    out = subprocess.run(["pdffonts", str(pdf_path)], capture_output=True,
                         text=True, check=False)
    return out.stdout.splitlines()
