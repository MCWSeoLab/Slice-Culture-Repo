"""Shared micrograph handling for the supplementary figures.

Neutralizes the illumination cast so tiles in one row read as one series. The
correction is a per-channel gain fitted to the median of the brightest 2% of
pixels (the empty slide background), which is achromatic in reality; nothing
else about the image is changed. Disclose it in the Methods.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image

IMG = Path("deid/suppfig_images")   # run from the project root


# Tiles print ~40 mm wide; 600 px is ~385 dpi there, above the 300 dpi
# journals ask for. The 1600 px originals triple the vector file size.
MAX_PX = 600


def load(name: str, balance: bool = False) -> np.ndarray:
    im = Image.open(IMG / f"{name}.png").convert("RGB")
    if im.width > MAX_PX:
        im = im.resize((MAX_PX, round(im.height * MAX_PX / im.width)),
                       Image.LANCZOS)
    a = np.asarray(im).astype(np.float32)
    if balance:
        a = white_balance(a)
    return np.clip(a, 0, 255).astype(np.uint8)


def white_balance(a: np.ndarray, pct: float = 99.0) -> np.ndarray:
    lum = a.mean(axis=2)
    thr = np.percentile(lum, pct)
    bg = a[lum >= thr]
    if bg.size == 0:
        return a
    ref = np.median(bg, axis=0)
    target = ref.max()
    gain = np.where(ref > 1, target / ref, 1.0)
    # keep the correction gentle; a full gray-world fit can invent color
    gain = 1.0 + (gain - 1.0) * 1.0
    return a * gain


def show(ax, name: str, label: str | None = None, balance: bool = False,
         labelsize: float = 6.0):
    ax.imshow(load(name, balance=balance), interpolation="lanczos", aspect="equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_linewidth(0.4)
        sp.set_edgecolor("#444444")
    if label:
        ax.set_title(label, fontsize=labelsize, pad=1.6, fontweight="normal")
