# -*- coding: utf-8 -*-
"""Rebuild Figure 4's panels B and F larger, in place.

Same magnification (250 um fields), same tiles, same two agents per case - just more of the
page.  The tiles go from 17.3 mm to about 20.6 mm by:

  * filling the white patch completely (the v3 panels were 238 x 116 pt inside a 254 x 122 pt
    patch), and extending the patch into the blank bands above and below it, which a
    row-by-row ink scan of the rendered page put at 40.6-45.0 mm and 83.0-89.1 mm from the
    top edge;
  * cutting the row-label gutter from ~11 mm to 5 mm by setting the agent labels vertically.

Nothing else on the page moves - panels A, C, D, E, G and H are untouched.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

PT = 72.0 / 25.4                      # points per mm
DPI = 600

# --- page and patch geometry, in PDF points, measured from fig4_two_cases_v3.pdf
PAGE_W, PAGE_H = 547.936, 643.3269
PATCH_TOP_MM, PATCH_BOT_MM = 41.2, 88.6          # from the top of the page
Y_TOP = PAGE_H - PATCH_TOP_MM * PT
Y_BOT = PAGE_H - PATCH_BOT_MM * PT
PATCHES = {                                       # x0, x1 in points
    "B": (4.0, 258.0),
    "F": (262.0, 520.0),
}

STAINS = ["H&E", "Cleaved caspase-3", "CD45", "EpCAM"]
PANELS = {
    "B": [("Gem/nab-Pac\n(1%)", ["_tiles_B/10_27_LiverBX_GemPac_1.png",
                                 "_tiles_B/L_Bx_10_27_GemNab_CC_2.png",
                                 "_tiles_B/10_27_LiverBx_GemPac_CD45_4.png",
                                 "_tiles_B/10_27_LiverBx_GemPac_EpCam_5.png"]),
          ("FOLFIRINOX\n(61%)", ["_tiles_B/10_27_LiverBx_Folfirinox_1.png",
                                 "_tiles_B/L_BX_10_27_Folfi_CC_2.png",
                                 "_tiles_B/10_27_LiverBx_Folfirinox_CD45_4.png",
                                 "_tiles_B/10_27_LiverBx_Folfirinox_EpCam_5.png"])],
    "F": [("Vemurafenib\n(91%)", ["_tiles_F/L_Bx_A3_Raf_1.png", "_tiles_F/L_Bx_A3_Raf_2.png",
                                  "_tiles_F/L_Bx_A3_Raf_4.png", "_tiles_F/L_Bx_A3_Raf_5.png"]),
          ("Trametinib\n(19%)", ["_tiles_F/L_Bx_A2_Mek_1.png", "_tiles_F/L_Bx_A2_Mek_2.png",
                                 "_tiles_F/L_Bx_A2_Mek_4.png", "_tiles_F/L_Bx_A2_Mek_5.png"])],
}

FIELD_UM, BAR_UM = 250.0, 50.0
GUTTER_MM = 5.2                       # rotated agent labels live here
TITLE_MM = 3.6
GAP_MM = 0.7

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "pdf.fonttype": 42, "svg.fonttype": "none",
})


def build(key):
    x0, x1 = PATCHES[key]
    w_mm = (x1 - x0) / PT
    h_mm = (Y_TOP - Y_BOT) / PT
    rows = PANELS[key]

    tile_w = min((w_mm - GUTTER_MM - 3 * GAP_MM) / 4,
                 (h_mm - TITLE_MM - GAP_MM) / 2)
    block_w = 4 * tile_w + 3 * GAP_MM
    x_left = GUTTER_MM
    y_top_block = h_mm - TITLE_MM

    fig = plt.figure(figsize=(w_mm / 25.4, h_mm / 25.4), dpi=DPI)
    fig.patch.set_alpha(0.0)

    for r, (label, files) in enumerate(rows):
        y = y_top_block - (r + 1) * tile_w - r * GAP_MM
        for c, f in enumerate(files):
            x = x_left + c * (tile_w + GAP_MM)
            ax = fig.add_axes([x / w_mm, y / h_mm, tile_w / w_mm, tile_w / h_mm])
            ax.imshow(mpimg.imread(f), interpolation="lanczos")
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(True); sp.set_linewidth(0.5); sp.set_color("0.15")
            if r == 0:
                ax.set_title(STAINS[c], fontsize=5.2, pad=1.8)
            if r == len(rows) - 1 and c == len(files) - 1:
                # 50 um bar on a 250 um field = 20.0% of the tile width, matching the v3 bar
                w_px = mpimg.imread(f).shape[1]
                bx1 = w_px * 0.955
                bx0 = bx1 - w_px * (BAR_UM / FIELD_UM)
                by = mpimg.imread(f).shape[0] * 0.93
                ax.plot([bx0, bx1], [by, by], color="black", lw=1.3,
                        solid_capstyle="butt", zorder=6)
                ax.text((bx0 + bx1) / 2, by * 0.985, f"{BAR_UM:.0f} µm", color="black",
                        fontsize=4.6, ha="center", va="bottom", zorder=6)
        fig.text((GUTTER_MM / 2) / w_mm, (y + tile_w / 2) / h_mm, label,
                 fontsize=4.8, rotation=90, va="center", ha="center", linespacing=1.10)

    # The panel letter lived inside the old overlay, so it is whited out with it and has to
    # be redrawn here.  Page-level letters sit ~3 pt inside the patch; match that.
    letter_x_pt = 7.1 if key == "B" else 269.9
    fig.text(((letter_x_pt - x0) / PT) / w_mm, 1.0, key, fontsize=9, fontweight="bold",
             va="top", ha="left")

    out = Path(f"panel{key}_big.pdf")
    fig.savefig(out, dpi=DPI, transparent=True)
    plt.close(fig)
    print(f"panel {key}: {w_mm:.1f} x {h_mm:.1f} mm, tile {tile_w:.2f} mm "
          f"(was 17.27), block {block_w:.1f} mm -> {out}")
    return tile_w


if __name__ == "__main__":
    for k in ("B", "F"):
        build(k)
