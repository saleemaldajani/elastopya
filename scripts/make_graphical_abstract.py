#!/usr/bin/env python3
"""Graphical abstract for the Scripta Materialia submission (2026-11-01).

Elsevier: one panel, at least 1328 x 531 px (w x h), readable at 13 x 5 cm. It is drawn at exactly 13 x 5 cm, so the
lettering sizes are the printed sizes, and exported as PDF (vector), PNG and TIFF (600 dpi, 3071 x 1181 px).

Left to right it follows the main-text figures, with the same data:
  1  unstable window: demixing free energy (PWSCF-PAW model of Fig. 1), shaded where its curvature is positive
  2  spinodal decomposition of Fe-36Cr in MD, random alloy -> fully separated (Figs. 1c, 3c)
  3  SAW frequency of the random and fully decomposed MD alloys at twelve compositions (Figs. 2a, 3)
  4  the result: shear stiffening of alpha', higher SAW frequency, lower K/G (Sec. 3)

Output: manuscript/figures/graphical_abstract/graphical_abstract.{pdf,png,tiff}
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_figures as mf  # noqa: E402
from visualize_decomposition import CR, FE, _cube_ax, md_grid3d, series  # noqa: E402

OUT = mf.ROOT / "manuscript" / "figures" / "graphical_abstract"
W, H = 13 / 2.54, 5 / 2.54                       # inches
BLUE, VERM, GREEN, GREY = "#0072B2", "#D55E00", "#009E73", "0.35"
WIN = "#E8E2F2"                                  # spinodal-window shading (light lilac, prints as light grey)
FS, FS_T = 7.0, 7.5                              # lettering, panel titles (pt)


def box(fig, x0, y0, w, h):
    """Axes from a position in inches."""
    return fig.add_axes([x0 / W, y0 / H, w / W, h / H])


def arrow(fig, x0, x1, y):
    fig.patches.append(FancyArrowPatch((x0 / W, y / H), (x1 / W, y / H), transform=fig.transFigure,
                                       arrowstyle="-|>", mutation_scale=9, lw=1.3, color=GREY))


def title(fig, x, y, num, text):
    fig.text(x / W, y / H, num, fontsize=FS_T + 0.5, fontweight="bold", color="white", ha="center", va="center",
             bbox=dict(boxstyle="circle,pad=0.18", fc=GREY, ec="none"))
    fig.text((x + 0.11) / W, y / H, text, fontsize=FS_T, fontweight="bold", ha="left", va="center")


def frame(ax):
    ax.tick_params(direction="in", length=2.2, width=0.6, labelsize=FS - 0.5, pad=1.5)
    for s in ax.spines.values():
        s.set_linewidth(0.6)


def shade_window(ax, ylo, yhi):
    p = mf.read_csv("sheet_levesque_parameters.csv")
    a0, b0, c0, d0 = (p[k][0] for k in ("alpha", "beta", "gamma", "delta"))
    xs = np.linspace(0.005, 0.995, 2000)
    _, curv = mf.levesque_model(xs, a0, b0, c0, d0)
    ax.fill_between(xs, ylo, yhi, where=curv > 0, color=WIN, lw=0, zorder=0)
    return xs, a0, b0, c0, d0


def main():
    plt.rcParams.update({"font.family": "Arial", "mathtext.fontset": "custom", "mathtext.rm": "Arial",
                         "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "font.size": FS,
                         "pdf.fonttype": 42, "axes.linewidth": 0.6})
    fig = plt.figure(figsize=(W, H))
    yb, hb, yt = 0.34, 1.18, 1.76                 # axes bottom, height; title row

    # 1 unstable window (Fig. 1)
    a = box(fig, 0.36, yb, 0.90, hb)
    xs, a0, b0, c0, d0 = shade_window(a, -0.09, 0.02)
    g, _ = mf.levesque_model(xs, a0, b0, c0, d0)
    a.plot(xs, g, color=GREY, lw=1.1, zorder=2)
    a.axhline(0, color="0.6", lw=0.5, zorder=1)
    a.plot([mf.X_FE36CR], [-0.083], "*", color=BLUE, ms=7, zorder=4, clip_on=False)
    a.text(0.485, 0.012, "unstable", fontsize=FS - 0.5, ha="center", va="top", color="0.25")
    a.text(mf.X_FE36CR + 0.07, -0.083, "Fe\u201336Cr", fontsize=FS - 1, color=BLUE, ha="left", va="center")
    a.text(0.53, -0.022, "DFT", fontsize=FS - 0.5, color=GREY, ha="left", va="center")
    a.set(xlim=(0, 1), ylim=(-0.09, 0.02), xticks=[0, 0.5, 1], yticks=[-0.06, 0])
    a.set_xticklabels(["0", "0.5", "1"]); a.set_yticklabels(["−0.06", "0"])
    a.set_xlabel(r"$x_{\rm Cr}$", labelpad=0.5)
    a.set_ylabel(r"$\Delta G_{\rm demix}$ (eV/atom)", labelpad=0.5, fontsize=FS - 0.5)
    frame(a)
    title(fig, 0.10, yt, "1", "Unstable window")

    # 2 spinodal decomposition of Fe-36Cr (MD)
    f36 = series()[0.36]
    for i, (f, lab) in enumerate(((f36[0], "random"), (f36[-1], r"$\alpha$ + $\alpha'$"))):
        c = box(fig, 1.43 + i * 0.54, yb + 0.16, 0.42, 0.80)
        _cube_ax(c, md_grid3d(f))
        fig.text((1.43 + i * 0.54 + 0.21) / W, (yb - 0.02) / H, lab, fontsize=FS, ha="center", va="bottom")
    arrow(fig, 1.86, 1.96, yb + 0.56)
    for k, (name, col) in enumerate((("Fe", FE), ("Cr", CR))):
        x0 = 1.62 + k * 0.38
        fig.patches.append(plt.Rectangle((x0 / W, (yb + 1.13) / H), 0.07 / W, 0.07 / H, transform=fig.transFigure,
                                         fc=col, ec="none"))
        fig.text((x0 + 0.10) / W, (yb + 1.165) / H, name, fontsize=FS - 0.5, va="center")
    title(fig, 1.50, yt, "2", "Decomposition")

    # 3 SAW frequency, random vs fully decomposed (Figs. 2a, 3)
    s = box(fig, 2.86, yb, 0.90, hb)
    shade_window(s, 0.43, 0.54)
    dec = mf.md_decomposition()
    xs3 = np.array(sorted(dec))
    f0 = np.array([dec[x]["f"][np.argmin(dec[x]["phi"])] for x in xs3])
    f1 = np.array([dec[x]["f"][np.argmax(dec[x]["phi"])] for x in xs3])
    s.plot(xs3, f0, "-s", color=GREEN, ms=2.6, lw=0.9, zorder=3, label="random")
    s.plot(xs3, f1, "--s", color=GREEN, mfc="white", mew=0.8, ms=2.8, lw=0.9, zorder=3, label="decomposed")
    i = int(np.argmax(f1 / f0))
    s.annotate("", xy=(xs3[i], f1[i] - 0.002), xytext=(xs3[i], f0[i] + 0.002),
               arrowprops=dict(arrowstyle="-|>", lw=0.9, color=VERM, mutation_scale=7), zorder=4)
    s.text(xs3[i] + 0.03, 0.4425, f"+{100 * (f1[i] / f0[i] - 1):.1f}%", color=VERM, fontsize=FS,
           fontweight="bold", va="center", ha="left")
    s.plot([mf.X_FE36CR], [0.4345], "*", color=BLUE, ms=7, zorder=4, clip_on=False)
    s.legend(loc="upper left", fontsize=FS - 1, frameon=False, handlelength=1.6, borderaxespad=0.2,
             labelspacing=0.15, handletextpad=0.4)
    s.set(xlim=(0, 1), ylim=(0.43, 0.54), xticks=[0, 0.5, 1], yticks=[0.45, 0.50])
    s.set_xticklabels(["0", "0.5", "1"])
    s.set_xlabel(r"$x_{\rm Cr}$", labelpad=0.5)
    s.set_ylabel(r"$f_{\rm SAW}$ (GHz)", labelpad=0.5, fontsize=FS - 0.5)
    frame(s)
    title(fig, 2.74, yt, "3", "Stiffer SAW")

    # 4 result
    x4, w4 = 3.96, 1.08
    fig.patches.append(FancyBboxPatch((x4 / W, (yb + 0.02) / H), w4 / W, (hb - 0.04) / H, transform=fig.transFigure,
                                      boxstyle="round,pad=0.0,rounding_size=0.012", fc="#F4F4F4", ec="0.6", lw=0.6))
    lines = [(r"$C'$ $\uparrow$", "shear stiffening\nof Cr-rich $\\alpha'$"),
             (r"$f_{\rm SAW}$ $\uparrow$", "non-destructive\nsignal"),
             (r"$K/G$ $\downarrow$", "embrittlement\n(Pugh)")]
    for k, (big, small) in enumerate(lines):
        y = yb + hb - 0.20 - k * 0.37
        fig.text((x4 + 0.06) / W, y / H, big, fontsize=FS + 0.5, fontweight="bold", color=VERM, ha="left", va="center")
        fig.text((x4 + 0.45) / W, y / H, small, fontsize=FS - 1, ha="left", va="center", linespacing=1.0)
    for k in range(2):
        y = yb + hb - 0.20 - k * 0.37
        fig.text((x4 + 0.17) / W, (y - 0.185) / H, "↓", fontsize=FS - 0.5, color=GREY, ha="center", va="center")
    title(fig, 4.02, yt, "4", "Indicator")

    for x0, x1 in ((1.29, 1.41), (2.42, 2.54), (3.79, 3.93)):
        arrow(fig, x0, x1, yb + hb / 2)

    OUT.mkdir(parents=True, exist_ok=True)
    meta = {"Title": "Graphical abstract", "Creator": "S. A. Al Dajani"}
    fig.savefig(OUT / "graphical_abstract.pdf", metadata=meta)
    fig.savefig(OUT / "graphical_abstract.png", dpi=600)
    fig.savefig(OUT / "graphical_abstract.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(OUT / "graphical_abstract_screen.png", dpi=96)    # what readers see on the journal page (13 x 5 cm)
    print("max stiffening:", xs3[i], round(100 * (f1[i] / f0[i] - 1), 2), "%")


if __name__ == "__main__":
    main()
