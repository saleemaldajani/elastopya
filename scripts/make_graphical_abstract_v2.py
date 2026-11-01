#!/usr/bin/env python3
"""Graphical abstract for the Scripta Materialia submission (2026-11-01).

Elsevier: one panel, at least 1328 x 531 px (w x h), readable at 13 x 5 cm. It is drawn at exactly 13 x 5 cm, so the
lettering sizes are the printed sizes, and exported as PDF (vector), PNG and TIFF (600 dpi, 3071 x 1181 px).

Left to right it follows the main-text figures, with the same data (v2, 2026-11-01: DFT in panels 1 and 3):
  1  unstable window: QE demixing free energy (Fig. 1a), shaded where the PWSCF-PAW model curvature is positive (Fig. 1b)
  2  spinodal decomposition of Fe-36Cr in MD, random alloy -> fully separated (Figs. 1c, 3c)
  3  DFT elasticity: QE tetragonal shear modulus C' (Fig. 2a tensors; Sec. 3); partitioning of Fe-36Cr into alpha/alpha'
  4  SAW frequency of the random and fully decomposed MD alloys at twelve compositions (Figs. 2-3)
  banner: the result (Sec. 3): C' up, f_SAW up, K/G down. v1 (no DFT elasticity panel): make_graphical_abstract_v1.py

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
    yb, hb, yt = 0.56, 0.98, 1.80                 # axes bottom, height; title row
    pw = 0.80                                     # plot width

    # 1 unstable window (Fig. 1): QE demixing free energy
    a = box(fig, 0.33, yb, pw, hb)
    shade_window(a, -0.08, 0.02)
    q = mf.qe_mixing()
    a.axhline(0, color="0.6", lw=0.5, zorder=1)
    a.plot(q["x"], q["demix"], "-o", color=BLUE, mfc="white", mew=0.9, ms=2.8, lw=0.8, zorder=3)
    a.plot([mf.X_FE36CR], [-0.074], "*", color=BLUE, ms=6.5, zorder=4)
    a.text(0.5, 0.013, "unstable", fontsize=FS - 1, ha="center", va="top", color="0.25")
    a.text(0.04, -0.045, "DFT\n(QE)", fontsize=FS - 1, color=BLUE, ha="left", va="center", linespacing=1.0)
    a.set(xlim=(0, 1), ylim=(-0.08, 0.02), xticks=[0, 0.5, 1], yticks=[-0.06, 0])
    a.set_xticklabels(["0", "0.5", "1"]); a.set_yticklabels(["\u22120.06", "0"])
    a.set_xlabel(r"$x_{\rm Cr}$", labelpad=0.5)
    a.set_ylabel(r"$\Delta G_{\rm demix}$ (eV/atom)", labelpad=0.5, fontsize=FS - 1)
    frame(a)
    title(fig, 0.08, yt, "1", "Unstable window")

    # 2 spinodal decomposition of Fe-36Cr (MD)
    f36 = series()[0.36]
    x2 = 1.30
    for i, (f, lab) in enumerate(((f36[0], "random"), (f36[-1], r"$\alpha$ + $\alpha'$"))):
        c = box(fig, x2 + i * 0.47, yb + 0.14, 0.34, 0.66)
        _cube_ax(c, md_grid3d(f))
        fig.text((x2 + i * 0.47 + 0.17) / W, (yb - 0.01) / H, lab, fontsize=FS - 0.5, ha="center", va="bottom")
    arrow(fig, x2 + 0.355, x2 + 0.455, yb + 0.47)
    for k, (name, col) in enumerate((("Fe", FE), ("Cr", CR))):
        x0 = x2 + 0.18 + k * 0.34
        fig.patches.append(plt.Rectangle((x0 / W, (yb + 0.93) / H), 0.07 / W, 0.07 / H, transform=fig.transFigure,
                                         fc=col, ec="none"))
        fig.text((x0 + 0.10) / W, (yb + 0.965) / H, name, fontsize=FS - 1, va="center")
    title(fig, 1.28, yt, "2", "Decomposition")

    # 3 DFT elasticity: C' of the alloy, and where decomposition takes Fe-36Cr
    e = box(fig, 2.60, yb, pw, hb)
    shade_window(e, 25, 175)
    for d, kw in ((mf.qe_dft(), dict(fmt="-o", mfc="white", ms=2.8)), (mf.qe_dft_aexp(), dict(fmt="-D", ms=2.4))):
        e.plot(d["x"], (d["c11"] - d["c12"]) / 2, kw["fmt"], color=BLUE, mfc=kw.get("mfc", BLUE), mew=0.9,
               ms=kw["ms"], lw=0.8, zorder=3)
    ya = 38
    e.annotate("", xy=(0.15, ya), xytext=(mf.X_FE36CR, ya),
               arrowprops=dict(arrowstyle="-|>", lw=1.0, color=VERM, mutation_scale=7, shrinkA=4))
    e.annotate("", xy=(0.85, ya), xytext=(mf.X_FE36CR, ya),
               arrowprops=dict(arrowstyle="-|>", lw=1.0, color=VERM, mutation_scale=7, shrinkA=4))
    e.plot([mf.X_FE36CR], [ya], "*", color=BLUE, ms=6.5, zorder=4)
    e.text(0.12, ya, r"$\alpha$", color=VERM, fontsize=FS, ha="right", va="center")
    e.text(0.87, ya, r"$\alpha'$", color=VERM, fontsize=FS, ha="left", va="center")
    e.text(0.05, 160, "DFT (QE)", fontsize=FS - 1, color=BLUE, ha="left", va="center")
    e.set(xlim=(0, 1), ylim=(25, 175), xticks=[0, 0.5, 1], yticks=[50, 100, 150])
    e.set_xticklabels(["0", "0.5", "1"])
    e.set_xlabel(r"$x_{\rm Cr}$", labelpad=0.5)
    e.set_ylabel(r"Shear modulus $C'$ (GPa)", labelpad=0.5, fontsize=FS - 1)
    frame(e)
    title(fig, 2.40, yt, "3", "DFT: stiff " + r"$\alpha'$")

    # 4 SAW frequency, random vs fully decomposed (Figs. 2-3)
    s = box(fig, 3.98, yb, pw + 0.10, hb)
    shade_window(s, 0.43, 0.545)
    dec = mf.md_decomposition()
    xs3 = np.array(sorted(dec))
    f0 = np.array([dec[x]["f"][np.argmin(dec[x]["phi"])] for x in xs3])
    f1 = np.array([dec[x]["f"][np.argmax(dec[x]["phi"])] for x in xs3])
    s.plot(xs3, f0, "-s", color=GREEN, ms=2.4, lw=0.8, zorder=3, label="random")
    s.plot(xs3, f1, "--s", color=GREEN, mfc="white", mew=0.8, ms=2.6, lw=0.8, zorder=3, label="decomposed")
    i = int(np.argmax(f1 / f0))
    s.annotate("", xy=(xs3[i], f1[i] - 0.002), xytext=(xs3[i], f0[i] + 0.002),
               arrowprops=dict(arrowstyle="-|>", lw=0.9, color=VERM, mutation_scale=7), zorder=4)
    s.text(xs3[i] + 0.03, 0.4425, f"+{100 * (f1[i] / f0[i] - 1):.1f}%", color=VERM, fontsize=FS,
           fontweight="bold", va="center", ha="left")
    s.plot([mf.X_FE36CR], [0.4355], "*", color=BLUE, ms=6.5, zorder=4)
    s.legend(loc="upper left", fontsize=FS - 1, frameon=False, handlelength=1.5, borderaxespad=0.2,
             labelspacing=0.1, handletextpad=0.3)
    s.set(xlim=(0, 1), ylim=(0.43, 0.545), xticks=[0, 0.5, 1], yticks=[0.45, 0.50])
    s.set_xticklabels(["0", "0.5", "1"])
    s.set_xlabel(r"$x_{\rm Cr}$", labelpad=0.5)
    s.set_ylabel(r"$f_{\rm SAW}$ (GHz)", labelpad=0.5, fontsize=FS - 1)
    frame(s)
    title(fig, 3.74, yt, "4", "Stiffer SAW")

    for x0, x1 in ((1.15, 1.27), (2.18, 2.30), (3.44, 3.56)):
        arrow(fig, x0, x1, yb + hb / 2)

    # banner: the result
    fig.patches.append(FancyBboxPatch((0.08 / W, 0.04 / H), (W - 0.16) / W, 0.20 / H, transform=fig.transFigure,
                                      boxstyle="round,pad=0.0,rounding_size=0.01", fc="#F4F4F4", ec="0.65", lw=0.5))
    fig.text(0.5, 0.14 / H, r"Spinodal decomposition:  $C'$ $\uparrow$  $\rightarrow$  $f_{\rm SAW}$ $\uparrow$,  "
             r"$K/G$ $\downarrow$   =   a non-destructive indicator of embrittlement", fontsize=FS, ha="center",
             va="center", color="0.1")

    OUT.mkdir(parents=True, exist_ok=True)
    meta = {"Title": "Graphical abstract", "Creator": "S. A. Al Dajani"}
    fig.savefig(OUT / "graphical_abstract.pdf", metadata=meta)
    fig.savefig(OUT / "graphical_abstract.png", dpi=600)
    fig.savefig(OUT / "graphical_abstract.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(OUT / "graphical_abstract_screen.png", dpi=96)    # what readers see on the journal page (13 x 5 cm)
    print("max stiffening:", xs3[i], round(100 * (f1[i] / f0[i] - 1), 2), "%")


if __name__ == "__main__":
    main()
