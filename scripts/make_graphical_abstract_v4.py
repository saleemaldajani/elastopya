#!/usr/bin/env python3
"""Graphical abstract for the Scripta Materialia submission (v4, 2026-11-01).

Elsevier: one panel, at least 1328 x 531 px (w x h), readable at 13 x 5 cm. Drawn at exactly 13 x 5 cm, so the lettering
sizes are the printed sizes; exported as PDF (vector), PNG and TIFF (600 dpi, 3071 x 1181 px).

Two halves, each a key result plotted above the structures it comes from (user, 2026-11-01: "show the cube trends for DFT
and MD with the key points in figures above"):
  MD   Fe-36Cr decomposing: 25^3 cells at phi = 0.5 ... 1.0, each under its phi tick, and their Rayleigh SAW frequency
       (Fig. 3; mean of five alloys)
  DFT  QE 40-atom SQS cells (relaxed, local moments as arrows) for alpha (5% Cr), Fe-36Cr (35%) and alpha' (95%), with C'
       at the experimental lattice parameter, joined to their compositions on a plot of the highest (P)SAW frequency over
       10^4 grain orientations: only the Cr-rich alpha' tensor enters the range unallowed for pristine CASS, and no MD
       tensor does (Fig. 2b); spinodal window of Fig. 1b shaded
Earlier versions: make_graphical_abstract_v1.py ... _v3.py (outputs in figures/graphical_abstract_v1 ... _v3).

Output: manuscript/figures/graphical_abstract/graphical_abstract.{pdf,png,tiff}, graphical_abstract_screen.png (96 dpi)
"""
import glob
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import ConnectionPatch, FancyArrowPatch  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_figures as mf  # noqa: E402
import visualize_sqs as VS  # noqa: E402
from visualize_decomposition import _cube_ax, md_grid3d, series  # noqa: E402

OUT = mf.ROOT / "manuscript" / "figures" / "graphical_abstract"
W, H = 13 / 2.54, 5 / 2.54                       # inches
BLUE, VERM, GREEN = "#0072B2", "#D55E00", "#009E73"
WIN = "#ECE6F5"                                  # spinodal-window shading
FS = 7.0                                         # lettering (pt)
LO, HI = 0.609, 0.833                            # unallowed range for pristine CASS at 6.4 um (Fig. 2b)


def box(fig, x0, y0, w, h):
    return fig.add_axes([x0 / W, y0 / H, w / W, h / H])


def frame(ax):
    ax.tick_params(direction="in", length=2.2, width=0.6, labelsize=FS - 0.5, pad=1.5)
    for s in ax.spines.values():
        s.set_linewidth(0.6)


def header(fig, x0, x1, y, text, col):
    fig.text(x0 / W, y / H, text, fontsize=FS + 0.5, fontweight="bold", color=col, ha="left", va="center")
    fig.add_artist(plt.Line2D([x0 / W, x1 / W], [(y - 0.11) / H] * 2, color=col, lw=1.1))


def top_mode_ghz(path):
    """Highest (P)SAW frequency over the 10^4 grain orientations, with the intensity cut of Fig. 2b."""
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    keep = d[:, 3] / d[:, 2] >= 2e-3
    return float(np.r_[d[:, 0], d[keep, 1]].max() / 6.4e-6 / 1e9)


def series_top(prefix):
    rows = sorted((float(re.search(r"_x([\d.]+)\.csv", p).group(1)), top_mode_ghz(p))
                  for p in glob.glob(str(mf.DATA / "saw" / f"{prefix}_x*.csv")))
    return np.array([r[0] for r in rows]), np.array([r[1] for r in rows])


def main():
    plt.rcParams.update({"font.family": "Arial", "mathtext.fontset": "custom", "mathtext.rm": "Arial",
                         "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "font.size": FS,
                         "pdf.fonttype": 42, "axes.linewidth": 0.6})
    fig = plt.figure(figsize=(W, H))
    yp, hp = 1.00, 0.62                           # plots: bottom, height (inches)
    yt = 1.86                                     # header row

    # ---------------- MD: Fe-36Cr decomposing, frequency above its cells
    x0, w = 0.48, 1.98
    a = box(fig, x0, yp, w, hp)
    d = mf.md_decomposition()[0.36]
    o = np.argsort(d["phi"]); ph, f = d["phi"][o], d["f"][o]
    a.plot(ph, f, "-s", color=GREEN, ms=3.0, lw=1.1, zorder=3)
    a.plot([0.5, 1.0], [f[0], f[0]], ls=":", color="0.6", lw=0.7, zorder=1)
    a.annotate("", xy=(1.0, f[-1] - 0.0006), xytext=(1.0, f[0]),
               arrowprops=dict(arrowstyle="-|>", lw=1.0, color=VERM, mutation_scale=7, shrinkA=0, shrinkB=0))
    a.text(0.975, (f[0] + f[-1]) / 2, f"+{100 * (f[-1] / f[0] - 1):.1f}%", color=VERM, fontsize=FS + 0.5,
           fontweight="bold", ha="right", va="center")
    a.set(xlim=(0.45, 1.05), ylim=(0.446, 0.478), xticks=ph, yticks=[0.45, 0.46, 0.47])
    a.set_xticklabels([f"{v:.1f}" for v in ph])
    a.set_ylabel(r"$f_{\rm SAW}$ (GHz)", fontsize=FS - 0.5, labelpad=1)
    a.text(1.05, 0.4425, r"$\phi$", fontsize=FS, ha="left", va="top", clip_on=False)
    frame(a)
    cw, ch = 0.29, 0.33                           # cells under their phi ticks
    for fpath, p in zip(series()[0.36], ph):
        xc = x0 + (p - 0.45) / 0.60 * w
        _cube_ax(box(fig, xc - cw / 2, 0.36, cw, ch), md_grid3d(fpath))
    fig.text((x0 + w / 2) / W, 0.31 / H, r"Fe$_{64}$Cr$_{36}$, MD 25$^3$ cells:  random $\rightarrow$ $\alpha$ + $\alpha'$",
             fontsize=FS - 0.5, ha="center", va="top", color="0.2")
    header(fig, 0.10, 2.50, yt, "MD: decomposition stiffens the SAW", GREEN)

    # ---------------- DFT: alpha, Fe-36Cr, alpha' cells under a plot of the top surface mode
    x1, w1 = 3.06, 1.84
    u = box(fig, x1, yp, w1, hp)
    p = mf.read_csv("sheet_levesque_parameters.csv")
    xs = np.linspace(0.005, 0.995, 2000)
    _, curv = mf.levesque_model(xs, *(p[k][0] for k in ("alpha", "beta", "gamma", "delta")))
    u.fill_between(xs, 0.5, 0.86, where=curv > 0, color=WIN, lw=0, zorder=0)
    u.axhspan(LO, HI, color="0.86", alpha=0.85, lw=0, zorder=0.5)
    u.axhline(LO, color="0.35", lw=0.6, ls="--", zorder=1)
    u.text(0.03, 0.825, "unallowed for pristine CASS", fontsize=FS - 1, color="0.25", style="italic",
           ha="left", va="top")
    qa = mf.qe_dft_aexp()
    xq, tq = series_top("qeaexp")
    m = np.isin(np.round(xq, 3), np.round(qa["x"], 3))           # the QE tensors of Fig. 2b (80% Cr excluded)
    xq, tq = xq[m], tq[m]
    xm, tm = series_top("md5r")
    u.plot(xm, tm, "-", color=GREEN, lw=1.0, zorder=2)
    u.plot(xq, tq, "-D", color=BLUE, ms=2.6, lw=0.9, zorder=3)
    u.text(0.82, 0.585, "MD", color=GREEN, fontsize=FS - 0.5, fontweight="bold", ha="center", va="top")
    u.text(0.16, 0.540, "DFT", color=BLUE, fontsize=FS - 0.5, fontweight="bold", ha="center", va="center")
    u.set(xlim=(0, 1), ylim=(0.50, 0.86), xticks=[0, 0.5, 1], yticks=[0.5, 0.6, 0.7, 0.8])
    u.set_xticklabels(["0", "0.5", "1"])
    u.set_ylabel("top (P)SAW (GHz)", fontsize=FS - 0.5, labelpad=1)
    u.text(1.03, 0.495, r"$x_{\rm Cr}$", fontsize=FS, ha="left", va="top", clip_on=False)
    u.text(0.48, 0.512, "spinodal window", fontsize=FS - 1.5, color="#6A5A8C", ha="center", va="bottom")
    frame(u)
    cp = dict(zip(np.round(qa["x"], 2), (qa["c11"] - qa["c12"]) / 2))
    cells = ((0.05, r"$\alpha$ (5% Cr)"), (0.35, "Fe–36Cr"), (0.95, r"$\alpha'$ (95% Cr)"))
    bw = 0.54; bh = bw * VS.ASPECT; yc = 0.42
    centers = [x1 + 0.28, x1 + w1 / 2, x1 + w1 - 0.28]
    for (xc, lab), cx in zip(cells, centers):
        VS.draw(box(fig, cx - bw / 2, yc, bw, bh), xc)
        fig.text(cx / W, (yc - 0.03) / H, lab, fontsize=FS - 0.5, ha="center", va="top",
                 color=BLUE if xc == 0.35 else VERM)
        fig.text(cx / W, (yc - 0.17) / H, rf"$C'$ = {cp[round(xc, 2)]:.0f} GPa", fontsize=FS - 1, ha="center",
                 va="top", color="0.15", fontweight="bold" if xc == 0.95 else "normal")
        fig.add_artist(ConnectionPatch(xyA=(xc, 0.50), coordsA=u.transData, xyB=(cx / W, (yc + bh + 0.01) / H),
                                       coordsB=fig.transFigure, color="0.55", lw=0.5, ls=(0, (2, 1.5))))
    for xa, xb in ((centers[1] - 0.29, centers[0] + 0.29), (centers[1] + 0.29, centers[2] - 0.29)):
        fig.patches.append(FancyArrowPatch((xa / W, (yc + bh / 2) / H), (xb / W, (yc + bh / 2) / H),
                                           transform=fig.transFigure, arrowstyle="-|>", mutation_scale=7, lw=1.0,
                                           color=VERM))
    header(fig, 2.74, 5.04, yt, r"DFT: stiff $\alpha'$ enters the unallowed range", BLUE)
    fig.text(5.06 / W, 0.05 / H, "QE cells; arrows: local magnetic moments", fontsize=FS - 1.5, color="0.4",
             ha="right", va="bottom")
    fig.add_artist(plt.Line2D([2.64 / W] * 2, [0.08 / H, 1.93 / H], color="0.8", lw=0.6))

    OUT.mkdir(parents=True, exist_ok=True)
    meta = {"Title": "Graphical abstract", "Creator": "S. A. Al Dajani"}
    fig.savefig(OUT / "graphical_abstract.pdf", metadata=meta)
    fig.savefig(OUT / "graphical_abstract.png", dpi=600)
    fig.savefig(OUT / "graphical_abstract.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(OUT / "graphical_abstract_screen.png", dpi=96)
    print("Fe-36Cr stiffening", round(100 * (f[-1] / f[0] - 1), 2), "%; MD top max", tm.max().round(3),
          "; QE top", dict(zip(xq, tq.round(3))))


if __name__ == "__main__":
    main()
