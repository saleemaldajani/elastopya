#!/usr/bin/env python3
"""Graphical abstract for the Scripta Materialia submission (v7, 2026-11-01: v6 with "spinodal decomposition" in the
MD title and banner, and alpha / alpha' labels on the cells and the DFT plot as in v5).

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
MD card (v6): stiffening on decomposition at all twelve MD compositions for phi = 0.7-1 (Fig. 3b), fully decomposed
cells at six compositions under their x_Cr (Fig. 1c). DFT card (v6): five QE cells with C' (exp. a) under the top-mode plot.
Earlier versions: make_graphical_abstract_v1.py ... _v6.py (outputs in figures/graphical_abstract_v1 ... _v6).

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
from matplotlib.patches import ConnectionPatch, FancyArrowPatch, FancyBboxPatch  # noqa: E402

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


def card(fig, x0, y0, w, h, fc):
    fig.patches.append(FancyBboxPatch((x0 / W, y0 / H), w / W, h / H, transform=fig.transFigure,
                                      boxstyle="round,pad=0,rounding_size=0.025", fc=fc, ec="none", zorder=-5))


def pill(fig, x, y, tag, text, col):
    t = fig.text(x / W, y / H, tag, fontsize=FS + 0.5, fontweight="bold", color="white", ha="left", va="center",
                 bbox=dict(boxstyle="round,pad=0.28,rounding_size=0.5", fc=col, ec="none"))
    fig.canvas.draw()
    bb = t.get_window_extent(fig.canvas.get_renderer())
    fig.text(bb.x1 / fig.dpi / W + 0.07 / W, y / H, text, fontsize=FS + 0.5, fontweight="bold", color="0.12",
             ha="left", va="center")


def open_axes(ax):
    for k in ("top", "right"):
        ax.spines[k].set_visible(False)
    for k in ("left", "bottom"):
        ax.spines[k].set_linewidth(0.6); ax.spines[k].set_color("0.25")
    ax.set_facecolor("none")
    ax.tick_params(direction="out", length=2.0, width=0.6, labelsize=FS - 0.5, pad=1.2, colors="0.2")


def main():
    plt.rcParams.update({"font.family": "Arial", "mathtext.fontset": "custom", "mathtext.rm": "Arial",
                         "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "font.size": FS,
                         "pdf.fonttype": 42, "axes.linewidth": 0.6})
    fig = plt.figure(figsize=(W, H))
    MDT, DFTT = "#E9F5EF", "#E7F0F8"                # card tints
    yp, hp = 1.04, 0.55                           # plots: bottom, height (inches)
    card(fig, 0.06, 0.31, 2.50, 1.62, MDT)
    card(fig, 2.62, 0.31, 2.44, 1.62, DFTT)

    # ---------------- MD: stiffening at every composition, decomposed cells under their compositions
    pill(fig, 0.16, 1.78, "MD", "Spinodal decomposition stiffens the SAW", GREEN)
    x0, w = 0.50, 1.86
    a = box(fig, x0, yp, w, hp)
    p = mf.read_csv("sheet_levesque_parameters.csv")
    xs = np.linspace(0.005, 0.995, 2000)
    _, curv = mf.levesque_model(xs, *(p[k][0] for k in ("alpha", "beta", "gamma", "delta")))
    a.fill_between(xs, -0.5, 8.2, where=curv > 0, color="#D9D0EA", alpha=0.55, lw=0, zorder=0)
    a.text(0.235, 8.1, "spinodal window", fontsize=FS - 1.5, color="#6A5A8C", ha="left", va="top")
    dec = mf.md_decomposition()
    xd = np.array(sorted(dec))
    rel = lambda ph: np.array([100 * (dec[x]["f"][np.argmin(abs(dec[x]["phi"] - ph))] /
                                      dec[x]["f"][np.argmin(dec[x]["phi"])] - 1) for x in xd])
    for ph, al in ((0.7, 0.30), (0.8, 0.50), (0.9, 0.72)):
        a.plot(xd, rel(ph), "-", color=GREEN, alpha=al, lw=0.9, zorder=2)
        a.text(xd[6] + 0.015, rel(ph)[6] + 0.25, rf"{ph:g}", color=GREEN, alpha=min(1, al + 0.25), fontsize=FS - 2,
               ha="left", va="bottom")
    r1 = rel(1.0)
    a.fill_between(xd, 0, r1, color=GREEN, alpha=0.15, lw=0, zorder=1)
    a.plot(xd, r1, "-o", color=GREEN, mfc="white", mew=1.0, ms=3.0, lw=1.4, zorder=3)
    k = int(np.argmax(r1))
    a.text(xd[k] + 0.04, r1[k] - 0.05, f"+{r1[k]:.1f}%", color=VERM, fontsize=FS + 1, fontweight="bold", ha="left",
           va="bottom")
    a.text(0.75, 4.4, r"$\phi$ = 1", color=GREEN, fontsize=FS - 1, ha="left", va="bottom", fontweight="bold")
    a.set(xlim=(0, 1), ylim=(-0.3, 8.2), xticks=[0, 0.5, 1], yticks=[0, 2, 4, 6, 8])
    a.set_xticklabels(["0", "0.5", "1"])
    open_axes(a)
    a.set_ylabel(r"$\Delta f_{\rm SAW}$ (%)", fontsize=FS - 0.5, labelpad=1, color="0.2")
    a.text(1.03, -0.9, r"$x_{\rm Cr}$", fontsize=FS, ha="left", va="top", clip_on=False, color="0.2")
    cw, ch = 0.25, 0.29
    mx = (0.05, 0.2, 0.36, 0.5, 0.7, 0.95)
    ser = series()
    for n, xc in enumerate(mx):
        cx = x0 + 0.08 + n * (w - 0.16) / (len(mx) - 1)
        _cube_ax(box(fig, cx - cw / 2, 0.50, cw, ch), md_grid3d(ser[xc][-1]))
        fig.text(cx / W, 0.47 / H, f"{round(100 * xc)}%", fontsize=FS - 1, ha="center", va="top",
                 color=BLUE if xc == 0.36 else "0.2", fontweight="bold" if xc == 0.36 else "normal")
        fig.add_artist(ConnectionPatch(xyA=(xc, -0.3), coordsA=a.transData, xyB=(cx / W, (0.50 + ch + 0.01) / H),
                                       coordsB=fig.transFigure, color="0.6", lw=0.5, ls=(0, (2, 1.5))))
    fig.text(0.14 / W, 0.645 / H, r"$\alpha$ + $\alpha'$" + "\n" + r"($\phi$ = 1)", fontsize=FS - 1.5, ha="left",
             va="center", color=VERM, linespacing=1.0)

    # ---------------- DFT: top surface mode at every QE composition, five cells with their C'
    pill(fig, 2.72, 1.78, "DFT", r"Stiff $\alpha'$ enters the unallowed range", BLUE)
    x1, w1 = 3.10, 1.74
    u = box(fig, x1, yp, w1, hp)
    u.fill_between(xs, 0.5, 0.86, where=curv > 0, color="#D9D0EA", alpha=0.55, lw=0, zorder=0)
    u.axhspan(LO, 0.86, color="#FBE0D2", lw=0, zorder=0.5)
    u.axhline(LO, color=VERM, lw=0.6, ls=(0, (2, 2)), zorder=1)
    u.text(0.03, 0.845, "unallowed for pristine CASS", fontsize=FS - 1, color="#9A3A10", style="italic",
           ha="left", va="top")
    qa = mf.qe_dft_aexp()
    xq, tq = series_top("qeaexp")
    m = np.isin(np.round(xq, 3), np.round(qa["x"], 3))           # the QE tensors of Fig. 2b (80% Cr excluded)
    xq, tq = xq[m], tq[m]
    xm, tm = series_top("md5r")
    u.plot(xm, tm, "-", color=GREEN, lw=1.3, zorder=2, solid_capstyle="round")
    u.plot(xq, tq, "-", color=BLUE, lw=1.1, zorder=3)
    u.plot(xq[:-1], tq[:-1], "D", color=BLUE, mec="white", mew=0.4, ms=3.0, zorder=4)
    u.plot(xq[-1:], tq[-1:], "o", ms=9, mfc="none", mec=VERM, mew=1.0, zorder=4)
    u.plot(xq[-1:], tq[-1:], "D", color=VERM, mec="white", mew=0.4, ms=3.6, zorder=5)
    u.text(xq[-1] - 0.06, tq[-1] + 0.005, r"$\alpha'$", color=VERM, fontsize=FS + 1, fontweight="bold",
           ha="right", va="center")
    u.text(xq[0] + 0.02, 0.528, r"$\alpha$", color=VERM, fontsize=FS + 1, fontweight="bold", ha="center",
           va="center")
    u.text(0.84, 0.585, "MD", color=GREEN, fontsize=FS - 0.5, fontweight="bold", ha="center", va="top")
    u.text(0.30, 0.527, "DFT", color=BLUE, fontsize=FS - 0.5, fontweight="bold", ha="center", va="center")
    u.set(xlim=(0, 1), ylim=(0.50, 0.86), xticks=[0, 0.5, 1], yticks=[0.5, 0.6, 0.7, 0.8])
    u.set_xticklabels(["0", "0.5", "1"])
    open_axes(u)
    u.set_ylabel("top (P)SAW (GHz)", fontsize=FS - 0.5, labelpad=1, color="0.2")
    u.text(1.03, 0.495, r"$x_{\rm Cr}$", fontsize=FS, ha="left", va="top", clip_on=False, color="0.2")
    cp = dict(zip(np.round(qa["x"], 2), (qa["c11"] - qa["c12"]) / 2))
    qx = (0.05, 0.2, 0.35, 0.5, 0.95)
    bw = 0.37; bh = bw * VS.ASPECT; yc = 0.62
    for n, xc in enumerate(qx):
        cx = x1 + 0.08 + n * (w1 - 0.10) / (len(qx) - 1)
        VS.draw(box(fig, cx - bw / 2, yc, bw, bh), xc)
        hi = xc == 0.95
        col = VERM if hi else (BLUE if xc == 0.35 else "0.2")
        lab = {0.05: r"$\alpha$  5%", 0.95: r"$\alpha'$  95%"}.get(xc, f"{round(100 * xc)}%")
        fig.text(cx / W, (yc - 0.02) / H, lab, fontsize=FS - 1, ha="center", va="top",
                 color=VERM if xc in (0.05, 0.95) else col, fontweight="bold" if xc in (0.05, 0.35, 0.95) else "normal")
        fig.text(cx / W, (yc - 0.15) / H, f"{cp[round(xc, 2)]:.0f}", fontsize=FS - 0.5, ha="center", va="top",
                 color=VERM if hi else "0.2", fontweight="bold" if hi else "normal")
        fig.add_artist(ConnectionPatch(xyA=(xc, 0.50), coordsA=u.transData, xyB=(cx / W, (yc + bh + 0.005) / H),
                                       coordsB=fig.transFigure, color="0.6", lw=0.5, ls=(0, (2, 1.5))))
    fig.text(2.68 / W, (yc - 0.15) / H, r"$C'$ (GPa)", fontsize=FS - 1.5, ha="left", va="top", color="0.2")
    fig.text(4.98 / W, 1.635 / H, "QE cells; arrows: local moments", fontsize=FS - 2, color="0.45",
             ha="right", va="center")

    # ---------------- banner
    fig.text(0.5, 0.16 / H, r"Spinodal decomposition $\rightarrow$ stiffer $\delta$-ferrite ($C'$ $\uparrow$) $\rightarrow$ "
             r"higher SAW frequency $\rightarrow$ non-destructive sign of embrittlement",
             fontsize=FS - 0.5, ha="center", va="center", color="0.12")

    OUT.mkdir(parents=True, exist_ok=True)
    meta = {"Title": "Graphical abstract", "Creator": "S. A. Al Dajani"}
    fig.savefig(OUT / "graphical_abstract.pdf", metadata=meta)
    fig.savefig(OUT / "graphical_abstract.png", dpi=600)
    fig.savefig(OUT / "graphical_abstract.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(OUT / "graphical_abstract_screen.png", dpi=96)
    print("MD stiffening at phi = 1:", dict(zip(xd, r1.round(2))), "; MD top max", tm.max().round(3),
          "; QE top", dict(zip(xq, tq.round(3))))


if __name__ == "__main__":
    main()
