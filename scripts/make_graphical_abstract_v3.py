#!/usr/bin/env python3
"""Graphical abstract for the Scripta Materialia submission (2026-11-01).

Elsevier: one panel, at least 1328 x 531 px (w x h), readable at 13 x 5 cm. It is drawn at exactly 13 x 5 cm, so the
lettering sizes are the printed sizes, and exported as PDF (vector), PNG and TIFF (600 dpi, 3071 x 1181 px).

Left to right it follows the main-text figures, with the same data (v3, 2026-11-01: MD and DFT side by side, and the
unallowed range of Fig. 2b):
  1  unstable window: QE demixing free energy (Fig. 1a), shaded where the PWSCF-PAW model curvature is positive (Fig. 1b)
  2  MD: Fe-36Cr decomposing, random alloy -> alpha + alpha' (Figs. 1c, 3c), and the Rayleigh SAW frequency of the random
     and fully decomposed MD alloys at twelve compositions (Fig. 3)
  3  DFT (QE, experimental lattice parameter): C' rises towards Cr (d-band filling, Sec. 3), and the highest (P)SAW
     frequency over 10^4 grain orientations enters the range unallowed for pristine CASS only for the Cr-rich tensors,
     whereas no MD tensor, random or decomposed, reaches it (Fig. 2b)
  banner: the result (Sec. 3). Earlier versions: make_graphical_abstract_v1.py, _v2.py
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


def top_mode_ghz(path):
    """Highest (P)SAW frequency over the 10^4 grain orientations, with the intensity cut of Fig. 2b."""
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    keep = d[:, 3] / d[:, 2] >= 2e-3
    return float(np.r_[d[:, 0], d[keep, 1]].max() / 6.4e-6 / 1e9)


def series_top(prefix):
    import glob
    import re
    rows = sorted((float(re.search(r"_x([\d.]+)\.csv", p).group(1)), top_mode_ghz(p))
                  for p in glob.glob(str(mf.DATA / "saw" / f"{prefix}_x*.csv")))
    return np.array([r[0] for r in rows]), np.array([r[1] for r in rows])


def header(fig, x0, x1, num, text, col):
    title(fig, x0, YT, num, text)
    fig.add_artist(plt.Line2D([x0 / W, x1 / W], [(YT - 0.12) / H] * 2, color=col, lw=1.2))


YT = 1.82


def main():
    plt.rcParams.update({"font.family": "Arial", "mathtext.fontset": "custom", "mathtext.rm": "Arial",
                         "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "font.size": FS,
                         "pdf.fonttype": 42, "axes.linewidth": 0.6})
    fig = plt.figure(figsize=(W, H))
    yb, hb = 0.58, 0.96                           # axes bottom, height
    pw = 0.70                                     # plot width
    xlab = lambda ax: ax.set_xlabel(r"$x_{\rm Cr}$", labelpad=0.3)

    # 1 unstable window (Fig. 1): QE demixing free energy
    a = box(fig, 0.40, yb, pw, hb)
    shade_window(a, -0.08, 0.02)
    q = mf.qe_mixing()
    a.axhline(0, color="0.6", lw=0.5, zorder=1)
    a.plot(q["x"], q["demix"], "-o", color=BLUE, mfc="white", mew=0.8, ms=2.5, lw=0.8, zorder=3)
    a.plot([mf.X_FE36CR], [-0.074], "*", color=BLUE, ms=6, zorder=4)
    a.text(0.5, 0.013, "unstable", fontsize=FS - 1, ha="center", va="top", color="0.25")
    a.text(0.05, -0.05, "DFT", fontsize=FS - 1, color=BLUE, ha="left", va="center")
    a.set(xlim=(0, 1), ylim=(-0.08, 0.02), xticks=[0, 0.5, 1], yticks=[-0.06, 0])
    a.set_xticklabels(["0", "0.5", "1"]); a.set_yticklabels(["−0.06", "0"])
    xlab(a)
    a.set_ylabel(r"$\Delta G_{\rm demix}$ (eV/atom)", labelpad=0.3, fontsize=FS - 1)
    frame(a)
    header(fig, 0.12, 1.10, "1", "Unstable window", "0.55")

    # 2 MD: decomposition of Fe-36Cr (cubes) and the SAW frequency it raises
    f36 = series()[0.36]
    xc = 1.21
    for i, (f, lab) in enumerate(((f36[0], "random"), (f36[-1], r"$\alpha$ + $\alpha'$"))):
        y0 = yb + 0.50 - i * 0.52
        c = box(fig, xc, y0, 0.34, 0.40)
        _cube_ax(c, md_grid3d(f))
        if i == 0:
            fig.text((xc + 0.17) / W, (y0 + 0.42) / H, lab, fontsize=FS - 1.5, ha="center", va="bottom")
        else:
            fig.text((xc + 0.17) / W, (y0 - 0.01) / H, lab, fontsize=FS - 1.5, ha="center", va="top")
    fig.patches.append(FancyArrowPatch(((xc + 0.17) / W, (yb + 0.49) / H), ((xc + 0.17) / W, (yb + 0.40) / H),
                                       transform=fig.transFigure, arrowstyle="-|>", mutation_scale=7, lw=1.0,
                                       color=GREY))
    s = box(fig, 1.93, yb, pw, hb)
    shade_window(s, 0.43, 0.545)
    dec = mf.md_decomposition()
    xs3 = np.array(sorted(dec))
    f0 = np.array([dec[x]["f"][np.argmin(dec[x]["phi"])] for x in xs3])
    f1 = np.array([dec[x]["f"][np.argmax(dec[x]["phi"])] for x in xs3])
    s.plot(xs3, f0, "-s", color=GREEN, ms=2.2, lw=0.8, zorder=3, label="random")
    s.plot(xs3, f1, "--s", color=GREEN, mfc="white", mew=0.7, ms=2.4, lw=0.8, zorder=3, label="decomposed")
    i = int(np.argmax(f1 / f0))
    s.annotate("", xy=(xs3[i], f1[i] - 0.002), xytext=(xs3[i], f0[i] + 0.002),
               arrowprops=dict(arrowstyle="-|>", lw=0.9, color=VERM, mutation_scale=6), zorder=4)
    s.text(xs3[i] + 0.03, 0.4415, f"+{100 * (f1[i] / f0[i] - 1):.1f}%", color=VERM, fontsize=FS,
           fontweight="bold", va="center", ha="left")
    s.plot([mf.X_FE36CR], [0.4355], "*", color=BLUE, ms=6, zorder=4)
    s.legend(loc="upper left", fontsize=FS - 1.5, frameon=False, handlelength=1.4, borderaxespad=0.15,
             labelspacing=0.05, handletextpad=0.3)
    s.set(xlim=(0, 1), ylim=(0.43, 0.545), xticks=[0, 0.5, 1], yticks=[0.45, 0.50])
    s.set_xticklabels(["0", "0.5", "1"])
    xlab(s)
    s.set_ylabel(r"$f_{\rm SAW}$ (GHz)", labelpad=0.3, fontsize=FS - 1)
    frame(s)
    header(fig, 1.20, 2.63, "2", "MD: SAW stiffening", GREEN)

    # 3 DFT (QE, experimental a): stiff Cr-rich alpha', whose surface modes reach the unallowed range
    qa = mf.qe_dft_aexp()
    e = box(fig, 2.99, yb, pw, hb)
    shade_window(e, 25, 175)
    e.plot(qa["x"], (qa["c11"] - qa["c12"]) / 2, "-D", color=BLUE, ms=2.3, lw=0.8, zorder=3)
    ya = 36
    for xe in (0.15, 0.85):
        e.annotate("", xy=(xe, ya), xytext=(mf.X_FE36CR, ya),
                   arrowprops=dict(arrowstyle="-|>", lw=0.9, color=VERM, mutation_scale=6, shrinkA=4))
    e.plot([mf.X_FE36CR], [ya], "*", color=BLUE, ms=6, zorder=4)
    e.text(0.12, ya, r"$\alpha$", color=VERM, fontsize=FS, ha="right", va="center")
    e.text(0.87, ya, r"$\alpha'$", color=VERM, fontsize=FS, ha="left", va="center")
    e.text(0.40, 112, r"$d$-band" + "\nfilling", fontsize=FS - 1.5, color="0.4", style="italic", ha="center",
           va="center", linespacing=0.95)
    e.annotate("", xy=(0.80, 128), xytext=(0.58, 108),
               arrowprops=dict(arrowstyle="-|>", lw=0.6, color="0.5", mutation_scale=5))
    e.set(xlim=(0, 1), ylim=(25, 175), xticks=[0, 0.5, 1], yticks=[50, 100, 150])
    e.set_xticklabels(["0", "0.5", "1"])
    xlab(e)
    e.set_ylabel(r"$C'$ (GPa)", labelpad=0.3, fontsize=FS - 1)
    frame(e)

    u = box(fig, 4.10, yb, pw + 0.06, hb)
    lo, hi = 0.609, 0.833
    u.axhspan(lo, hi, color="0.88", lw=0, zorder=0)
    u.axhline(lo, color="0.4", lw=0.6, ls="--", zorder=1)
    u.text(0.04, 0.80, "unallowed\n(pristine)", fontsize=FS - 1.5, color="0.3", style="italic", ha="left", va="top",
           linespacing=0.95)
    xm, tm = series_top("md5r")
    xd, td = series_top("mddec")
    xq, tq = series_top("qeaexp")
    m = np.isin(np.round(xq, 3), np.round(qa["x"], 3))     # the QE tensors of Fig. 2b (80% Cr AF branch excluded)
    xq, tq = xq[m], tq[m]
    u.plot(xm, tm, "-s", color=GREEN, ms=2.2, lw=0.8, zorder=3)
    u.plot(xd, td, "--s", color=GREEN, mfc="white", mew=0.7, ms=2.4, lw=0.8, zorder=3)
    u.plot(xq, tq, "-D", color=BLUE, ms=2.3, lw=0.8, zorder=4)
    u.text(0.72, 0.570, "MD", color=GREEN, fontsize=FS - 1, fontweight="bold", ha="center", va="center")
    u.text(0.30, 0.535, "DFT", color=BLUE, fontsize=FS - 1, fontweight="bold", ha="center", va="center")
    u.text(0.90, tq[-1], r"$\alpha'$", color=VERM, fontsize=FS, ha="right", va="center")
    u.set(xlim=(0, 1), ylim=(0.50, 0.86), xticks=[0, 0.5, 1], yticks=[0.5, 0.6, 0.7, 0.8])
    u.set_xticklabels(["0", "0.5", "1"])
    xlab(u)
    u.set_ylabel("top (P)SAW (GHz)", labelpad=0.3, fontsize=FS - 1)
    frame(u)
    header(fig, 2.80, 4.92, "3", r"DFT: stiff $\alpha'$ reaches the unallowed range", BLUE)

    for x0, x1 in ((1.12, 1.19), (1.57, 1.64), (2.65, 2.72), (3.71, 3.78)):
        arrow(fig, x0, x1, yb + hb / 2)

    # banner: the result
    fig.patches.append(FancyBboxPatch((0.06 / W, 0.03 / H), (W - 0.12) / W, 0.21 / H, transform=fig.transFigure,
                                      boxstyle="round,pad=0.0,rounding_size=0.01", fc="#F4F4F4", ec="0.65", lw=0.5))
    fig.text(0.5, 0.135 / H, r"Spinodal decomposition:  $C'$ $\uparrow$  $\rightarrow$  $f_{\rm SAW}$ $\uparrow$,  "
             r"$K/G$ $\downarrow$   =   a non-destructive indicator of embrittlement", fontsize=FS, ha="center",
             va="center", color="0.1")

    OUT.mkdir(parents=True, exist_ok=True)
    meta = {"Title": "Graphical abstract", "Creator": "S. A. Al Dajani"}
    fig.savefig(OUT / "graphical_abstract.pdf", metadata=meta)
    fig.savefig(OUT / "graphical_abstract.png", dpi=600)
    fig.savefig(OUT / "graphical_abstract.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(OUT / "graphical_abstract_screen.png", dpi=96)    # what readers see on the journal page (13 x 5 cm)
    print("max stiffening:", xs3[i], round(100 * (f1[i] / f0[i] - 1), 2), "%")
    print("top (P)SAW: MD random max", tm.max().round(3), "MD dec max", td.max().round(3), "QE aexp", dict(zip(xq, tq.round(3))))


if __name__ == "__main__":
    main()
