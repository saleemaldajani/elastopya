#!/usr/bin/env python3
"""Graphical abstract for the Scripta Materialia submission (v11, 2026-11-01).

Modelled on the graphical abstract of Al Dajani et al., Acta Mater. 2023 (three panels, one strong image each, a bold
claim underneath), as its theoretical counterpart (user request; no spinodal shading, which the user found confusing):
  1  MD: fully decomposed Fe-36Cr (25^3 cell, phi = 1) separated into Fe-rich alpha and Cr-rich alpha' (cf. the TEM map)
  2  (P)SAW mode frequencies over 10^4 grain orientations at 6.4 um (Fig. 2b data; cf. the measured SAW spectrum): homogeneous
     Fe-36Cr (QE, experimental a) and decomposed MD Fe-36Cr stay below the range unallowed for pristine CASS (0.609 GHz);
     the QE alpha' tensor (95% Cr) reaches 0.69 GHz; inset zooms on the range, as in the Acta inset
  3  MD at all twelve compositions: on full decomposition f_SAW rises and Pugh's K/G falls (VRH moduli of the same tensors;
     at Fe-36Cr G +11%, K -1.2%, K/G -11%, Sec. 3.2)
Elsevier: at least 1328 x 531 px, readable at 13 x 5 cm; drawn at 13 x 5 cm, exported PDF/PNG/TIFF (600 dpi).
Versions with their scripts: manuscript/figures/graphical_abstract_versions/ (VERSIONS.md).
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyBboxPatch, Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_figures as mf  # noqa: E402
from visualize_decomposition import CR, FE, _cube_ax, md_grid3d, series  # noqa: E402

OUT = mf.ROOT / "manuscript" / "figures" / "graphical_abstract"
W, H = 13 / 2.54, 5 / 2.54
BLUE, VERM, GREEN, GREY = "#0072B2", "#D55E00", "#009E73", "0.45"
FS = 7.0
LO, HI = 0.609, 0.833                            # unallowed range for pristine CASS at 6.4 um (Fig. 2b)


def box(fig, x0, y0, w, h):
    return fig.add_axes([x0 / W, y0 / H, w / W, h / H])


def frame(ax):
    ax.tick_params(direction="in", length=2.2, width=0.6, labelsize=FS - 1, pad=1.5, colors="0.15")
    for s in ax.spines.values():
        s.set_linewidth(0.6)


def caption(fig, x0, w, text):
    """Bold claim in a light grey strip under each panel, as in the Acta graphical abstract."""
    fig.patches.append(Rectangle((x0 / W, 0.03 / H), w / W, 0.36 / H, transform=fig.transFigure, fc="0.92",
                                 ec="none", zorder=-1))
    fig.text((x0 + w / 2) / W, 0.21 / H, text, fontsize=FS - 0.5, fontweight="bold", ha="center", va="center",
             linespacing=1.15, color="0.1")


def modes(tag):
    d = np.loadtxt(mf.DATA / "saw" / f"{tag}.csv", delimiter=",", skiprows=1)
    keep = d[:, 3] / d[:, 2] >= 2e-3                     # pseudo-SAW kept as in Fig. 2b
    return np.r_[d[:, 0], d[keep, 1]] / 6.4e-6 / 1e9


def kg(x, ph):
    d = mf.read_csv("md_decomposition.csv")
    m = (np.abs(d["x_cr"] - x) < 1e-6) & (np.abs(d["phi"] - ph) < 1e-6)
    c11, c12, c44 = (float(np.mean(d[k][m])) for k in ("c11", "c12", "c44"))
    K, cp = (c11 + 2 * c12) / 3, (c11 - c12) / 2
    G = ((2 * cp + 3 * c44) / 5 + 5 * (2 * cp) * c44 / (4 * c44 + 3 * (2 * cp))) / 2
    return K / G


def main():
    plt.rcParams.update({"font.family": "Arial", "mathtext.fontset": "custom", "mathtext.rm": "Arial",
                         "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "font.size": FS,
                         "pdf.fonttype": 42, "axes.linewidth": 0.6})
    fig = plt.figure(figsize=(W, H))
    yb, ht = 0.52, 1.36                            # panel area above the caption strips

    # ---------------- 1: MD decomposition of Fe-36Cr
    x0, w0 = 0.05, 1.16
    c = box(fig, x0 + 0.16, yb + 0.06, 0.86, 1.04)
    _cube_ax(c, md_grid3d(series()[0.36][-1]))
    for k, (name, col, lab) in enumerate((("Fe", FE, r"Fe-rich $\alpha$"), ("Cr", CR, r"Cr-rich $\alpha'$"))):
        y = yb + ht - 0.02
        xk = x0 + 0.04 + k * 0.58
        fig.patches.append(Rectangle((xk / W, (y - 0.04) / H), 0.08 / W, 0.08 / H,
                                     transform=fig.transFigure, fc=col, ec="none"))
        fig.text((xk + 0.10) / W, y / H, lab, fontsize=FS - 1.5, va="center", color="0.15")
    fig.text((x0 + w0 - 0.02) / W, (yb + 0.02) / H, r"MD, 25$^3$ cell", fontsize=FS - 1.5, ha="right", va="bottom",
             color=GREY)
    caption(fig, x0, w0, "Simulated spinodal\ndecomposition of Fe–36Cr")

    # ---------------- 2: surface-mode spectrum and the unallowed range
    x1, w1 = 1.30, 2.10
    s = box(fig, x1 + 0.12, yb + 0.20, w1 - 0.20, ht - 0.26)
    s.axvspan(LO, HI, color="#FBE0D2", lw=0, zorder=0)
    s.axvline(LO, color=VERM, lw=0.6, ls=(0, (2, 2)), zorder=1)
    bins = np.arange(0.35, 0.80, 0.005)
    series_ = ((modes("qeaexp_x0.35"), BLUE, "Fe–36Cr, homogeneous (DFT)", "-"),
               (modes("mddec_x0.36"), GREEN, "Fe–36Cr, decomposed (MD)", "-"),
               (modes("qeaexp_x0.95"), VERM, r"Cr-rich $\alpha'$, 95% Cr (DFT)", "-"))
    hists = []
    for m, col, lab, ls in series_:
        h, e = np.histogram(m, bins=bins, density=True)
        xc = (e[:-1] + e[1:]) / 2
        hists.append((xc, h, col))
        s.fill_between(xc, 0, h, color=col, alpha=0.18, lw=0, step="mid")
        s.plot(xc, h, color=col, lw=1.0, drawstyle="steps-mid", label=lab)
    s.text(LO + 0.012, s.get_ylim()[1] * 0.93 if False else 0, "", fontsize=1)
    s.set(xlim=(0.35, 0.80), xticks=[0.4, 0.5, 0.6, 0.7, 0.8], yticks=[])
    s.set_xlabel(r"(P)SAW frequency (GHz), $\lambda$ = 6.4 $\mu$m", fontsize=FS - 1, labelpad=1)
    s.set_ylabel("grains", fontsize=FS - 1, labelpad=2)
    top = max(h.max() for _, h, _ in hists)
    s.set_ylim(0, top * 1.42)
    s.text(0.795, top * 1.36, "unallowed for\npristine CASS", fontsize=FS - 1.5, color="#9A3A10", style="italic",
           ha="right", va="top", linespacing=1.0)
    s.legend(loc="upper left", fontsize=FS - 1.5, frameon=False, handlelength=1.0, borderaxespad=0.2,
             labelspacing=0.15, handletextpad=0.4)
    frame(s)
    # inset: zoom on the unallowed range, as in the Acta inset
    ins = s.inset_axes([0.60, 0.20, 0.38, 0.40])
    ins.axvspan(LO, HI, color="#FBE0D2", lw=0, zorder=0)
    ins.axvline(LO, color=VERM, lw=0.5, ls=(0, (2, 2)))
    for xc, h, col in hists:
        ins.plot(xc, h, color=col, lw=0.9, drawstyle="steps-mid")
        ins.fill_between(xc, 0, h, color=col, alpha=0.18, lw=0, step="mid")
    ymax = max(h[(xc > LO)].max() for xc, h, _ in hists) * 1.6
    ins.set(xlim=(0.58, 0.72), ylim=(0, ymax), xticks=[0.6, 0.7], yticks=[])
    ins.tick_params(labelsize=FS - 2.5, length=1.5, width=0.5, pad=1, direction="in")
    for sp in ins.spines.values():
        sp.set_linewidth(0.6); sp.set_color(VERM)
    s.indicate_inset_zoom(ins, edgecolor=VERM, lw=0.6, alpha=0.8)
    caption(fig, x1, w1, "DFT: stiff Cr-rich α′ places surface waves in the\nrange unallowed for pristine ferrite")

    # ---------------- 3: f_SAW up, K/G down at every composition
    x2, w2 = 3.49, 1.58
    k = box(fig, x2 + 0.26, yb + 0.20, w2 - 0.32, ht - 0.26)
    dec = mf.md_decomposition()
    xd = np.array(sorted(dec))
    df = np.array([100 * (dec[x]["f"][np.argmax(dec[x]["phi"])] / dec[x]["f"][np.argmin(dec[x]["phi"])] - 1)
                   for x in xd])
    dkg = np.array([100 * (kg(x, 1.0) / kg(x, 0.5) - 1) for x in xd])
    k.axhline(0, color="0.6", lw=0.6)
    k.fill_between(xd, 0, df, color=GREEN, alpha=0.18, lw=0)
    k.fill_between(xd, 0, dkg, color=VERM, alpha=0.15, lw=0)
    k.plot(xd, df, "-o", color=GREEN, ms=2.4, lw=1.2, mfc="white", mew=0.8)
    k.plot(xd, dkg, "-o", color=VERM, ms=2.4, lw=1.2, mfc="white", mew=0.8)
    k.text(0.5, df.max() + 1.0, r"$f_{\rm SAW}$ $\uparrow$  (signal)", color=GREEN, fontsize=FS - 0.5,
           fontweight="bold", ha="center", va="bottom")
    k.text(0.5, dkg.min() - 1.0, r"$K/G$ $\downarrow$  (embrittlement)", color=VERM, fontsize=FS - 0.5, fontweight="bold",
           ha="center", va="top")
    i36 = int(np.argmin(abs(xd - mf.X_FE36CR)))
    k.plot([mf.X_FE36CR] * 2, [dkg[i36], df[i36]], color=BLUE, lw=0.6, ls=(0, (2, 2)))
    k.text(mf.X_FE36CR - 0.02, -3.0, "Fe–36Cr", color=BLUE, fontsize=FS - 1.5, rotation=90, ha="right",
           va="center")
    k.set(xlim=(0, 1), ylim=(-16, 11), xticks=[0, 0.5, 1], yticks=[-10, 0, 10])
    k.set_xticklabels(["0", "0.5", "1"])
    k.set_xlabel(r"$x_{\rm Cr}$", fontsize=FS - 0.5, labelpad=1)
    k.set_ylabel("change on decomposition (%)", fontsize=FS - 1.5, labelpad=1)
    frame(k)
    caption(fig, x2, w2, "Shear stiffening raises $f_{\\rm SAW}$\nand lowers $K/G$: a signature\nof embrittlement")

    OUT.mkdir(parents=True, exist_ok=True)
    meta = {"Title": "Graphical abstract", "Creator": "S. A. Al Dajani"}
    fig.savefig(OUT / "graphical_abstract.pdf", metadata=meta)
    fig.savefig(OUT / "graphical_abstract.png", dpi=600)
    fig.savefig(OUT / "graphical_abstract.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    fig.savefig(OUT / "graphical_abstract_screen.png", dpi=96)
    print("df", dict(zip(xd, df.round(1))), "dK/G", dict(zip(xd, dkg.round(1))))


if __name__ == "__main__":
    main()
