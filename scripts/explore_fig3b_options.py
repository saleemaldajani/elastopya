#!/usr/bin/env python3
"""Exploration (2026-10-16): ways to show all twelve MD compositions in Fig. 3b without a tangle.
Four candidate panels at the size of Fig. 3b, side by side for comparison (not a paper figure):
  A  absolute f_SAW vs phi, every composition, coloured by x_Cr (colour bar), Acta 50^3 in black
  B  relative change (f - f(phi=0.5))/f(phi=0.5) vs phi, same colouring: the curves fan out by stiffening
  C  map of the relative change over (x_Cr, phi), spinodal window (CALPHAD, 748 K) and Fe-36Cr marked
  D  stiffening at phi = 1 vs x_Cr (five-alloy mean +- s.d.), with the phi = 0.7 and 0.9 steps
  explore_fig3b_options.py -> manuscript/figures/explore/fig3b_options.png
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_figures_kdense as mk  # noqa: E402,F401  (Scripta style on make_figures)
from make_figures_kdense import mf, plt  # noqa: E402

OUT = mf.ROOT / "manuscript" / "figures" / "explore"


def main():
    dec = mf.md_decomposition()
    acta = mf.md()
    xs = sorted(dec)
    cmap = plt.get_cmap("viridis")
    norm = plt.Normalize(0, 1)
    rel = lambda s: 100 * (s["f"] / s["f"][np.argmin(s["phi"])] - 1)
    fig, axs = plt.subplots(1, 4, figsize=(14.4, 3.0))
    a, b, c, d = axs
    for x in xs:
        s = dec[x]; o = np.argsort(s["phi"])
        a.plot(s["phi"][o], s["f"][o], "-s", color=cmap(norm(x)), ms=3, lw=0.9)
        b.plot(s["phi"][o], rel(s)[o], "-s", color=cmap(norm(x)), ms=3, lw=0.9)
    a.plot(acta["phi"], acta["f"], "-o", color="black", ms=4, lw=0.9, label=r"Fe$_{50}$Cr$_{50}$, Acta 50$^3$")
    ra = 100 * (acta["f"] / acta["f"][np.argmin(acta["phi"])] - 1)
    b.plot(acta["phi"], ra, "-o", color="black", ms=4, lw=0.9, label=r"Fe$_{50}$Cr$_{50}$, Acta 50$^3$")
    for ax, yl, t in ((a, r"$f_{\rm SAW}$ (GHz), 6.4 $\mu$m", "(A) absolute"), (b, r"$\Delta f_{\rm SAW}/f_{\rm SAW}(\phi=0.5)$ (%)", "(B) relative")):
        ax.set(xlabel=r"$\phi$", ylabel=yl); ax.legend(frameon=False, fontsize=6, loc="upper left"); mf.style(ax, t)
        for x in xs:   # end labels for the extremes and Fe-36Cr only
            if round(x, 2) in (0.05, 0.36, 0.5, 0.95):
                s = dec[x]; o = np.argsort(s["phi"]); y = (s["f"] if ax is a else rel(s))[o][-1]
                ax.text(1.01, y, f"{round(100 * x)}%", color=cmap(norm(x)), fontsize=6, va="center")
        ax.set_xlim(0.48, 1.08)
    fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=b, pad=0.01, fraction=0.05).set_label(r"$x_{\rm Cr}$")
    # C: map
    phis = np.unique(np.round(np.concatenate([dec[x]["phi"] for x in xs]), 2))
    Z = np.array([[float(np.interp(p, np.sort(dec[x]["phi"]), rel(dec[x])[np.argsort(dec[x]["phi"])])) for x in xs] for p in phis])
    im = c.pcolormesh(np.array(xs), phis, Z, shading="nearest", cmap="magma_r")
    fig.colorbar(im, ax=c, pad=0.01, fraction=0.05).set_label(r"$\Delta f/f$ (%)")
    c.axvline(0.30, color="white", lw=0.8, ls="--"); c.axvline(0.71, color="white", lw=0.8, ls="--")
    c.plot([0.36], [1.0], "*", color="#56B4E9", ms=9, mec="white")
    c.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"$\phi$"); mf.style(c, "(C) map")
    # D: stiffening vs composition
    for p, mkr, lab in ((0.7, "^", r"$\phi$ = 0.7"), (0.9, "v", r"$\phi$ = 0.9"), (1.0, "s", r"$\phi$ = 1")):
        y, e = [], []
        for x in xs:
            s = dec[x]; i = int(np.argmin(abs(s["phi"] - p))); i0 = int(np.argmin(s["phi"]))
            y.append(100 * (s["f"][i] / s["f"][i0] - 1))
            e.append(100 * s["f_sd"][i] / s["f"][i0] if s.get("f_sd") is not None else 0)
        d.errorbar(xs, y, yerr=e, fmt=mkr + "-", ms=4, lw=0.8, capsize=1.5, label=lab,
                   color={0.7: "0.65", 0.9: "0.4", 1.0: mf.GREEN}[p])
    d.axvspan(0.30, 0.71, color="#FFD27F", alpha=0.35, lw=0)
    d.plot([0.36], [0], "*", color="#0072B2", ms=9)
    d.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"$\Delta f_{\rm SAW}/f_{\rm SAW}(\phi=0.5)$ (%)"); d.legend(frameon=False, fontsize=6)
    mf.style(d, "(D) vs composition")
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / "fig3b_options.png", dpi=200)
    print("wrote", OUT / "fig3b_options.png")


if __name__ == "__main__":
    main()


def option_b2():
    """B refined: Fe-rich (x <= 0.5) solid circles, Cr-rich dashed squares, so mirror compositions with nearly equal
    stiffening stay distinguishable; five-alloy s.d. as error bars; every curve labelled at its end."""
    dec = mf.md_decomposition(); acta = mf.md(); xs = sorted(dec)
    cmap = plt.get_cmap("viridis"); norm = plt.Normalize(0, 1)
    rel = lambda s: 100 * (s["f"] / s["f"][np.argmin(s["phi"])] - 1)
    fig, b = plt.subplots(figsize=(3.6, 3.0))
    ends = []
    for x in xs:
        s = dec[x]; o = np.argsort(s["phi"]); y = rel(s)[o]
        sd = 100 * s["f_sd"][o] / s["f"][np.argmin(s["phi"])] if s.get("f_sd") is not None else None
        fe = x <= 0.5 + 1e-9
        b.errorbar(s["phi"][o], y, yerr=sd, fmt="-o" if fe else "--s", color=cmap(norm(x)), ms=3, lw=0.9,
                   mfc=cmap(norm(x)) if fe else "white", elinewidth=0.5, capsize=1)
        ends.append((y[-1], x))
    ra = 100 * (acta["f"] / acta["f"][np.argmin(acta["phi"])] - 1)
    b.plot(acta["phi"], ra, "-", color="black", lw=1.4, zorder=0, alpha=0.35)
    b.plot(acta["phi"], ra, "o", mfc="none", mec="black", ms=6, mew=1.0, label=r"Fe$_{50}$Cr$_{50}$, Acta 50$^3$")
    b.plot([], [], "-o", color="0.4", ms=3, label=r"Fe-rich, 25$^3$ ($x_{\rm Cr}\leq0.5$)")
    b.plot([], [], "--s", color="0.4", mfc="white", ms=3, label=r"Cr-rich, 25$^3$")
    ends.sort(); placed = []
    for yv, x in ends:
        yl = yv if not placed or yv - placed[-1] > 0.28 else placed[-1] + 0.28
        placed.append(yl)
        b.text(1.015, yl, f"{round(100 * x)}", color=cmap(norm(x)), fontsize=6, va="center", fontweight="bold")
    b.text(1.015, max(placed) + 0.45, r"% Cr", fontsize=6, color="0.3")
    b.set(xlabel=r"Spinodal conversion fraction, $\phi$", ylabel=r"$\Delta f_{\rm SAW}/f_{\rm SAW}(\phi=0.5)$ (%)",
          xlim=(0.48, 1.09), ylim=(-0.3, 7.6))
    b.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    b.legend(frameon=False, fontsize=6, loc="upper left")
    mf.style(b, "(b)")
    fig.tight_layout()
    fig.savefig(OUT / "fig3b_option_B2.png", dpi=250)
    print("wrote", OUT / "fig3b_option_B2.png")


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "b2":
    option_b2()


def option_a2_d2():
    """A refined: absolute f_SAW vs phi, x_Cr as a colour bar (no end labels). D refined: stiffening vs x_Cr at every
    phi (0.6-1.0), phi as a colour bar; spinodal window (CALPHAD 748 K) shaded, Fe-36Cr marked."""
    dec = mf.md_decomposition(); acta = mf.md(); xs = sorted(dec)
    fig, (a, d) = plt.subplots(1, 2, figsize=(7.4, 3.0))
    cx = plt.get_cmap("viridis"); nx = plt.Normalize(0, 1)
    for x in xs:
        s = dec[x]; o = np.argsort(s["phi"])
        a.errorbar(s["phi"][o], s["f"][o], yerr=s["f_sd"][o] if s.get("f_sd") is not None else None, fmt="-s",
                   color=cx(nx(x)), ms=3, lw=0.9, elinewidth=0.5, capsize=1)
    a.plot(acta["phi"], acta["f"], "-o", color="black", ms=4, lw=0.9, label=r"Fe$_{50}$Cr$_{50}$, Acta 50$^3$")
    a.plot([], [], "-s", color="0.4", ms=3, label=r"present, 25$^3$ (5 alloys)")
    a.set(xlabel=r"Spinodal conversion fraction, $\phi$", ylabel=r"$f_{\rm SAW}$ (GHz), $\lambda$ = 6.4 $\mu$m", xlim=(0.48, 1.02))
    a.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    a.legend(frameon=False, fontsize=6, loc="upper left")
    fig.colorbar(plt.cm.ScalarMappable(cmap=cx, norm=nx), ax=a, pad=0.02, fraction=0.05).set_label(r"$x_{\rm Cr}$")
    mf.style(a, "(A2) absolute, colour bar")
    cp = plt.get_cmap("plasma"); npn = plt.Normalize(0.55, 1.0)
    phis = [0.6, 0.7, 0.8, 0.9, 1.0]
    for p in phis:
        y, e = [], []
        for x in xs:
            s = dec[x]; i = int(np.argmin(abs(s["phi"] - p))); i0 = int(np.argmin(s["phi"]))
            y.append(100 * (s["f"][i] / s["f"][i0] - 1))
            e.append(100 * s["f_sd"][i] / s["f"][i0] if s.get("f_sd") is not None else 0)
        d.errorbar(xs, y, yerr=e, fmt="-o", color=cp(npn(p)), ms=3, lw=0.9, elinewidth=0.5, capsize=1)
    ra = 100 * (acta["f"] / acta["f"][np.argmin(acta["phi"])] - 1)
    for p, r in zip(acta["phi"], ra):
        if p > 0.55:
            d.plot([0.5], [r], "o", mfc="none", mec=cp(npn(p)), ms=6, mew=1.0)
    d.plot([], [], "o", mfc="none", mec="0.3", ms=6, ls="none", label=r"Fe$_{50}$Cr$_{50}$, Acta 50$^3$")
    d.plot([], [], "-o", color="0.3", ms=3, label=r"present, 25$^3$ (5 alloys)")
    d.axvspan(0.30, 0.71, color="#FFD27F", alpha=0.35, lw=0, label="spinodal, 748 K")
    d.plot([0.36], [-0.25], "*", color="#0072B2", ms=8, clip_on=False)
    d.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"$\Delta f_{\rm SAW}/f_{\rm SAW}(\phi=0.5)$ (%)", xlim=(0, 1), ylim=(-0.4, 7.9))
    d.legend(frameon=False, fontsize=6, loc="upper right", bbox_to_anchor=(1.0, 1.0))
    fig.colorbar(plt.cm.ScalarMappable(cmap=cp, norm=npn), ax=d, pad=0.02, fraction=0.05, ticks=phis).set_label(r"$\phi$")
    mf.style(d, "(D2) all $\\phi$")
    fig.tight_layout()
    fig.savefig(OUT / "fig3b_options_A2_D2.png", dpi=250)
    print("wrote", OUT / "fig3b_options_A2_D2.png")


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "a2d2":
    option_a2_d2()
