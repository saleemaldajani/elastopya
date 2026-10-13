#!/usr/bin/env python3
"""Fig. 5 (draft): from elastic stiffening to embrittlement, after Pugh, Phil. Mag. 45, 823 (1954).

Pugh's relations for pure metals: flow resistance ~ G b (dislocation energy G b^2, stored mostly in the Hookean far
field), fracture strength ~ K a (the non-Hookean peak of the stress-separation curve), plastic range ~ K/G
(strain at fracture e_F ~ (K/G)^2). Spinodal decomposition raises G at nearly constant K, so K/G falls; the TGS
frequency measures G through v_R ~ sqrt(G/rho).
 (a) schematic: stress vs atomic separation (Hookean slope, non-Hookean peak) and the dislocation energy G b^2
 (b) Pugh-Pettifor map G/K vs Cauchy pressure C12 - C44 for every tensor (literature EMTO-CPA, VASP, QE, MD),
     with the MD random -> fully decomposed trajectories at 20, 36 and 50 % Cr
 (c) K/G vs x_Cr with the empirical ductile/brittle boundary K/G = 1.75
 (d) the chain for Fe-36Cr: d f/f -> d G/G -> d(K/G)/(K/G) -> d e_F/e_F, MD values and the electronic lever rule
Output: manuscript/figures/kdense_scripta/fig5_embrittlement_v0.{png,pdf,svg}
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_figures as mf  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.lines  # noqa: E402,F401

ROOT = HERE.parent
SER = {"lit": ("0.12", "o", "EMTO-CPA (Razumovskiy)"), "dft": ("#D55E00", "D", "VASP (this work)"),
       "qe": ("#8F3A00", "D", "QE (this work)"), "md": ("#009E73", "s", "MD random")}
DEC = "#0B5D3B"


def vrh(c11, c12, c44):
    b = (c11 + 2 * c12) / 3; cp = (c11 - c12) / 2
    g = ((2 * cp + 3 * c44) / 5 + 5 * (2 * cp) * c44 / (4 * c44 + 3 * (2 * cp))) / 2
    return b, g


def tensors():
    out = {k: [] for k in SER}
    for tag, (x, c11, c12, c44) in mf.saw_tensors().items():
        pre = tag.split("_x")[0]
        if pre in out:
            b, g = vrh(c11, c12, c44)
            out[pre].append((x, b, g, c12 - c44))
    return {k: np.array(sorted(v)) for k, v in out.items() if v}


def md_decomp():
    dd = mf.read_csv("md_decomposition.csv")
    res = {}
    for x0 in (0.2, 0.36, 0.5):
        m = np.abs(dd["x_cr"] - x0) < 1e-9
        o = np.argsort(dd["phi"][m])
        rows = [(ph, *vrh(c11, c12, c44), c12 - c44) for ph, c11, c12, c44 in
                zip(dd["phi"][m][o], dd["c11"][m][o], dd["c12"][m][o], dd["c44"][m][o])]
        res[x0] = np.array(rows)          # phi, K, G, Cauchy
    return res


def panel_a(ax):
    r = np.linspace(0.85, 2.2, 400)
    s = np.where(r < 1, 6 * (r - 1), 6 * (r - 1) * np.exp(-2.2 * (r - 1)))      # schematic stress-separation curve
    ax.plot(r, s, color="0.1", lw=1.3)
    ax.plot([0.9, 1.25], [6 * (0.9 - 1), 6 * 0.25], "--", color="#2A9DF4", lw=1.0)
    rp = 1 + 1 / 2.2; sp = 6 / 2.2 * np.exp(-1)
    ax.plot([rp], [sp], "o", color="#FF3B00", ms=4)
    ax.axvspan(0.97, 1.08, color="#2A9DF4", alpha=0.12, lw=0)
    ax.axvspan(1.25, 2.2, color="#FF3B00", alpha=0.07, lw=0)
    ax.text(1.02, -0.9, "Hookean:\nslope $\\propto K$\n(and $G$ in shear)", fontsize=5.4, color="#1B6FB5", ha="left", va="top")
    ax.text(rp + 0.03, sp + 0.08, "non-Hookean peak:\nfracture $\\sigma_F\\propto K a$", fontsize=5.4, color="#C02F00", va="bottom")
    ax.text(2.17, 2.5, "dislocation energy $\\propto G b^2$,\nmostly in the Hookean far field\n$\\Rightarrow$ flow stress $\\propto G b$", fontsize=5.2,
            color="0.2", va="top", ha="right")
    ax.axhline(0, color="0.6", lw=0.5)
    ax.set(xlim=(0.85, 2.2), ylim=(-1.7, 2.6), xticks=[], yticks=[], xlabel="atomic separation $r/r_0$", ylabel="stress")
    ax.text(0.99, 0.02, "schematic", transform=ax.transAxes, ha="right", fontsize=4.8, style="italic", color="0.5")


def panel_b(ax, T, D):
    ax.axvspan(-20, 0, color="0.93", lw=0); ax.axhspan(1 / 1.75, 0.75, color="0.93", lw=0)
    ax.text(-18, 0.36, "brittle\n(Pettifor,\nPugh)", fontsize=5.0, color="0.4", va="bottom")
    for k, (col, mk, lab) in SER.items():
        if k in T:
            a = T[k]
            ax.plot(a[:, 3], a[:, 2] / a[:, 1], mk, color=col, mfc="white" if k == "qe" else col, mec=col, ms=3.5, ls="none", label=lab)
    for x0, a in D.items():
        r, d = a[0], a[-1]
        ax.annotate("", xy=(d[3], d[2] / d[1]), xytext=(r[3], r[2] / r[1]),
                    arrowprops=dict(arrowstyle="-|>", color="#FF3B00", lw=1.0, mutation_scale=7))
        ax.plot([d[3]], [d[2] / d[1]], "s", mfc="white", mec=DEC, ms=4, mew=1)
    ax.plot([], [], "s", mfc="white", mec=DEC, ms=4, label="MD decomposed ($\\phi$=1)")
    ax.set(xlabel=r"Cauchy pressure $C_{12}-C_{44}$ (GPa)", ylabel=r"$G/K$", xlim=(-20, 70), ylim=(0.35, 0.75))
    ax.legend(fontsize=4.8, loc="center right", bbox_to_anchor=(1.0, 0.74), frameon=False)


def panel_c(ax, T, D):
    ax.axhspan(1.0, 1.75, color="0.93", lw=0)
    ax.axhline(1.75, color="0.4", lw=0.8, ls="--"); ax.text(0.98, 1.72, "brittle", ha="right", va="top", fontsize=5.4, color="0.4")
    ax.text(0.98, 1.78, "ductile", ha="right", va="bottom", fontsize=5.4, color="0.4")
    for k, (col, mk, lab) in SER.items():
        if k in T:
            a = T[k]
            ax.plot(a[:, 0], a[:, 1] / a[:, 2], mk + ("-" if k in ("lit", "md") else ""), color=col, mfc="white" if k == "qe" else col,
                    mec=col, ms=3.2, lw=0.7)
    for x0, a in D.items():
        ax.annotate("", xy=(x0, a[-1, 1] / a[-1, 2]), xytext=(x0, a[0, 1] / a[0, 2]),
                    arrowprops=dict(arrowstyle="-|>", color="#FF3B00", lw=1.0, mutation_scale=7))
    ax.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"Pugh ratio $K/G$", xlim=(0, 1), ylim=(1.4, 2.8))


def panel_d(ax, D, legend_at=(1.0, 0.47), legend_loc="center right"):
    a = D[0.36]
    r, d = a[0], a[-1]
    dG = d[2] / r[2] - 1; dK = d[1] / r[1] - 1; dKG = (d[1] / d[2]) / (r[1] / r[2]) - 1
    dec = mf.md_decomposition()[0.36]
    df = float(dec["f"][np.argmax(dec["phi"])] / dec["f"][np.argmin(dec["phi"])]) - 1
    deF = (1 + dKG) ** 2 - 1
    lever_e = 0.077                                         # electronic lever rule (Fig. 4e)
    rows = [(r"$f$ (TGS)", df, lever_e), (r"$G$", dG, 2 * lever_e), (r"$K$", dK, None),
            (r"$K/G$", dKG, (1 + dK) / (1 + 2 * lever_e) - 1), (r"$\varepsilon_F\propto(K/G)^2$", deF,
                                                                               ((1 + dK) / (1 + 2 * lever_e)) ** 2 - 1)]
    y = np.arange(len(rows))[::-1]
    for yi, (lab, v, ve) in zip(y, rows):
        ax.barh(yi + 0.17, 100 * v, height=0.32, color="#009E73")
        ax.text(100 * v + (1 if v >= 0 else -1), yi + 0.17, f"{100 * v:+.1f}%", va="center", ha="left" if v >= 0 else "right", fontsize=5.3)
        if ve is not None:
            ax.barh(yi - 0.17, 100 * ve, height=0.32, color="0.15")
            ax.text(100 * ve + (1 if ve >= 0 else -1), yi - 0.17, f"{100 * ve:+.1f}%", va="center", ha="left" if ve >= 0 else "right", fontsize=5.3)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=5.6)
    ax.axvline(0, color="0.5", lw=0.6)
    ax.set(xlim=(-45, 32), xlabel="change on full decomposition, Fe$-$36Cr (%)")
    from matplotlib.patches import Patch
    short = legend_loc != "center right"
    ax.legend(handles=[Patch(color="#009E73", label="MD" if short else "MD (computed tensors)"),
                       Patch(color="0.15", label="electronic (est.)" if short else "electronic lever rule (estimate)")],
              fontsize=4.8, loc=legend_loc, bbox_to_anchor=legend_at, frameon=False, handlelength=1.0)
    ax.text(0.02, 0.97, r"$\Delta(K/G)/(K/G)\approx-2\,\Delta f/f$" "\n($K$ nearly unchanged)", transform=ax.transAxes, ha="left", va="top",
            fontsize=5.4, color="#C02F00")


def main():
    sys.path.insert(0, str(HERE / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    T, D = tensors(), md_decomp()
    fig, axs = plt.subplots(1, 4, figsize=(190 / 25.4, 2.45), gridspec_kw={"width_ratios": [1.05, 1, 1, 1.15]})
    panel_a(axs[0]); panel_b(axs[1], T, D); panel_c(axs[2], T, D); panel_d(axs[3], D)
    for ax, (l, t) in zip(axs, (("a", "Pugh 1954"), ("b", "bonding map"), ("c", "plastic range"), ("d", "TGS shift to fracture strain"))):
        ax.set_title(f"({l}) {t}", fontsize=7, loc="left")
        ax.tick_params(labelsize=5.5); ax.xaxis.label.set_size(6); ax.yaxis.label.set_size(6)
    fig.tight_layout(w_pad=0.6)
    out = ROOT / "manuscript" / "figures" / "kdense_scripta" / "fig5_embrittlement_v0"
    fig.savefig(f"{out}.png", dpi=600); fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.svg")
    print("wrote", out)
    a = D[0.36]
    print("MD Fe-36Cr random -> decomposed: K", a[0, 1].round(1), a[-1, 1].round(1), " G", a[0, 2].round(1), a[-1, 2].round(1),
          " K/G", (a[0, 1] / a[0, 2]).round(3), (a[-1, 1] / a[-1, 2]).round(3))




# ----------------------------------------------------------------------------- v1: atomistic mechanisms (top row)
FE_C, CR_C = "#F06060", "#6060F0"


def _lattice(ax, nx, ny, xfrac, seed=3, s=26):
    rng = np.random.default_rng(seed)
    pts = np.array([(i, j) for j in range(ny) for i in range(nx)], float)
    el = rng.random(len(pts)) < xfrac
    return pts, el


def atom_a(ax):
    """Two atoms and their bond: Hookean stretch vs bond breaking, beside Pugh's stress-separation curve."""
    r = np.linspace(0.85, 2.2, 300)
    sig = np.where(r < 1, 6 * (r - 1), 6 * (r - 1) * np.exp(-2.2 * (r - 1)))
    ax.plot(r, sig, color="0.1", lw=1.2)
    ax.plot([0.9, 1.25], [-0.6, 1.5], "--", color="#2A9DF4", lw=0.9)
    rp = 1 + 1 / 2.2; sp = 6 / 2.2 * np.exp(-1)
    ax.plot([rp], [sp], "o", color="#FF3B00", ms=3.5)
    for (x0, y0, d, col, lab) in ((1.08, -1.25, 0.10, "#1B6FB5", "stretch\n($K$, $G$)"),
                                  (1.85, -1.25, 0.22, "#C02F00", "break\n($\\sigma_F\\propto Ka$)")):
        ax.add_patch(plt.Circle((x0 - d, y0), 0.06, color=FE_C, ec="0.2", lw=0.4, zorder=4))
        ax.add_patch(plt.Circle((x0 + d, y0), 0.06, color=CR_C, ec="0.2", lw=0.4, zorder=4))
        ax.plot([x0 - d, x0 + d], [y0, y0], color=col, lw=1.6 if d < 0.15 else 0.8, ls="-" if d < 0.15 else (0, (1.5, 1.5)), zorder=3)
        ax.text(x0, y0 - 0.22, lab, ha="center", va="top", fontsize=5.0, color=col)
    ax.axhline(0, color="0.7", lw=0.5)
    ax.set(xlim=(0.82, 2.2), ylim=(-2.3, 2.4), xticks=[], yticks=[], xlabel="atomic separation $r/r_0$", ylabel="stress")


def atom_b(ax):
    """Edge dislocation in a mixed Fe/Cr lattice (isotropic displacement field): non-Hookean core vs Hookean far field."""
    nu = 0.29
    pts, el = _lattice(ax, 15, 11, 0.36)
    x, y = pts[:, 0] - 7.0, pts[:, 1] - 5.5
    keep = ~((x > 0.01) & (y < 0)) | (np.abs(x) > 99)                    # extra half-plane above: remove x>0 below? keep all
    th = np.arctan2(y, x); r2 = x ** 2 + y ** 2 + 1e-9
    ux = (1 / (2 * np.pi)) * (th + x * y / (2 * (1 - nu) * r2))
    uy = -(1 / (2 * np.pi)) * ((1 - 2 * nu) / (4 * (1 - nu)) * np.log(r2) + (x ** 2 - y ** 2) / (4 * (1 - nu) * r2))
    X, Y = x + ux, y + 0.6 * uy
    ax.add_patch(plt.Circle((0, 0), 1.6, color="#FF3B00", alpha=0.15, lw=0))
    ax.add_patch(plt.Circle((0, 0), 5.6, fill=False, ec="#2A9DF4", lw=0.8, ls="--"))
    ax.scatter(X, Y, s=11, c=[CR_C if e else FE_C for e in el], edgecolors="0.25", linewidths=0.25, zorder=3)
    ax.text(0, -0.15, "$\\perp$", ha="center", va="center", fontsize=9, fontweight="bold", zorder=5)
    BB = dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85)
    ax.text(0, -6.0, "core (red): non-Hookean", ha="center", va="top", fontsize=4.8, color="#C02F00", bbox=BB)
    ax.text(0, 5.9, "far field: Hookean, energy $\\propto Gb^2\\ln(R/r_0)$", ha="center", va="bottom", fontsize=4.8, color="#1B6FB5", bbox=BB)
    ax.set_xlim(-7.6, 7.6); ax.set_ylim(-7.4, 7.4); ax.set_aspect("equal"); ax.axis("off")


def atom_c(ax):
    """Crack tip: cleave the tip bond (~ K a) or emit a dislocation on the slip plane (~ G b)."""
    pts, el = _lattice(ax, 15, 11, 0.36, seed=5)
    x, y = pts[:, 0] - 4.0, pts[:, 1] - 5.0
    crack = (x < 0) & (np.abs(y) < 0.6 + 0.12 * (-x))
    x2, y2 = x[~crack], y[~crack]
    open_ = np.where(x2 < 0.5, np.sign(y2) * 0.35 * np.sqrt(np.clip(0.5 - x2, 0, None)), 0)
    ax.scatter(x2, y2 + open_, s=11, c=[CR_C if e else FE_C for e in el[~crack]], edgecolors="0.25", linewidths=0.25, zorder=3)
    ax.plot([0, 0], [0.62, -0.62], color="#FF3B00", lw=1.0, ls=(0, (1.2, 1.2)), zorder=4)
    ax.annotate("", xy=(3.2, 0), xytext=(0.6, 0), arrowprops=dict(arrowstyle="-|>", color="#C02F00", lw=1.0))
    BB = dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.9)
    ax.text(3.3, 0.15, "cleave ($\\propto Ka$)", fontsize=4.8, color="#C02F00", va="bottom", bbox=BB)
    ax.annotate("", xy=(2.4, 3.2), xytext=(0.3, 0.5), arrowprops=dict(arrowstyle="-|>", color="#1B6FB5", lw=1.0))
    ax.plot([0.3, 3.0], [0.5, 4.0], color="#1B6FB5", lw=0.5, ls="--")
    ax.text(2.5, 4.1, "emit $\\perp$ ($\\propto Gb$)", fontsize=4.8, color="#1B6FB5", va="bottom", bbox=BB)
    ax.text(3.0, -6.4, "$\\alpha'$ stiffer in shear ($G\\uparrow$): emission harder,\ncleavage favoured", fontsize=4.8,
            color="0.15", va="top", ha="center")
    ax.set_xlim(-4.3, 10.3); ax.set_ylim(-8.0, 6.2); ax.set_aspect("equal"); ax.axis("off")


def atom_d(ax):
    """Dislocation pinned by Cr-rich alpha' domains (real decomposed MD plane, Fe-36Cr, phi = 1)."""
    from make_figures_kdense import md_grid
    from visualize_decomposition import series
    from matplotlib.colors import ListedColormap
    g = md_grid(series()[0.36][-1])
    ax.imshow(g, cmap=ListedColormap([FE_C, CR_C]), origin="lower", interpolation="nearest", alpha=0.85)
    n = g.shape[0]
    yy = np.linspace(0, n - 1, 200)
    xx = 11 + 2.2 * np.sin(np.pi * yy / (n - 1) * 3) ** 2          # line bowing between obstacles
    ax.plot(xx, yy, color="black", lw=2.2); ax.plot(xx, yy, color="#FFD400", lw=1.0)
    ax.annotate("", xy=(17, n * 0.5), xytext=(13.5, n * 0.5), arrowprops=dict(arrowstyle="-|>", color="black", lw=1.0))
    ax.text(17.3, n * 0.5, "$\\tau$", fontsize=7, va="center")
    T = tensors()
    pc = mf.read_csv("md_phase_compositions.csv"); m = (np.abs(pc["x_cr"] - 0.36) < 1e-6) & (pc["phi"] == 1.0)
    xa, xb = float(pc["x_alpha"][m][0]), float(pc["x_alpha_prime"][m][0])
    gdiff = {k: float(np.interp(xb, T[k][:, 0], T[k][:, 2]) - np.interp(xa, T[k][:, 0], T[k][:, 2])) for k in ("lit", "md")}
    eta = 100 * (xb - xa) * (2.8846 - 2.8665) / 2.8665
    dg = gdiff["lit"]
    ax.text(0.5, -1.0, f"$G_{{\\alpha'}}-G_\\alpha$ = {gdiff['lit']:+.0f} GPa (EMTO-CPA),\n{gdiff['md']:+.0f} GPa (MD); lattice misfit {eta:.2f}%\n"
            f"$\\Rightarrow$ modulus + misfit pinning: hardening",
            fontsize=4.7, va="top", ha="left", color="0.15")
    ax.set_xlim(-0.5, n - 0.5); ax.set_ylim(-6.5, n - 0.5); ax.axis("off")
    return dg


def tgs_summary(ax):
    """Pristine vs aged Fe-36Cr: real MD planes under a TGS spot, the measured peak shift, and the consequence
    for the stress-strain curve (Pugh: higher flow stress ~ Gb, shorter plastic range ~ (K/G)^2)."""
    from make_figures_kdense import md_grid
    from visualize_decomposition import series
    from matplotlib.colors import ListedColormap
    import matplotlib.patches as mpa
    ser = series()[0.36]
    cm = ListedColormap([FE_C, CR_C])
    dec = mf.md_decomposition()[0.36]
    f0 = float(dec["f"][np.argmin(dec["phi"])]) * 1e3; f1 = float(dec["f"][np.argmax(dec["phi"])]) * 1e3
    for k, (x0, grid, lab, f, col) in enumerate(((0.03, md_grid(ser[0]), "pristine", f0, "#009E73"),
                                                  (0.55, md_grid(ser[-1]), "aged 475 °C", f1, "#0B5D3B"))):
        w = 0.42
        ax.imshow(grid, extent=(x0, x0 + w, 0.40, 0.66), cmap=cm, interpolation="nearest", aspect="auto", alpha=0.9, zorder=1)
        ax.add_patch(mpa.Rectangle((x0, 0.40), w, 0.26, fill=False, ec="0.2", lw=0.6, zorder=2))
        xc = x0 + w / 2
        for side in (-1, 1):                                  # crossed pump beams
            ax.plot([xc + side * 0.17, xc], [0.86, 0.665], color="#FFB000", lw=1.6, alpha=0.75, zorder=3)
        spot = mpa.Ellipse((xc, 0.66), 0.22, 0.05, fc="#FFB000", ec="#FF8C1A", lw=0.6, alpha=0.55, zorder=4)
        ax.add_patch(spot)
        for xf in np.linspace(xc - 0.08, xc + 0.08, 5):        # grating fringes in the spot
            ax.plot([xf, xf], [0.645, 0.675], color="#C02F00", lw=0.5, zorder=5)
        fx = np.linspace(0, 1, 100)                           # measured peak, frequency axis shared
        fr = 430 + 60 * fx
        pk = np.exp(-0.5 * ((fr - f) / 4.0) ** 2)
        ax.plot(x0 + w * fx, 0.885 + 0.075 * pk, color=col, lw=1.0, zorder=5)
        ax.plot([x0, x0 + w], [0.885, 0.885], color="0.5", lw=0.4)
        ax.text(x0 + w * (f - 430) / 60, 0.97, f"{f:.0f} MHz", ha="center", va="bottom", fontsize=4.9, color=col)
        ax.text(xc, 0.37, lab, ha="center", va="top", fontsize=5.4, fontweight="bold")
        # stress-strain sketch: aged has higher flow stress (~ G b) and shorter plastic range (~ (K/G)^2)
        eps = np.linspace(0, 1, 200)
        sy, ef = (0.10, 0.95) if k == 0 else (0.10 * (1 + 0.111), 0.95 * (1 - 0.209))
        e = eps[eps <= ef]
        sig = np.where(e < 0.08, e / 0.08 * sy, sy + 0.05 * (e - 0.08) ** 0.5)
        ax.plot(x0 + w * e, 0.04 + sig * 1.25, color=col, lw=1.1)
        ax.plot([x0, x0 + w], [0.04, 0.04], color="0.5", lw=0.4); ax.plot([x0, x0], [0.04, 0.27], color="0.5", lw=0.4)
        ax.text(x0 - 0.01, 0.27, "$\\sigma$", fontsize=4.8, ha="right", va="top", color="0.35")
        ax.text(x0 + w, 0.035, "$\\varepsilon$", fontsize=4.8, ha="right", va="top", color="0.35")
        ax.plot([x0 + w * e[-1]], [0.04 + sig[-1] * 1.25], "x", color="#C02F00" if k else col, ms=4, mew=1.0)
        if k:
            ax.text(xc, 0.0, "flow +11%,  $\\varepsilon_F$ $-$21%", fontsize=4.9, color="#C02F00", va="top", ha="center")
        else:
            ax.text(xc, 0.0, "ductile", fontsize=4.9, color=col, va="top", ha="center")
    ax.annotate("", xy=(0.56, 0.53), xytext=(0.445, 0.53), arrowprops=dict(arrowstyle="-|>", color="#FF3B00", lw=1.2), zorder=6)
    ax.text(0.5, 0.995, "TGS peak", ha="center", va="top", fontsize=4.9, color="0.35")
    ax.set_xlim(0, 1); ax.set_ylim(-0.07, 1.04); ax.axis("off")


def main_v1():
    sys.path.insert(0, str(HERE / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    T, D = tensors(), md_decomp()
    fig = plt.figure(figsize=(190 / 25.4, 4.9))
    gt = fig.add_gridspec(1, 4, left=0.045, right=0.995, top=0.95, bottom=0.56, wspace=0.12, width_ratios=[1.05, 1, 1, 1])
    gb = fig.add_gridspec(1, 4, left=0.07, right=0.995, top=0.44, bottom=0.09, wspace=0.42, width_ratios=[1, 1, 1.05, 0.95])
    top = [fig.add_subplot(gt[0, i]) for i in range(4)]
    bot = [fig.add_subplot(gb[0, i]) for i in range(4)]
    atom_a(top[0]); atom_b(top[1]); atom_c(top[2]); dg = atom_d(top[3])
    for ax, t in zip(top, ("(a) bonds: Hookean vs breaking", "(b) dislocation: core vs far field",
                           "(c) crack tip: cleave or emit", "(d) $\\alpha'$ domains pin dislocations")):
        ax.set_title(t, fontsize=7, loc="left")
    panel_b(bot[0], T, D); panel_c(bot[1], T, D); panel_d(bot[2], D)
    tgs_summary(bot[3])
    for ax, t in zip(bot[:3], ("(e) bonding map", "(f) plastic range", "(g) TGS shift to fracture strain")):
        ax.set_title(t, fontsize=7, loc="left")
        ax.tick_params(labelsize=5.5); ax.xaxis.label.set_size(6); ax.yaxis.label.set_size(6)
    bot[3].set_title("(h) what TGS sees", fontsize=7, loc="left")
    fig.text(0.995, 0.515, "atoms and defects (schematic; MD plane real)", ha="right", fontsize=5, style="italic", color="0.5")
    out = ROOT / "manuscript" / "figures" / "kdense_scripta" / "fig5_embrittlement_v1"
    fig.savefig(f"{out}.png", dpi=600); fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.svg")
    print("wrote", out, "dG alpha'-alpha =", round(dg, 1))

def dft_cells(fig, rect):
    """Atomic (DFT) view: the relaxed 40-atom QE SQS of an alpha-like and an alpha'-like composition with their QE
    moments, and their computed shear modulus and Pugh ratio."""
    import visualize_sqs as VS
    T = tensors()
    qe = {round(float(r[0]), 2): r for r in T["qe"]}
    w = rect[2] / 2 - 0.004
    for i, x in enumerate((0.20, 0.95)):
        ax = fig.add_axes([rect[0] + i * (w + 0.008), rect[1] + rect[3] * 0.30, w, rect[3] * 0.55])
        VS.draw(ax, x); ax.set_aspect("equal", adjustable="datalim")
        r = qe[x]
        ax.text(0.5, -0.12, f"{round(100 * x)}% Cr SQS (QE)\n$G$ = {r[2]:.0f} GPa, $K/G$ = {r[1] / r[2]:.2f}", transform=ax.transAxes,
                ha="center", va="top", fontsize=5.0)
    fig.text(rect[0], rect[1] + rect[3] * 0.92, "40 atoms; arrows: QE local moments", fontsize=5.0, color="0.35")


def md_cubes(fig, rect):
    """Molecular (MD) view: the 25^3 Fe-36Cr cells, random (phi = 0.5) and fully decomposed (phi = 1)."""
    from visualize_decomposition import _cube_ax, md_grid3d, series
    ser = series()[0.36]
    w = rect[2] / 2 - 0.004
    for i, (f, lab) in enumerate(((ser[0], "random ($\\phi$ = 0.5)"), (ser[-1], "decomposed ($\\phi$ = 1)"))):
        ax = fig.add_axes([rect[0] + i * (w + 0.008), rect[1] + rect[3] * 0.18, w, rect[3] * 0.70])
        _cube_ax(ax, md_grid3d(f))
        ax.text(0.5, -0.04, lab, transform=ax.transAxes, ha="center", va="top", fontsize=5.0)
    fig.text(rect[0], rect[1] + rect[3] * 0.92, "microstructure: 31,250 atoms (25$^3$ bcc), EAM", fontsize=5.0, color="0.35")


def md_phi(ax):
    """MD Fe-36Cr: G, K (VRH) and K/G versus conversion fraction phi."""
    a = md_decomp()[0.36]
    ph, K, G = a[:, 0], a[:, 1], a[:, 2]
    ax.plot(ph, 100 * (G / G[0] - 1), "-s", color="#0B5D3B", ms=3, lw=0.9, label="$G$")
    ax.plot(ph, 100 * (K / K[0] - 1), "-o", color="0.45", ms=3, lw=0.9, label="$K$")
    ax.plot(ph, 100 * ((K / G) / (K[0] / G[0]) - 1), "-^", color="#C02F00", ms=3, lw=0.9, label="$K/G$")
    ax.axhline(0, color="0.6", lw=0.5)
    ax.set(xlabel="spinodal conversion $\\phi$", ylabel="change vs random (%)", xlim=(0.48, 1.02))
    ax.legend(fontsize=5, frameon=False, loc="center left")


def main_v2():
    """v1 plus an explicit atomic (DFT) row and molecular (MD) row."""
    sys.path.insert(0, str(HERE / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    T, D = tensors(), md_decomp()
    W, H = 190 / 25.4, 7.2
    fig = plt.figure(figsize=(W, H))
    g1 = fig.add_gridspec(1, 3, left=0.055, right=0.745, top=0.955, bottom=0.715, wspace=0.12, width_ratios=[1.05, 1, 1])
    g3 = fig.add_gridspec(1, 4, left=0.085, right=0.995, top=0.30, bottom=0.06, wspace=0.42, width_ratios=[1, 1, 1.05, 0.95])
    a = [fig.add_subplot(g1[0, i]) for i in range(3)]
    atom_a(a[0]); atom_b(a[1]); atom_c(a[2])
    for ax, t in zip(a, ("(a) bonds: Hookean vs breaking", "(b) dislocation: core vs far field", "(c) crack tip: cleave or emit")):
        ax.set_title(t, fontsize=7, loc="left")
    dft_cells(fig, [0.755, 0.70, 0.24, 0.27])
    fig.text(0.757, 0.962, "(d) DFT cells: $\\alpha$-like vs $\\alpha'$-like", fontsize=7)
    # molecular row
    md_cubes(fig, [0.045, 0.375, 0.29, 0.27])
    fig.text(0.05, 0.66, "(e) MD cells: Fe$-$36Cr before and after decomposition", fontsize=7)
    pin = fig.add_axes([0.37, 0.38, 0.17, 0.25]); dg = atom_d(pin)
    fig.text(0.37, 0.66, "(f) $\\alpha'$ domains pin dislocations", fontsize=7)
    mp = fig.add_axes([0.62, 0.40, 0.35, 0.22]); md_phi(mp)
    fig.text(0.62, 0.66, "(g) MD: $G$ rises, $K$ flat, $K/G$ falls with $\\phi$", fontsize=7)
    mp.tick_params(labelsize=5.5); mp.xaxis.label.set_size(6); mp.yaxis.label.set_size(6)
    for y, lab in ((0.84, "ATOMIC (DFT)"), (0.52, "MOLECULAR (MD)"), (0.18, "DATA & TGS")):
        fig.text(0.006, y, lab, fontsize=7.5, fontweight="heavy", color="0.1", va="center", ha="left", rotation=90)
    for y in (0.695, 0.345):
        fig.add_artist(matplotlib.lines.Line2D([0.01, 0.99], [y, y], color="0.75", lw=0.6, ls=(0, (4, 3))))
    b = [fig.add_subplot(g3[0, i]) for i in range(4)]
    panel_b(b[0], T, D); panel_c(b[1], T, D); panel_d(b[2], D); tgs_summary(b[3])
    for ax, t in zip(b, ("(h) bonding map", "(i) plastic range", "(j) TGS shift to fracture strain", "(k) what TGS sees")):
        ax.set_title(t, fontsize=7, loc="left")
        if t[1] != "k":
            ax.tick_params(labelsize=5.5); ax.xaxis.label.set_size(6); ax.yaxis.label.set_size(6)
    out = ROOT / "manuscript" / "figures" / "kdense_scripta" / "fig5_embrittlement_v2"
    fig.savefig(f"{out}.png", dpi=600); fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.svg")
    print("wrote", out)


def main_v3():
    """v2 with the panel (j) legend moved clear of the bar labels (upper left, below the formula)."""
    global panel_d
    orig = panel_d
    panel_d = lambda ax, D: orig(ax, D, legend_at=(1.0, 0.37), legend_loc="right")
    try:
        import types
        src_out = ROOT / "manuscript" / "figures" / "kdense_scripta" / "fig5_embrittlement_v2"
        main_v2.__globals__["panel_d"] = panel_d
        main_v2()
    finally:
        panel_d = orig
        main_v2.__globals__["panel_d"] = orig
    for ext in ("png", "pdf", "svg"):
        v2 = ROOT / "manuscript" / "figures" / "kdense_scripta" / f"fig5_embrittlement_v2.{ext}"
        v3 = ROOT / "manuscript" / "figures" / "kdense_scripta" / f"fig5_embrittlement_v3.{ext}"
        v2.replace(v3)
    main_v2()                                   # restore v2 unchanged
    print("wrote fig5_embrittlement_v3")


ETA_MISFIT = (2.8846 - 2.8665) / 2.8665      # (1/a) da/dx, Vegard, room temperature
OMEGA_A3 = 11.8                               # atomic volume (A^3)


def young_planar(K, G):
    """Cahn's isotropic coherency modulus Y = E/(1 - nu) = 18 K G / (3K + 4G)."""
    return 18 * K * G / (3 * K + 4 * G)


def psi_split(ax, D):
    """Energy stored per unit strain: volumetric (K/2 per unit dilatation^2) vs shear (G per unit shear^2), Fe-36Cr."""
    a = D[0.36]
    r, d = a[0], a[-1]
    lever = 1.154                                   # electronic lever-rule G factor (2 x 7.7 %)
    cats = ["$K/2$\n(volumetric)", "$G$\n(shear)", "$Y$\n(coherency)"]
    rnd = [r[1] / 2, r[2], young_planar(r[1], r[2])]
    dec = [d[1] / 2, d[2], young_planar(d[1], d[2])]
    x = np.arange(3)
    ax.bar(x - 0.2, rnd, 0.38, color="#009E73", label="random (MD)")
    ax.bar(x + 0.2, dec, 0.38, color="#0B5D3B", label="decomposed (MD)")
    for xi, a0, b0 in zip(x, rnd, dec):
        ax.text(xi + 0.2, b0 + 6, f"{100 * (b0 / a0 - 1):+.1f}%", ha="center", fontsize=5.4, color="#C02F00")
    ax.set_xticks(x); ax.set_xticklabels(cats, fontsize=5.6)
    ax.set_ylabel("stiffness in $\\Psi$ (GPa)", fontsize=6)
    ax.legend(fontsize=5, frameon=False, loc="upper left")
    ax.set_ylim(0, 450)
    ax.text(0.5, 0.80, "$\\Psi=\\frac{1}{2}K(\\mathrm{tr}\\,\\varepsilon)^2+G\\,\\varepsilon':\\varepsilon'$", transform=ax.transAxes,
            ha="center", va="center", fontsize=6.0)


def coherency(ax, T):
    """Coherency strain energy 2 eta^2 Y Omega (meV/atom per unit (dx)^2) vs composition, against the chemical term."""
    ev = 6.241509e18
    for k, col, mk, lab in (("lit", "0.12", "o", "EMTO-CPA"), ("qe", "#8F3A00", "D", "QE"), ("md", "#009E73", "s", "MD")):
        if k in T:
            a = T[k]
            y = young_planar(a[:, 1], a[:, 2]) * 1e9 * 2 * ETA_MISFIT ** 2 * OMEGA_A3 * 1e-30 * ev * 1e3
            ax.plot(a[:, 0], y, mk + ("-" if k != "qe" else ""), color=col, mfc="white" if k == "qe" else col, ms=3, lw=0.7, label=lab)
    ax.set(xlim=(0, 1), ylim=(0, 3.8), xlabel=r"$x_{\rm Cr}$", ylabel=r"$2\eta^2Y\Omega$ (meV/atom)")
    ax.legend(fontsize=5, frameon=False, loc="lower right")
    ax.text(0.03, 0.97, "chemical curvature $|\\partial^2\\Delta G/\\partial x^2|\\approx$ 1400 meV/atom\n"
            f"$\\Rightarrow$ coherency shifts the spinodal\n    by ~0.1% ($\\eta$ = {100 * ETA_MISFIT:.2f}%)",
            transform=ax.transAxes, va="top", fontsize=4.6, color="0.15")


def psi_hub(ax, D):
    """How it comes together: Psi(K, G) and its five consequences, with Fe-36Cr numbers (MD, random -> decomposed)."""
    a = D[0.36]; r, d = a[0], a[-1]
    dec = mf.md_decomposition()[0.36]
    df = float(dec["f"][np.argmax(dec["phi"])] / dec["f"][np.argmin(dec["phi"])]) - 1
    dG = d[2] / r[2] - 1; dK = d[1] / r[1] - 1; dKG = (d[1] / d[2]) / (r[1] / r[2]) - 1
    dY = young_planar(d[1], d[2]) / young_planar(r[1], r[2]) - 1
    ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.05, 1.05); ax.set_aspect("equal"); ax.axis("off")
    ax.add_patch(plt.Circle((0, 0), 0.30, fc="#FFF1D6", ec="#FF8C1A", lw=1.2, zorder=3))
    ax.text(0, 0.03, "$\\Psi(K,G)$", ha="center", va="center", fontsize=8, fontweight="bold", zorder=4)
    ax.text(0, -0.13, "Hookean energy", ha="center", va="center", fontsize=4.8, zorder=4)
    spokes = [(90, "TGS / SAW", "$f=\\eta_R(\\nu)\\sqrt{G/\\rho}/\\lambda$", df),
              (18, "slip", "$E_\\perp\\propto Gb^2$", dG),
              (-54, "cleavage", "$\\sigma_F\\propto Ka$", dK),
              (-126, "ductility", "$\\varepsilon_F\\propto(K/G)^2$", (1 + dKG) ** 2 - 1),
              (162, "coherency", "$Y=\\frac{18KG}{3K+4G}$", dY)]
    for ang, name, eq, v in spokes:
        t = np.radians(ang); x, y = 0.78 * np.cos(t), 0.78 * np.sin(t)
        ax.plot([0.30 * np.cos(t), 0.66 * np.cos(t)], [0.30 * np.sin(t), 0.66 * np.sin(t)], color="#FF8C1A", lw=1.0, zorder=1)
        col = "#C02F00" if (name == "ductility") else "0.1"
        ax.text(x, y, f"{name}\n{eq}\n{100 * v:+.1f}%", ha="center", va="center", fontsize=5.0, color=col, linespacing=1.2,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.75", lw=0.5), zorder=2)


def psi_flow(ax, D):
    """(n) as a mechanism: decomposition -> which elastic channel of Psi changes (shear G vs volume K) -> what each
    channel controls, with Fe-36Cr numbers (MD, random -> fully decomposed); the dashed path is the TGS inference."""
    a = D[0.36]; r, d = a[0], a[-1]
    dec = mf.md_decomposition()[0.36]
    df = float(dec["f"][np.argmax(dec["phi"])] / dec["f"][np.argmin(dec["phi"])]) - 1
    dG = d[2] / r[2] - 1; dK = d[1] / r[1] - 1; dKG = (d[1] / d[2]) / (r[1] / r[2]) - 1
    dY = young_planar(d[1], d[2]) / young_planar(r[1], r[2]) - 1
    deF = (1 + dKG) ** 2 - 1
    SH, VO, BAD, INK = "#0B7A55", "#0072B2", "#C02F00", "0.15"
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    def box(x, y, txt, col, w=None, fc="white", bold_first=True, fs=5.0):
        t = ax.text(x, y, txt, ha="center", va="center", fontsize=fs, color=INK, linespacing=1.25, zorder=3,
                    bbox=dict(boxstyle="round,pad=0.35", fc=fc, ec=col, lw=0.9))
        return t

    def arrow(src, dst, col, ls="-", rad=0.0, lab=None, lab_xy=None, lw=0.9):
        ax.annotate("", xy=(0, 0.5), xycoords=dst, xytext=(1, 0.5), textcoords=src,
                    arrowprops=dict(arrowstyle="-|>", color=col, lw=lw, ls=ls, mutation_scale=6,
                                    shrinkA=1, shrinkB=1, connectionstyle=f"arc3,rad={rad}"), zorder=2)
        if lab:
            ax.text(*lab_xy, lab, fontsize=4.3, color=col, ha="center", va="center", style="italic")

    # column headers
    for x, h in ((0.11, "cause"), (0.40, "elastic channel of $\\Psi$"), (0.74, "what it controls")):
        ax.text(x, 0.995, h, ha="center", va="top", fontsize=5.2, color="0.4", fontweight="bold")
    # cause, and the (negligible) elastic feedback on it
    c = box(0.11, 0.62, "spinodal\ndecomposition\nFe$-$36Cr\n$\\to\\alpha+\\alpha'$", "0.35", fc="#F2F2F2")
    coh = box(0.11, 0.13, f"coherency\n$Y(K,G)$ {100 * dY:+.1f}%\n$\\approx$2 meV/atom\n$\\ll$ chemical", "0.45")
    # channels, with K/G between them
    g = box(0.40, 0.82, f"shear: $G\\,\\varepsilon':\\varepsilon'$\n$G$ {100 * dG:+.1f}%\n($\\alpha'$ stiffer in shear)", SH, fc="#E8F5EF")
    kg = box(0.40, 0.50, f"$K/G$ {100 * dKG:+.1f}%", BAD, fc="#FBEAE5")
    k = box(0.40, 0.18, f"volume: $\\frac{{1}}{{2}}K(\\mathrm{{tr}}\\,\\varepsilon)^2$\n$K$ {100 * dK:+.1f}%\n(mean $e/a$ unchanged)", VO, fc="#E7F1F9")
    arrow(c, g, SH, rad=-0.2); arrow(c, k, VO, rad=0.2)
    for src, col in ((g, SH), (k, VO)):
        ax.annotate("", xy=(0.5, 1.0 if src is g else 0.0), xycoords=kg, xytext=(0.5, 0.0 if src is g else 1.0),
                    textcoords=src, arrowprops=dict(arrowstyle="-|>", color=col, lw=0.9, mutation_scale=6, shrinkA=1, shrinkB=1))
    ax.annotate("", xy=(1.0, 0.5), xycoords=coh, xytext=(0.0, 0.5), textcoords=k,      # Y(K, G); its feedback is tiny
                arrowprops=dict(arrowstyle="-|>", color="0.45", lw=0.6, ls=(0, (2, 1.5)), mutation_scale=5, shrinkA=1, shrinkB=1))
    ax.annotate("", xy=(0.5, 0.0), xycoords=c, xytext=(0.5, 1.0), textcoords=coh,
                arrowprops=dict(arrowstyle="-|>", color="0.55", lw=0.6, ls=(0, (2, 1.5)), mutation_scale=5, shrinkA=1, shrinkB=1))
    ax.text(0.095, 0.375, "no feedback (~0.1%)", fontsize=4.2, color="0.45", ha="right", va="center", style="italic", rotation=90)
    # consequences
    tgs = box(0.74, 0.93, f"TGS: $f\\propto\\sqrt{{G/\\rho}}$  {100 * df:+.1f}%", SH)
    slip = box(0.74, 0.71, f"slip: $E_\\perp\\propto Gb^2$  {100 * dG:+.1f}%", SH)
    duc = box(0.74, 0.50, f"ductility: $\\varepsilon_F\\propto(K/G)^2$\n{100 * deF:+.1f}% (embrittlement)", BAD, fc="#FBEAE5")
    clv = box(0.74, 0.18, f"cleavage: $\\sigma_F\\propto Ka$  {100 * dK:+.1f}%", VO)
    arrow(g, tgs, SH); arrow(g, slip, SH); arrow(kg, duc, BAD); arrow(k, clv, VO)
    # the TGS inference: f measures the shear channel; with K fixed, d(K/G)/(K/G) = -2 df/f
    ax.annotate("", xy=(1.0, 0.5), xycoords=duc, xytext=(1.0, 0.5), textcoords=tgs,
                arrowprops=dict(arrowstyle="-|>", color=BAD, lw=0.8, ls=(0, (3, 2)), mutation_scale=6,
                                connectionstyle="arc,angleA=0,angleB=0,armA=7,armB=7,rad=3"), zorder=2)
    ax.text(0.985, 0.72, "TGS reads it: $\\Delta(K/G)/(K/G)\\approx-2\\,\\Delta f/f$", fontsize=4.3, color=BAD,
            ha="center", va="center", rotation=-90)


def main_v4(panel_h=None, name="fig5_embrittlement_v4", panel_n=None, row4=None):
    """v3 plus a fourth row: strain energy Psi(K, G) and how the chain comes together."""
    sys.path.insert(0, str(HERE / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    T, D = tensors(), md_decomp()
    W, H = 190 / 25.4, 9.4
    sc = 7.2 / 9.4                                   # rows 1-3 keep their v3 size in inches
    fig = plt.figure(figsize=(W, H))
    off = 1 - sc                                     # rows 1-3 occupy the top sc of the figure
    Y = lambda y: off + y * sc
    g1 = fig.add_gridspec(1, 3, left=0.055, right=0.745, top=Y(0.955), bottom=Y(0.715), wspace=0.12, width_ratios=[1.05, 1, 1])
    g3 = fig.add_gridspec(1, 4, left=0.085, right=0.995, top=Y(0.30), bottom=Y(0.06), wspace=0.42, width_ratios=[1, 1, 1.05, 0.95])
    a = [fig.add_subplot(g1[0, i]) for i in range(3)]
    atom_a(a[0]); atom_b(a[1]); atom_c(a[2])
    for ax, t in zip(a, ("(a) bonds: Hookean vs breaking", "(b) dislocation: core vs far field", "(c) crack tip: cleave or emit")):
        ax.set_title(t, fontsize=7, loc="left")
    dft_cells(fig, [0.755, Y(0.70), 0.24, 0.27 * sc])
    fig.text(0.757, Y(0.962), "(d) DFT cells: $\\alpha$-like vs $\\alpha'$-like", fontsize=7)
    md_cubes(fig, [0.045, Y(0.375), 0.29, 0.27 * sc])
    fig.text(0.05, Y(0.66), "(e) MD cells: Fe$-$36Cr before and after decomposition", fontsize=7)
    pin = fig.add_axes([0.37, Y(0.38), 0.17, 0.25 * sc]); atom_d(pin)
    fig.text(0.37, Y(0.66), "(f) $\\alpha'$ domains pin dislocations", fontsize=7)
    mp = fig.add_axes([0.62, Y(0.40), 0.35, 0.22 * sc]); md_phi(mp)
    fig.text(0.62, Y(0.66), "(g) MD: $G$ rises, $K$ flat, $K/G$ falls with $\\phi$", fontsize=7)
    mp.tick_params(labelsize=5.5); mp.xaxis.label.set_size(6); mp.yaxis.label.set_size(6)
    b = [fig.add_subplot(g3[0, i]) for i in range(4)]
    (panel_h or panel_b)(b[0], T, D); panel_c(b[1], T, D)
    panel_d(b[2], D, legend_at=(1.0, 0.37), legend_loc="right"); tgs_summary(b[3])
    for ax, t in zip(b, ("(h) bonding map", "(i) plastic range", "(j) TGS shift to fracture strain", "(k) what TGS sees")):
        ax.set_title(t, fontsize=7, loc="left")
        if t[1] != "k":
            ax.tick_params(labelsize=5.5); ax.xaxis.label.set_size(6); ax.yaxis.label.set_size(6)
    # row 4: strain energy
    rl, rm, rn = row4 or ([0.085, 0.04, 0.24, 0.155], [0.40, 0.04, 0.27, 0.155], [0.70, 0.0, 0.30, 0.215])
    l = fig.add_axes(rl); psi_split(l, D)
    m = fig.add_axes(rm); coherency(m, T)
    n = fig.add_axes(rn); (panel_n or psi_hub)(n, D)
    for ax, t in ((l, "(l) where $\\Psi$ stores energy (Fe$-$36Cr)"), (m, "(m) coherency strain energy"), (n, "(n) how it comes together")):
        ax.set_title(t, fontsize=7, loc="left")
    for ax in (l, m):
        ax.tick_params(labelsize=5.5); ax.xaxis.label.set_size(6); ax.yaxis.label.set_size(6)
    for y, lab in ((Y(0.84), "ATOMIC (DFT)"), (Y(0.52), "MOLECULAR (MD)"), (Y(0.18), "DATA & TGS"), (0.12, "STRAIN ENERGY")):
        fig.text(0.006, y, lab, fontsize=7.5, fontweight="heavy", color="0.1", va="center", ha="left", rotation=90)
    for y in (Y(0.695), Y(0.345), 0.238):
        fig.add_artist(matplotlib.lines.Line2D([0.01, 0.99], [y, y], color="0.75", lw=0.6, ls=(0, (4, 3))))
    out = ROOT / "manuscript" / "figures" / "kdense_scripta" / name
    fig.savefig(f"{out}.png", dpi=600); fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.svg")
    print("wrote", out)


def panel_h_v5(ax, T, D):
    """Pettifor-Pugh map read as composition trajectories: each series joined Fe -> Cr, the two criteria drawn as
    lines with the doubly brittle quadrant shaded, and the MD random -> decomposed paths labelled by at.% Cr."""
    gk = 1 / 1.75
    ax.fill_between([-20, 0], gk, 0.70, color="0.90", lw=0)
    ax.axvline(0, color="0.45", lw=0.6, ls="--"); ax.axhline(gk, color="0.45", lw=0.6, ls="--")
    ax.text(-19, 0.692, "brittle", fontsize=5.4, color="0.3", va="top", fontweight="bold")
    ax.text(69, gk - 0.006, "ductile", fontsize=5.4, color="0.3", ha="right", va="top", fontweight="bold")
    ax.text(-1.5, 0.345, "Pettifor: $C_{12}-C_{44}=0$", fontsize=4.4, color="0.4", ha="right", va="bottom", rotation=90)
    ax.text(69, gk + 0.003, "Pugh: $K/G=1.75$", fontsize=4.4, color="0.4", ha="right", va="bottom")
    # end labels (series that span the composition range); offsets in points
    ends = {"lit": {0: ("Fe", (8, 0)), -1: ("Cr", (0, -7))}, "md": {0: ("Fe", (-7, 0)), -1: ("Cr", (7, 2))},
            "qe": {0: ("20", (-7, 2)), -1: ("95", (2, 6))}}
    for k, (col, mk, lab) in SER.items():
        if k not in T:
            continue
        a = T[k]; cp, r = a[:, 3], a[:, 2] / a[:, 1]
        ax.plot(cp, r, "-", color=col, lw=0.6, alpha=0.6)
        ax.plot(cp, r, mk, color=col, mfc="white" if k == "qe" else col, mec=col, ms=3.2, ls="none", label=lab)
        for i, (txt, off) in ends.get(k, {}).items():
            ax.annotate(txt, (cp[i], r[i]), xytext=off, textcoords="offset points", fontsize=4.6, color=col,
                        ha="center", va="center", fontweight="bold")
    for x0, a in D.items():
        cp, r = a[:, 3], a[:, 2] / a[:, 1]
        ax.plot(cp, r, "-", color=DEC, lw=1.0)
        ax.annotate("", xy=(cp[-1], r[-1]), xytext=(cp[-2], r[-2]),
                    arrowprops=dict(arrowstyle="-|>", color=DEC, lw=1.0, mutation_scale=6))
        ax.plot([cp[-1]], [r[-1]], "s", mfc="white", mec=DEC, ms=3.4, mew=0.9)
        ax.annotate(f"{round(100 * x0)}", (cp[0], r[0]), xytext={0.2: (0, -6), 0.36: (-7, -4), 0.5: (3, -6)}[x0], textcoords="offset points", fontsize=4.4,
                    color=DEC, ha="center", va="center")
    ax.plot([], [], "s-", mfc="white", mec=DEC, color=DEC, ms=3.4, lw=1.0, label="MD random $\\to$ decomposed\n(20, 36, 50% Cr)")
    ax.set(xlabel=r"Cauchy pressure $C_{12}-C_{44}$ (GPa)", ylabel=r"$G/K$", xlim=(-20, 70), ylim=(0.34, 0.70))
    ax.legend(fontsize=4.4, loc="upper right", frameon=False, handlelength=1.4, borderaxespad=0.2, labelspacing=0.22)


def main_v5():
    """v4 with a readable bonding map (h): composition trajectories, criteria as lines, labelled decomposition paths."""
    main_v4(panel_h=panel_h_v5, name="fig5_embrittlement_v5")


def main_v6():
    """v5 with a mechanistic (n): decomposition -> shear vs volume channel of Psi -> consequences, TGS inference."""
    main_v4(panel_h=panel_h_v5, name="fig5_embrittlement_v6", panel_n=psi_flow,
            row4=([0.085, 0.04, 0.20, 0.155], [0.345, 0.04, 0.20, 0.155], [0.585, 0.0, 0.40, 0.212]))


if __name__ == "__main__":
    main()
    main_v1()
    main_v2()
    main_v3()
    main_v4()
