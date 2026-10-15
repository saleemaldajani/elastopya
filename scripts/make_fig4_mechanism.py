#!/usr/bin/env python3
"""Fig. 4 (draft): mechanism of the anomalous elasticity and of decomposition stiffening in bcc Fe-Cr.

FRUSTRATE -> NEST -> STIFFEN -> DETECT, left to right, every inset from data:
 (a) atoms & moments: Fe-rich / frustrated / decomposed motifs; QE <m_Cr>(x); QE and Levesque dG_demix(x)
 (b) electrons: Fermi-surface section of NM Cr (QE) with the SDW nesting vector Q; Fermi-level DOS of the SQS
 (c) electron -> lattice: f_SAW(x) (literature EMTO-CPA, MD, QE) with the lever-rule chord of Fe-36Cr decomposition
 (d) measurement: (P)SAW frequencies of Fe-36Cr random vs fully decomposed (MD, 10^4 grains) against the
     unallowed band of pristine CASS
Output: manuscript/figures/kdense_scripta/fig4_mechanism_v0.{png,pdf}
"""
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "dft_qe"))
import make_figures as mf  # noqa: E402
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap, to_rgb  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402

ROOT = HERE.parent
ENG = ROOT / "dft_qe" / "engaging"
FE, CR = "#F06060", "#6060F0"          # as Figs. 1-3 strips
QE_COL, LIT, MD = "#8F3A00", "0.12", "#009E73"
X0, XA, XB = 0.36, None, None           # Fe-36Cr and its MD alpha / alpha' compositions (phi = 1)


def glow(ax, x, y, color, lw=1.2, label=None, **kw):
    for w, a in ((7, 0.06), (4.5, 0.12), (2.6, 0.25)):
        ax.plot(x, y, color=color, lw=w * lw / 1.2, alpha=a, solid_capstyle="round",
                **{k: v for k, v in kw.items() if k not in ("marker", "ms")})
    return ax.plot(x, y, color=color, lw=lw, label=label, **kw)


def motif(ax, x0, y0, spins, title, s=0.13):
    """3 x 3 bcc (001)-projected motif: corner sites (i, j) and body centres (i + 1/2, j + 1/2)."""
    sites = [(i, j) for i in range(3) for j in range(3)] + [(i + 0.5, j + 0.5) for i in range(2) for j in range(2)]
    for k, (i, j) in enumerate(sites):
        el, m = spins[k]
        xc, yc = x0 + i * s, y0 + j * s
        ax.add_patch(plt.Circle((xc, yc), 0.04, color=FE if el == "Fe" else CR, ec="0.15", lw=0.4, zorder=3))
        if m:
            ax.annotate("", xy=(xc, yc + 0.06 * m), xytext=(xc, yc - 0.06 * m),
                        arrowprops=dict(arrowstyle="-|>", lw=0.8, color="black", mutation_scale=5), zorder=4)
        elif el == "Cr":
            ax.text(xc, yc, "?", ha="center", va="center", fontsize=6, color="white", fontweight="bold", zorder=4)
    ax.text(x0 + s, y0 - 0.07, title, ha="center", va="top", fontsize=6.0, linespacing=1.0)


def panel_a(fig, rect):
    ax = fig.add_axes([rect[0], rect[1] + rect[3] * 0.50, rect[2], rect[3] * 0.50])
    ax.set_xlim(-0.02, 1.0); ax.set_ylim(0.38, 1.05); ax.set_aspect("equal", adjustable="datalim"); ax.axis("off")
    fe_rich = [("Fe", 1)] * 13; fe_rich[4] = ("Cr", -1)
    frus = [("Fe", 1), ("Cr", 0), ("Fe", 1), ("Cr", 0), ("Fe", 1), ("Cr", 0), ("Fe", 1), ("Fe", 1), ("Cr", 0),
            ("Cr", 0), ("Fe", 1), ("Cr", 0), ("Fe", 1)]
    dec = [("Fe", 1), ("Fe", 1), ("Fe", 1), ("Fe", 1), ("Fe", 1), ("Fe", 1), ("Cr", 1), ("Cr", -1), ("Cr", 1),
           ("Fe", 1), ("Fe", 1), ("Cr", -1), ("Cr", -1)]
    motif(ax, -0.02, 0.62, fe_rich, "Fe-rich: Cr\nantiparallel")
    motif(ax, 0.36, 0.62, frus, "Cr-Cr\nfrustrated")
    ax.text(1.0, 1.04, "schematic", fontsize=4.6, color="0.5", ha="right", style="italic")
    motif(ax, 0.74, 0.62, dec, "decomposed\n$\\alpha$ FM | $\\alpha'$ AF")
    # inset: QE moments and demixing free energy
    ins = fig.add_axes([rect[0] + 0.055, rect[1] + 0.09, rect[2] * 0.58, rect[3] * 0.34])
    qm = {}
    import visualize_sqs as VS
    for x in (0.2, 0.35, 0.5, 0.65, 0.95):
        qm[x] = VS.mean_moments(x)
    xs = np.array(sorted(qm)); mcr = np.array([qm[x][1] for x in xs])
    glow(ins, xs, mcr, CR, lw=1.1, marker="o", ms=3)
    ins.axhline(0, color="0.6", lw=0.5)
    ins.set_ylabel(r"$\langle m_{\rm Cr}\rangle$ ($\mu_{\rm B}$)", fontsize=6, color=CR)
    ins.set_xlim(0, 1); ins.set_ylim(-0.7, 0.15); ins.tick_params(labelsize=5.5)
    ins.set_xlabel(r"$x_{\rm Cr}$", fontsize=6, labelpad=1)
    i2 = ins.twinx()
    q = mf.qe_mixing()
    if q is not None:
        i2.plot(q["x"], q["demix"], "s", mfc="white", mec=QE_COL, ms=3, mew=0.8)
    p = mf.read_csv("sheet_levesque_parameters.csv")
    xx = np.linspace(0.01, 0.99, 200)
    g, _ = mf.levesque_model(xx, p["alpha"][0], p["beta"][0], p["gamma"][0], p["delta"][0])
    i2.plot(xx, g, "-", color="0.35", lw=0.8)
    i2.axvspan(0.21, 0.75, color="#FFD27F", alpha=0.35, lw=0)
    i2.text(0.48, -0.005, "spinodal", ha="center", fontsize=5.5, color="#8A5A00")
    i2.set_ylabel(r"$\Delta G_{\rm demix}$ (eV)", fontsize=6); i2.tick_params(labelsize=5.5); i2.set_ylim(-0.09, 0.03)
    return ax


N_CDW = 32


def cdw_schematic(ax):
    """SDW (moments ~ cos Q.z, AF between neighbouring layers) and the CDW + strain wave at 2Q. delta exaggerated."""
    P, N = 16, N_CDW                                      # SDW period in layers (real Cr: ~40), layers shown
    n = np.arange(N); z = n.astype(float)
    env = np.cos(2 * np.pi * (n + 0.5) / P)
    zz = np.linspace(-0.5, N - 0.5, 600)
    for sgn in (1, -1):
        glow(ax, zz, 2.35 + sgn * 0.55 * np.cos(2 * np.pi * (zz + 0.5) / P), "#FFB000", lw=0.7, ls="--")
    ax.scatter(z, np.full(N, 2.35), s=9, color=CR, edgecolor="0.15", linewidth=0.3, zorder=3)
    for k in n:
        m = 0.55 * (-1) ** k * env[k]
        if abs(m) > 0.06:
            ax.annotate("", xy=(z[k], 2.35 + m), xytext=(z[k], 2.35 - m),
                        arrowprops=dict(arrowstyle="-|>", lw=0.6, color="black", mutation_scale=3.5), zorder=4)
    rho = np.cos(4 * np.pi * (zz + 0.5) / P)
    ax.fill_between(zz, 1.05, 1.05 + 0.32 * rho, where=rho > 0, color="#7FDBFF", alpha=0.75, lw=0)
    ax.fill_between(zz, 1.05, 1.05 + 0.32 * rho, where=rho < 0, color="#0B1026", alpha=0.55, lw=0)
    glow(ax, zz, 1.05 + 0.32 * rho, "#2A9DF4", lw=0.8)
    u = 0.28 * np.sin(4 * np.pi * (n + 0.5) / P)
    for k in n:
        ax.plot([z[k], z[k]], [-0.45, -0.15], color="0.7", lw=0.4)
    ax.scatter(z + u, np.full(N, -0.30), s=7, color=CR, edgecolor="0.15", linewidth=0.3, zorder=3)
    for y, lab in ((2.35, r"SDW  $m\propto\cos\mathbf{Q}\!\cdot\!\mathbf{r}$"), (1.05, r"CDW  $\Delta\rho\propto\cos 2\mathbf{Q}\!\cdot\!\mathbf{r}$"),
                   (-0.30, r"strain wave at $2\mathbf{Q}$")):
        ax.text(-1.2, y + 0.48, lab, fontsize=5.6, va="bottom", ha="left", color="0.1")
    ax.annotate("", xy=(P - 0.5, 3.2), xytext=(-0.5, 3.2), arrowprops=dict(arrowstyle="<->", lw=0.6, color="0.3"))
    ax.text(P / 2 - 0.5, 3.25, r"$\lambda_{\rm SDW}=a/\delta$ ($\approx 20a$; $\delta$ exaggerated)", ha="center", va="bottom", fontsize=5.3, color="0.3")
    ax.set_xlim(-1.5, N - 0.3); ax.set_ylim(-0.6, 3.8); ax.axis("off")


def ett_strip(fig, x0, y0, h, note_below=False):
    """Schematic Lifshitz electronic topological transition: two Fermi-surface pockets touch at x_c and open a
    neck as alloying moves E_F (Korzhavyi 2009: Fe majority sheet, 5-10 % Cr); elastic signature: sqrt-type
    anomaly in B and C' (Olsson 2003)."""
    g = np.linspace(-2.0, 2.0, 241); X, Y = np.meshgrid(g, g)
    F = (X ** 2 - 1) ** 2 + 2.2 * Y ** 2
    fw, fh = fig.get_size_inches()
    w = h * fh / fw                                            # square axes in figure units
    for i, (c, lab) in enumerate(((0.55, "pockets"), (1.0, "$x_c$: touch"), (1.6, "neck"))):
        ax = fig.add_axes([x0 + i * (w + 0.006), y0, w, h])
        ax.add_patch(plt.Rectangle((-2, -2), 4, 4, fc="#0B1026", ec="none"))
        cs = ax.contour(X, Y, F, levels=[c], colors="#7FDBFF", linewidths=0.8)
        for seg in cs.allsegs[0]:
            ax.plot(seg[:, 0], seg[:, 1], color="#7FDBFF", lw=3, alpha=0.25)
        ax.set_xlim(-2, 2); ax.set_ylim(-2, 2); ax.set_aspect("equal"); ax.axis("off")
        ax.text(0, -2.25, lab, ha="center", va="top", fontsize=5.0, color="0.15")
    xe = x0 + 3 * (w + 0.006)
    fig.text(x0, y0 + h + 0.008, r"ETT (Lifshitz): adding Cr moves $E_{\rm F}$ $\rightarrow$", fontsize=5.4, color="0.1", va="bottom")
    if note_below:
        fig.text(x0, y0 - 0.03, "Fe$\\uparrow$ sheet, 5$-$10% Cr $\\Rightarrow$ \u221a-cusp in $B$, $C'$",
                 fontsize=5.4, color="#8A5A00", va="top")
        fig.text(xe - 0.006, y0 - 0.058, "schematic", fontsize=4.6, color="0.5", ha="right", va="top", style="italic")
        return
    fig.text(xe + 0.003, y0 + h * 0.5, "Fe$\\uparrow$ sheet,\n5$-$10% Cr:\n\u221a-type cusp\nin $B$, $C'$",
             fontsize=5.0, color="#8A5A00", va="center", linespacing=1.05)
    fig.text(xe + 0.003, y0 - 0.012, "schematic", fontsize=4.6, color="0.5", ha="left", va="top", style="italic")


def draw_fs(ax, label=True):
    """Fermi-surface section (k_z = 0) of NM bcc Cr (QE) with the SDW nesting vector."""
    import fermi_slices as fs
    k, ev, ef = fs.read_xml(ENG / "fermi" / "Cr_nm" / "fs_data.xml")
    n = int(round(np.sqrt(len(k)))); kx, ky = k[:, 0].reshape(n, n), k[:, 1].reshape(n, n)
    ax.add_patch(plt.Polygon([(1, 0), (0, 1), (-1, 0), (0, -1)], closed=True, fc="#0B1026", ec="0.6", lw=0.6, zorder=0))
    for b in range(ev.shape[2]):
        e = ev[0, :, b].reshape(n, n) - ef
        if e.min() < 0 < e.max():
            cs = ax.contour(kx, ky, e, levels=[0], colors="#7FDBFF", linewidths=0.9, zorder=3)
            for seg in cs.allsegs[0]:
                for w, a in ((5, 0.10), (3, 0.22)):
                    ax.plot(seg[:, 0], seg[:, 1], color="#7FDBFF", lw=w, alpha=a, zorder=2)
    ax.annotate("", xy=(0.97, -0.02), xytext=(0.02, -0.02),
                arrowprops=dict(arrowstyle="-|>", lw=1.6, color="#FFB000", mutation_scale=9), zorder=5)
    ax.text(0.5, -0.12, r"$\mathbf{Q}_{\rm SDW}\approx(1-\delta)\,2\pi/a$", color="#FFB000", ha="center", va="top", fontsize=6.5,
            path_effects=[pe.withStroke(linewidth=2, foreground="#0B1026")], zorder=6)
    for lab, (x, y) in {"Γ": (0, 0), "H": (1, 0)}.items():
        ax.text(x + 0.03, y + 0.05, lab, fontsize=6.5, color="white", zorder=6)
    if label:
        ax.text(0.0, 1.07, "Cr: jack ($\\Gamma$) nests onto octahedron (H)\n$\\Rightarrow$ spin- and charge-density waves",
                ha="center", va="bottom", fontsize=5.8, color="0.15", zorder=6)
    ax.set_xlim(-1.02, 1.02); ax.set_ylim(-1.02, 1.3 if label else 1.02); ax.set_aspect("equal"); ax.axis("off")


def draw_nef(ins):
    """Fermi-level DOS per atom of the relaxed SQS (QE, projwfc total DOS within +-0.1 eV of E_F)."""
    xs, nef = [], []
    for d in sorted((ENG / "dos").glob("x*")):
        if "_" in d.name or not (d / "nscf.out").exists():
            continue
        pre = next((p for p in ("pdosW", "pdosT") if (d / f"{p}.pdos_tot").exists()), None)
        if pre is None:
            continue
        efd = float(re.findall(r"the Fermi energy is\s+([-\d.]+) ev", (d / "nscf.out").read_text())[-1])
        a = np.loadtxt(d / f"{pre}.pdos_tot"); m = np.abs(a[:, 0] - efd) < 0.1
        xs.append(float(d.name[1:])); nef.append((a[m, 1].mean() + a[m, 2].mean()) / 40)
    glow(ins, xs, nef, "#2A9DF4", lw=1.0, marker="o", ms=3)
    ins.axvspan(0.05, 0.10, color="#FFB000", alpha=0.35, lw=0)
    ins.text(0.075, 0.97, "ETT", rotation=90, fontsize=5, va="top", ha="center", color="#8A5A00", transform=ins.get_xaxis_transform())
    ins.set(xlim=(0, 1)); ins.tick_params(labelsize=5.5)
    ins.set_xlabel(r"$x_{\rm Cr}$", fontsize=6, labelpad=1); ins.set_ylabel(r"$N(E_{\rm F})$ (eV$^{-1}$)", fontsize=6)


def panel_b(fig, rect, cdw=False):
    draw_fs(fig.add_axes([rect[0] + 0.01, rect[1] + rect[3] * (0.62 if cdw else 0.52), rect[2] * 0.95,
                          rect[3] * (0.31 if cdw else 0.40)]))
    if cdw:
        ett_strip(fig, rect[0] + 0.012, rect[1] + rect[3] * 0.465, rect[3] * 0.095)
        cax = fig.add_axes([rect[0] + 0.005, rect[1] + 0.015, rect[2] * 0.98, rect[3] * 0.40])
        cdw_schematic(cax)
        cax.text(N_CDW - 1, -0.2, "schematic", fontsize=4.6, color="0.5", ha="right", style="italic")
        return
    ett_strip(fig, rect[0] + 0.012, rect[1] + rect[3] * 0.355, rect[3] * 0.11)
    draw_nef(fig.add_axes([rect[0] + 0.05, rect[1] + 0.09, rect[2] * 0.75, rect[3] * 0.19]))


def panel_c(fig, rect, eq_at=None):
    global XA, XB
    ax = fig.add_axes([rect[0] + 0.045, rect[1] + 0.09, rect[2] * 0.78, rect[3] * 0.71])
    lit, mdr, qe = mf.literature(), mf.md_random(), mf.qe_dft()
    pc = mf.read_csv("md_phase_compositions.csv")
    m = (np.abs(pc["x_cr"] - X0) < 1e-6) & (pc["phi"] == 1.0)
    XA, XB = float(pc["x_alpha"][m][0]), float(pc["x_alpha_prime"][m][0])
    w = (X0 - XA) / (XB - XA)
    for s, col, lab in ((lit, LIT, "electronic (EMTO-CPA)"), (mdr, MD, "classical (MD/EAM)")):
        glow(ax, s["x"], s["f"], col, lw=1.1, label=lab)
        fa, fb, f0 = (float(np.interp(v, s["x"], s["f"])) for v in (XA, XB, X0))
        ax.plot([XA, XB], [fa, fb], ":", color=col, lw=0.9)
        fl = (1 - w) * fa + w * fb
        ax.annotate("", xy=(X0, fl), xytext=(X0, f0), arrowprops=dict(arrowstyle="-|>", color="#FF3B00", lw=1.3, mutation_scale=8))
        ax.text(X0 + 0.02, 0.5 * (fl + f0), f"+{100 * (fl / f0 - 1):.1f}%", color="#FF3B00", fontsize=6.5, fontweight="bold", va="center")
    if qe is not None:
        ax.plot(qe["x"], qe["f"], "o", mfc="white", mec=QE_COL, ms=3.5, mew=0.9, label="QE (this work)")
    for v, lab in ((XA, r"$\alpha$"), (XB, r"$\alpha'$")):
        ax.axvline(v, color="0.75", lw=0.5, ls="--"); ax.text(v, 0.443, lab, ha="center", fontsize=6.5)
    ax.plot([X0], [0.443], "*", color="#0072B2", ms=8)
    ax.axvspan(0.05, 0.10, color="#FFB000", alpha=0.30, lw=0, zorder=0)
    ax.text(0.075, 0.535, "ETT", rotation=90, ha="center", va="center", fontsize=5.5, color="#8A5A00")
    ax.set(xlim=(0, 1), ylim=(0.44, 0.62)); ax.tick_params(labelsize=5.5)
    ax.set_xlabel(r"$x_{\rm Cr}$", fontsize=6, labelpad=1); ax.set_ylabel(r"$f_{\rm SAW}$ (GHz, 6.4 $\mu$m)", fontsize=6)
    ax.legend(fontsize=5.3, loc="upper left", bbox_to_anchor=(0.12, 1.0), frameon=False)
    eq = r"$C = C_{\rm ion} + C_{\rm band}$,  $C_{\rm band}\propto-\!\sum_k(\partial\varepsilon_k/\partial\eta)^2\,\delta(\varepsilon_k-E_{\rm F})$"
    if eq_at is None:
        ax.text(0.5, 1.03, eq, transform=ax.transAxes, ha="center", va="bottom", fontsize=5.6, color="0.2")
    else:
        fig.text(*eq_at, eq, ha="right", va="top", fontsize=5.8, color="0.2",
                 bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="0.8", lw=0.5))


def saw_top(prefix, x):
    d = np.loadtxt(mf.DATA / "saw" / f"{prefix}_x{x:.2f}.csv", delimiter=",", skiprows=1)
    keep = d[:, 3] / d[:, 2] >= 2e-3
    return np.r_[d[:, 0], d[keep, 1]] / 6.4e-6 / 1e6


def saw_schematic(ax):
    """TGS: crossed pump beams write a grating of period lambda; the Rayleigh wave runs along the surface over the
    random (left) and decomposed (right) Fe-36Cr microstructure (real MD planes); f = v_R / lambda."""
    from make_figures_kdense import md_grid
    from visualize_decomposition import series
    ser = series()[X0]
    lam = 1.0
    cm = LinearSegmentedColormap.from_list("fc", [FE, CR])
    ax.imshow(md_grid(ser[0]), extent=(0, 1.0, -1.05, 0), cmap=cm, interpolation="nearest", aspect="auto", alpha=0.85)
    ax.imshow(md_grid(ser[-1]), extent=(1.0, 2.0, -1.05, 0), cmap=cm, interpolation="nearest", aspect="auto", alpha=0.85)
    ax.axvline(1.0, ymin=0.0, ymax=0.58, color="white", lw=1.2)
    xx = np.linspace(0, 2, 400)
    glow(ax, xx, 0.07 * np.sin(2 * np.pi * xx / lam), "#FF3B00", lw=1.3)
    for x0 in (0.25, 0.75, 1.25, 1.75):                       # retrograde elliptical particle orbits decaying with depth
        for dz, sc in ((0.0, 1.0), (-0.3, 0.55), (-0.6, 0.25)):
            ax.add_patch(matplotlib.patches.Ellipse((x0, dz), 0.10 * sc, 0.16 * sc, fill=False, ec="white", lw=0.7, zorder=4))
    for xb in (0.30, 1.70):                                   # pump beams
        ax.annotate("", xy=(1.0, 0.1), xytext=(xb, 0.80), arrowprops=dict(arrowstyle="-", lw=2.5, color="#FFB000", alpha=0.6))
    for xg in np.arange(0, 2.01, lam / 2):
        ax.plot([xg, xg], [0.12, 0.42], color="#FFB000", lw=1.0, alpha=0.5)
    ax.annotate("", xy=(1.5, 0.48), xytext=(0.5, 0.48), arrowprops=dict(arrowstyle="<->", lw=0.6, color="0.2"))
    ax.text(1.0, 0.5, r"$\lambda$ = 6.4 $\mu$m grating", ha="center", va="bottom", fontsize=5.6)
    ax.text(0.5, -1.12, "random", ha="center", va="top", fontsize=5.8); ax.text(1.5, -1.12, "decomposed", ha="center", va="top", fontsize=5.8)
    ax.text(-0.02, 0.97, "schematic (MD planes real)", fontsize=4.6, color="0.5", ha="left", va="top", style="italic")
    ax.set_xlim(-0.02, 2.02); ax.set_ylim(-1.3, 1.0); ax.axis("off")


def grain_top(prefix, x):
    d = np.loadtxt(mf.DATA / "saw" / f"{prefix}_x{x:.2f}.csv", delimiter=",", skiprows=1)
    keep = d[:, 3] / d[:, 2] >= 2e-3
    return np.where(keep, np.maximum(d[:, 0], d[:, 1]), d[:, 0]) / 6.4e-6 / 1e6


def panel_d(fig, rect, saw=False):
    if saw:
        sax = fig.add_axes([rect[0] + 0.01, rect[1] + rect[3] * 0.50, rect[2] * 0.96, rect[3] * 0.46])
        saw_schematic(sax)
        dec = mf.md_decomposition()[X0]
        f0, f1 = float(dec["f"][np.argmin(dec["phi"])]) * 1e3, float(dec["f"][np.argmax(dec["phi"])]) * 1e3
        fig.text(rect[0] + 0.012, rect[1] + rect[3] * 0.465,
                 r"$C_{ij}\!\uparrow\Rightarrow G\!\uparrow\Rightarrow v_R=\eta(\nu)\sqrt{G/\rho}\!\uparrow\Rightarrow f=v_R/\lambda\!\uparrow$" "\n"
                 f"MD: {f0:.0f} $\\to$ {f1:.0f} MHz (+{100 * (f1 / f0 - 1):.1f}%); electronic: +7.7%",
                 fontsize=5.4, va="top", color="0.1", linespacing=1.3)
    ax = fig.add_axes([rect[0] + 0.06, rect[1] + 0.09, rect[2] * 0.70, rect[3] * (0.27 if saw else 0.71)])
    ax.axhspan(609, 833, color="0.90", lw=0); ax.axhline(609, color="0.4", lw=0.8, ls="--")
    ax.text(0.5, 720, "unallowed in\npristine CASS", ha="center", va="center", fontsize=6, color="0.35")
    for pos, pref, xx, col, lab in ((0, "md", X0, MD, "MD\nrandom"), (1, "mddec", X0, "#0B5D3B", "MD\ndecomposed"),
                                    (2, "qe", 0.35, QE_COL, "DFT (QE)\nrandom")):
        f = saw_top(pref, xx)
        frac = 100 * np.mean(grain_top(pref, xx) > 609)          # % of the 10^4 grains with a mode above 609 MHz
        ax.text(pos, f.max() + 12, f"{f.max():.0f}" + (f"\n{frac:.0f}% in band" if frac > 0 else ""), ha="center", va="bottom",
                fontsize=5.0, color="#FF3B00" if frac > 0 else "0.3", fontweight="bold" if frac > 0 else "normal")
        v = ax.violinplot(f, positions=[pos], widths=0.75, showextrema=False, bw_method=0.12)
        for b in v["bodies"]:
            b.set_facecolor(col); b.set_edgecolor(col); b.set_alpha(0.55)
        ax.plot([pos - 0.22, pos + 0.22], [f.max()] * 2, color="#FF3B00", lw=1.1)
    ax.set(xlim=(-0.6, 2.6), ylim=(330, 860)); ax.tick_params(labelsize=5.5)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(["MD\nrandom", "MD\ndecomposed", "DFT (QE)\nrandom"], fontsize=5.2, linespacing=1.0)
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel("(P)SAW (MHz)" if saw else "(P)SAW frequency (MHz), Fe$-$36Cr", fontsize=6)
    ax.text(0.99, 1.01, r"$10^4$ grains, $\lambda$ = 6.4 $\mu$m; QE at 35% Cr", transform=ax.transAxes, ha="right", va="bottom", fontsize=4.8, color="0.4")
    if saw:
        for t in ax.texts:
            if "unallowed" in t.get_text():
                t.set_text("unallowed"); t.set_fontsize(5.5)


def build_full():
    """All components, two rows: atoms & electrons (top), lattice & measurement (bottom)."""
    sys.path.insert(0, str(HERE / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    W, H = 190 / 25.4, 6.1
    fig = plt.figure(figsize=(W, H))
    bg = fig.add_axes([0, 0, 1, 1], zorder=-5); bg.axis("off")
    cmap = LinearSegmentedColormap.from_list("fecr", [(*to_rgb(FE), 1), (1, 1, 1, 1), (*to_rgb(CR), 1)])
    bg.imshow(np.linspace(0, 1, 512)[None, :], aspect="auto", cmap=cmap, alpha=0.10, extent=(0, 1, 0, 1))
    bg.axhline(0.505, xmin=0.01, xmax=0.99, color="0.75", lw=0.6, ls=(0, (4, 3)))
    top = [(0.005, "FRUSTRATE", "(a) atoms & moments"), (0.255, "NEST", "(b) Cr Fermi surface"),
           (0.505, "TOPOLOGY", "(c) Fe-rich: ETT"), (0.755, "MODULATE", "(d) spin & charge waves")]
    bot = [(0.005, "STIFFEN", "(e) electron $\\to$ lattice: decomposition"), (0.505, "PROPAGATE", "(f) Rayleigh wave in TGS"),
           (0.755, "DETECT", "(g) (P)SAW vs pristine CASS")]
    for row, y in ((top, 0.985), (bot, 0.49)):
        for i, (x, verb, sub) in enumerate(row):
            fig.text(x + 0.012, y, verb, fontsize=10.5, fontweight="heavy", color="0.08", va="top")
            fig.text(x + 0.012, y - 0.035, sub, fontsize=6.6, color="0.3", va="top")
            nxt = row[i + 1][0] if i + 1 < len(row) else None
            if nxt is not None:
                fig.patches.append(matplotlib.patches.FancyArrow(nxt - 0.045, y - 0.016, 0.03, 0, width=0.006, head_width=0.018,
                                                                 head_length=0.010, color="#FF3B00", transform=fig.transFigure, figure=fig))
    fig.patches.append(matplotlib.patches.FancyArrow(0.97, 0.53, 0, -0.035, width=0.004, head_width=0.012, head_length=0.012,
                                                     color="#FF3B00", transform=fig.transFigure, figure=fig))
    # top row
    panel_a(fig, [0.005, 0.505, 0.245, 0.42])
    draw_fs(fig.add_axes([0.262, 0.53, 0.23, 0.37]))
    ett_strip(fig, 0.515, 0.80, 0.075, note_below=True)
    draw_nef(fig.add_axes([0.545, 0.565, 0.18, 0.15]))
    cax = fig.add_axes([0.758, 0.535, 0.24, 0.36]); cdw_schematic(cax)
    cax.text(N_CDW - 1, -0.2, "schematic", fontsize=4.6, color="0.5", ha="right", style="italic")
    # bottom row
    panel_c(fig, [0.03, 0.0, 0.46, 0.425], eq_at=(0.49, 0.432))
    sax = fig.add_axes([0.51, 0.17, 0.235, 0.24]); saw_schematic(sax)
    dec = mf.md_decomposition()[X0]
    f0, f1 = float(dec["f"][np.argmin(dec["phi"])]) * 1e3, float(dec["f"][np.argmax(dec["phi"])]) * 1e3
    fig.text(0.515, 0.115, r"$C_{ij}\!\uparrow\Rightarrow G\!\uparrow\Rightarrow v_R=\eta(\nu)\sqrt{G/\rho}\!\uparrow\Rightarrow f=v_R/\lambda\!\uparrow$" "\n"
             f"Fe$-$36Cr, MD: {f0:.0f} $\\to$ {f1:.0f} MHz (+{100 * (f1 / f0 - 1):.1f}%)\nelectronic lever rule: +7.7%",
             fontsize=5.6, va="top", color="0.1", linespacing=1.4)
    panel_d(fig, [0.755, 0.0, 0.245, 0.43])
    out = ROOT / "manuscript" / "figures" / "kdense_scripta" / "fig4_mechanism_v1_full"
    fig.savefig(f"{out}.png", dpi=600); fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.svg")
    print("wrote", out)


N_EF_CR = 0.70     # states/eV/atom, NM bcc Cr (Moruzzi, Janak & Williams 1978, p. 76) -> rigid-band shift


def draw_fs_change(ax, y_fe=0.10, legend_below=False):
    """Cr Fermi-surface section at E_F (pure Cr) and at the rigid-band E_F with y_fe Fe (2 extra electrons per Fe):
    the electron jack (Gamma) and hole octahedron (H) change size in opposite senses, so Q no longer nests both."""
    import fermi_slices as fs
    k, ev, ef = fs.read_xml(ENG / "fermi" / "Cr_nm" / "fs_data.xml")
    n = int(round(np.sqrt(len(k)))); kx, ky = k[:, 0].reshape(n, n), k[:, 1].reshape(n, n)
    de = 2 * y_fe / N_EF_CR
    ax.add_patch(plt.Polygon([(1, 0), (0, 1), (-1, 0), (0, -1)], closed=True, fc="#0B1026", ec="0.6", lw=0.6, zorder=0))
    for b in range(ev.shape[2]):
        e = ev[0, :, b].reshape(n, n) - ef
        if e.min() < 0 < e.max() or e.min() < de < e.max():
            for lev, col, ls, lw in ((0.0, "#7FDBFF", "-", 0.9), (de, "#FF8C1A", "--", 0.9)):
                cs = ax.contour(kx, ky, e, levels=[lev], colors=col, linewidths=lw, linestyles=ls, zorder=3)
                for seg in cs.allsegs[0]:
                    ax.plot(seg[:, 0], seg[:, 1], color=col, lw=3, alpha=0.18, zorder=2)
    ax.annotate("", xy=(0.97, -0.02), xytext=(0.02, -0.02),
                arrowprops=dict(arrowstyle="-|>", lw=1.6, color="#FFB000", mutation_scale=9), zorder=5)
    ax.text(0.5, -0.12, r"$\mathbf{Q}_{\rm SDW}$", color="#FFB000", ha="center", va="top", fontsize=6.5,
            path_effects=[pe.withStroke(linewidth=2, foreground="#0B1026")], zorder=6)
    for lab, (x, y) in {"Γ": (0, 0), "H": (1, 0)}.items():
        ax.text(x + 0.03, y + 0.05, lab, fontsize=6.5, color="white", zorder=6)
    ax.plot([], [], "-", color="#7FDBFF", label="pure Cr: jack and octahedron nest")
    ax.plot([], [], "--", color="#FF8C1A", label=f"+{round(100 * y_fe)}% Fe (rigid band, $\\Delta E_F$ = +{de:.2f} eV):\nsheets resize oppositely $\\Rightarrow$ nesting lost")
    if legend_below:
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 0.0), fontsize=5.2, frameon=False, handlelength=1.8)
    else:
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.22), fontsize=5.2, frameon=False, handlelength=1.8)
    ax.set_xlim(-1.02, 1.02); ax.set_ylim(-1.02, 1.02); ax.set_aspect("equal"); ax.axis("off")


def draw_kohn(ax):
    """Gamma-H phonons of Cr: DFPT (NM, k = 24^3) and neutron data (Shaw & Muhlestein 1971, 300 K). The longitudinal
    branch softens towards H (Kohn anomaly at the nesting vector), the lattice imprint of the SDW/CDW."""
    gp = np.loadtxt(ENG / "phonon" / "Cr_nm_v5_k24" / "freq.gp")
    m = gp[:, 0] <= 1.0 + 1e-9
    z, w = gp[m, 0], gp[m, 1:] * 0.0299792458
    ax.plot(z, w[:, 2], "-", color="#2A9DF4", lw=1.1, label="QE DFPT (NM), L")
    ax.plot(z, w[:, 0], "-", color="0.55", lw=0.8, label="QE DFPT, T")
    import csv
    rows = [r for r in csv.reader(l for l in open(ROOT / "dft_qe" / "refdata" / "cr_phonon_shaw1971.csv") if not l.startswith("#"))][1:]
    for br, mk, lab in (("L", "o", "neutron (Shaw 1971), L"), ("T", "s", "neutron, T")):
        pts = np.array([[float(r[2]), float(r[3])] for r in rows if r[0] == "[00z]" and r[1] == br])
        if len(pts):
            ax.plot(pts[:, 0], pts[:, 1], mk, mfc="white", mec="0.15", ms=2.6, mew=0.6, label=lab)
    ax.axvspan(0.93, 0.97, color="#FFB000", alpha=0.35, lw=0)
    ax.text(0.95, 0.5, r"$\mathbf{Q}$", ha="center", fontsize=6, color="#8A5A00")
    ax.set(xlim=(0, 1), ylim=(0, 11.5)); ax.tick_params(labelsize=5.5)
    ax.set_xticks([0, 0.5, 1]); ax.set_xticklabels(["Γ", "", "H"])
    ax.set_ylabel(r"$\nu$ (THz)", fontsize=6, labelpad=1)
    ax.legend(fontsize=4.6, loc="lower right", frameon=False, handlelength=1.2)
    ax.text(-0.2, 1.03, "Kohn softening at $\\mathbf{Q}$:\nlattice imprint of the SDW/CDW", transform=ax.transAxes, fontsize=5.2,
            va="bottom", color="0.15", linespacing=1.1)


def build_full_v2():
    """As build_full, with the Fermi-surface change (b) and the measured/computed Kohn softening (d)."""
    sys.path.insert(0, str(HERE / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    W, H = 190 / 25.4, 6.1
    fig = plt.figure(figsize=(W, H))
    bg = fig.add_axes([0, 0, 1, 1], zorder=-5); bg.axis("off")
    cmap = LinearSegmentedColormap.from_list("fecr", [(*to_rgb(FE), 1), (1, 1, 1, 1), (*to_rgb(CR), 1)])
    bg.imshow(np.linspace(0, 1, 512)[None, :], aspect="auto", cmap=cmap, alpha=0.10, extent=(0, 1, 0, 1))
    bg.axhline(0.505, xmin=0.01, xmax=0.99, color="0.75", lw=0.6, ls=(0, (4, 3)))
    top = [(0.005, "FRUSTRATE", "(a) atoms & moments"), (0.255, "NEST", "(b) Cr Fermi surface vs Fe content"),
           (0.505, "TOPOLOGY", "(c) Fe-rich: ETT"), (0.755, "MODULATE", "(d) SDW/CDW soften the lattice")]
    bot = [(0.005, "STIFFEN", "(e) electron $\\to$ lattice: decomposition"), (0.505, "PROPAGATE", "(f) Rayleigh wave in TGS"),
           (0.755, "DETECT", "(g) (P)SAW vs pristine CASS")]
    for row, y in ((top, 0.985), (bot, 0.49)):
        for i, (x, verb, sub) in enumerate(row):
            fig.text(x + 0.012, y, verb, fontsize=10.5, fontweight="heavy", color="0.08", va="top")
            fig.text(x + 0.012, y - 0.035, sub, fontsize=6.6, color="0.3", va="top")
            nxt = row[i + 1][0] if i + 1 < len(row) else None
            if nxt is not None:
                fig.patches.append(matplotlib.patches.FancyArrow(nxt - 0.045, y - 0.016, 0.03, 0, width=0.006, head_width=0.018,
                                                                 head_length=0.010, color="#FF3B00", transform=fig.transFigure, figure=fig))
    fig.patches.append(matplotlib.patches.FancyArrow(0.97, 0.53, 0, -0.035, width=0.004, head_width=0.012, head_length=0.012,
                                                     color="#FF3B00", transform=fig.transFigure, figure=fig))
    panel_a(fig, [0.005, 0.505, 0.245, 0.42])
    draw_fs_change(fig.add_axes([0.268, 0.535, 0.22, 0.30]))
    ett_strip(fig, 0.515, 0.80, 0.075, note_below=True)
    draw_nef(fig.add_axes([0.545, 0.565, 0.18, 0.15]))
    cax = fig.add_axes([0.758, 0.70, 0.24, 0.215]); cdw_schematic(cax)
    cax.text(N_CDW - 1, -0.2, "schematic", fontsize=4.6, color="0.5", ha="right", style="italic")
    draw_kohn(fig.add_axes([0.80, 0.555, 0.185, 0.115]))
    panel_c(fig, [0.03, 0.0, 0.46, 0.425], eq_at=(0.49, 0.432))
    sax = fig.add_axes([0.51, 0.17, 0.235, 0.24]); saw_schematic(sax)
    dec = mf.md_decomposition()[X0]
    f0, f1 = float(dec["f"][np.argmin(dec["phi"])]) * 1e3, float(dec["f"][np.argmax(dec["phi"])]) * 1e3
    fig.text(0.515, 0.115, r"$C_{ij}\!\uparrow\Rightarrow G\!\uparrow\Rightarrow v_R=\eta(\nu)\sqrt{G/\rho}\!\uparrow\Rightarrow f=v_R/\lambda\!\uparrow$" "\n"
             f"Fe$-$36Cr, MD: {f0:.0f} $\\to$ {f1:.0f} MHz (+{100 * (f1 / f0 - 1):.1f}%)\nelectronic lever rule: +7.7%",
             fontsize=5.6, va="top", color="0.1", linespacing=1.4)
    panel_d(fig, [0.755, 0.0, 0.245, 0.43])
    out = ROOT / "manuscript" / "figures" / "kdense_scripta" / "fig4_mechanism_v2_full"
    fig.savefig(f"{out}.png", dpi=600); fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.svg")
    print("wrote", out)


def draw_fs_spin(ax_up, ax_dn):
    """Fe (FM, QE) Fermi-surface sections by spin. Majority: a band level at H lies 0.15 eV above E_F, where the
    H pockets connect (the critical point of the ETT Korzhavyi et al. find at 5-10 % Cr); minority: no critical
    point within +-0.3 eV. Rigid-band filling moves E_F down, away from it, so the ETT is driven by the Cr-induced
    shift of the majority states (magnetism), not by band filling."""
    import fermi_slices as fs
    k, ev, ef = fs.read_xml(ENG / "fermi" / "Fe_fm" / "fs_data.xml")
    n = int(round(np.sqrt(len(k)))); kx, ky = k[:, 0].reshape(n, n), k[:, 1].reshape(n, n)
    crit = None
    for s, ax, col in ((0, ax_up, "#FF8C1A"), (1, ax_dn, "#7FDBFF")):
        ax.add_patch(plt.Rectangle((-1, -1), 2, 2, fc="#0B1026", ec="none", zorder=0))
        for b in range(ev.shape[2]):
            e = ev[s, :, b].reshape(n, n) - ef
            if e.min() < 0 < e.max():
                cs = ax.contour(kx, ky, e, levels=[0], colors=col, linewidths=0.8, zorder=3)
                for seg in cs.allsegs[0]:
                    ax.plot(seg[:, 0], seg[:, 1], color=col, lw=2.6, alpha=0.2, zorder=2)
        ax.plot([1, 0, -1, 0, 1], [0, 1, 0, -1, 0], color="0.55", lw=0.4, zorder=1)
        ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_aspect("equal"); ax.axis("off")
    # majority critical level at H (degenerate band level just above E_F)
    iH = np.argmin((k[:, 0] - 1) ** 2 + k[:, 1] ** 2)
    lv = ev[0, iH, :] - ef
    crit = float(lv[(lv > 0) & (lv < 0.5)].min())
    for b in range(ev.shape[2]):
        e = ev[0, :, b].reshape(n, n) - ef
        if e.min() < crit - 0.005 < e.max():
            ax_up.contour(kx, ky, e, levels=[crit - 0.005], colors="white", linewidths=0.7, linestyles=":", zorder=4)
    for (x, y) in ((1, 0), (-1, 0), (0, 1), (0, -1)):
        ax_up.add_patch(plt.Circle((x, y), 0.13, fill=False, ec="#FFD400", lw=0.9, zorder=5))
    ax_up.text(0, 1.04, "spin $\\uparrow$ (majority)", ha="center", va="bottom", fontsize=5.6, color="#C25E00")
    ax_dn.text(0, 1.04, "spin $\\downarrow$ (minority)", ha="center", va="bottom", fontsize=5.6, color="#1B6FB5")
    ax_up.text(0, -1.05, f"H level {crit:+.2f} eV: pockets\nconnect (critical point)", ha="center", va="top", fontsize=4.8, color="0.15")
    ax_dn.text(0, -1.05, "no critical point\nwithin $\\pm$0.3 eV", ha="center", va="top", fontsize=4.8, color="0.15")
    return crit


def build_full_v3():
    """v2 with panel (c) showing the spin-resolved Fe Fermi surface behind the ETT (majority vs minority)."""
    sys.path.insert(0, str(HERE / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    W, H = 190 / 25.4, 6.1
    fig = plt.figure(figsize=(W, H))
    bg = fig.add_axes([0, 0, 1, 1], zorder=-5); bg.axis("off")
    cmap = LinearSegmentedColormap.from_list("fecr", [(*to_rgb(FE), 1), (1, 1, 1, 1), (*to_rgb(CR), 1)])
    bg.imshow(np.linspace(0, 1, 512)[None, :], aspect="auto", cmap=cmap, alpha=0.10, extent=(0, 1, 0, 1))
    bg.axhline(0.505, xmin=0.01, xmax=0.99, color="0.75", lw=0.6, ls=(0, (4, 3)))
    top = [(0.005, "FRUSTRATE", "(a) atoms & moments"), (0.255, "NEST", "(b) Cr Fermi surface vs Fe content"),
           (0.505, "TOPOLOGY", "(c) Fe-rich ETT: spin $\\uparrow$ vs $\\downarrow$"), (0.755, "MODULATE", "(d) SDW/CDW soften the lattice")]
    bot = [(0.005, "STIFFEN", "(e) electron $\\to$ lattice: decomposition"), (0.505, "PROPAGATE", "(f) Rayleigh wave in TGS"),
           (0.755, "DETECT", "(g) (P)SAW vs pristine CASS")]
    for row, y in ((top, 0.985), (bot, 0.49)):
        for i, (x, verb, sub) in enumerate(row):
            fig.text(x + 0.012, y, verb, fontsize=10.5, fontweight="heavy", color="0.08", va="top")
            fig.text(x + 0.012, y - 0.035, sub, fontsize=6.6, color="0.3", va="top")
            nxt = row[i + 1][0] if i + 1 < len(row) else None
            if nxt is not None:
                fig.patches.append(matplotlib.patches.FancyArrow(nxt - 0.045, y - 0.016, 0.03, 0, width=0.006, head_width=0.018,
                                                                 head_length=0.010, color="#FF3B00", transform=fig.transFigure, figure=fig))
    fig.patches.append(matplotlib.patches.FancyArrow(0.97, 0.53, 0, -0.035, width=0.004, head_width=0.012, head_length=0.012,
                                                     color="#FF3B00", transform=fig.transFigure, figure=fig))
    panel_a(fig, [0.005, 0.505, 0.245, 0.42])
    draw_fs_change(fig.add_axes([0.268, 0.535, 0.22, 0.30]))
    ett_strip(fig, 0.515, 0.80, 0.06, note_below=False)
    wsq = 0.105; hsq = wsq * W / H
    crit = draw_fs_spin(fig.add_axes([0.515, 0.585, wsq, hsq]), fig.add_axes([0.632, 0.585, wsq, hsq]))
    fig.text(0.515, 0.548, "Cr shifts the majority states (band filling alone moves\n$E_{\\rm F}$ away) $\\Rightarrow$ ETT at 5$-$10% Cr in spin $\\uparrow$ only\n"
             "(Korzhavyi 2009); unfolded SQS Fermi surfaces pending",
             fontsize=4.7, color="#8A5A00", va="top", linespacing=1.15)
    cax = fig.add_axes([0.758, 0.70, 0.24, 0.215]); cdw_schematic(cax)
    cax.text(N_CDW - 1, -0.2, "schematic", fontsize=4.6, color="0.5", ha="right", style="italic")
    draw_kohn(fig.add_axes([0.80, 0.555, 0.185, 0.115]))
    panel_c(fig, [0.03, 0.0, 0.46, 0.425], eq_at=(0.49, 0.432))
    sax = fig.add_axes([0.51, 0.17, 0.235, 0.24]); saw_schematic(sax)
    dec = mf.md_decomposition()[X0]
    f0, f1 = float(dec["f"][np.argmin(dec["phi"])]) * 1e3, float(dec["f"][np.argmax(dec["phi"])]) * 1e3
    fig.text(0.515, 0.115, r"$C_{ij}\!\uparrow\Rightarrow G\!\uparrow\Rightarrow v_R=\eta(\nu)\sqrt{G/\rho}\!\uparrow\Rightarrow f=v_R/\lambda\!\uparrow$" "\n"
             f"Fe$-$36Cr, MD: {f0:.0f} $\\to$ {f1:.0f} MHz (+{100 * (f1 / f0 - 1):.1f}%)\nelectronic lever rule: +7.7%",
             fontsize=5.6, va="top", color="0.1", linespacing=1.4)
    panel_d(fig, [0.755, 0.0, 0.245, 0.43])
    out = ROOT / "manuscript" / "figures" / "kdense_scripta" / "fig4_mechanism_v3_full"
    fig.savefig(f"{out}.png", dpi=600); fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.svg")
    print("wrote", out, "majority critical level at H:", round(crit, 3), "eV")


def build_full_v4(out_name="fig4_mechanism_v4_full"):
    """Everything in one figure: v3 plus the Fermi-level DOS back in panel (c); taller top row, no overlaps."""
    sys.path.insert(0, str(HERE / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    W, H = 190 / 25.4, 7.0
    fig = plt.figure(figsize=(W, H))
    bg = fig.add_axes([0, 0, 1, 1], zorder=-5); bg.axis("off")
    cmap = LinearSegmentedColormap.from_list("fecr", [(*to_rgb(FE), 1), (1, 1, 1, 1), (*to_rgb(CR), 1)])
    bg.imshow(np.linspace(0, 1, 512)[None, :], aspect="auto", cmap=cmap, alpha=0.10, extent=(0, 1, 0, 1))
    bg.axhline(0.425, xmin=0.01, xmax=0.99, color="0.75", lw=0.6, ls=(0, (4, 3)))
    top = [(0.005, "FRUSTRATE", "(a) atoms & moments"), (0.255, "NEST", "(b) Cr Fermi surface vs Fe content"),
           (0.505, "TOPOLOGY", "(c) Fe-rich ETT: spin $\\uparrow$ vs $\\downarrow$"), (0.755, "MODULATE", "(d) SDW/CDW soften the lattice")]
    if _MOD_V6:
        top[-1] = (0.755, "MODULATE", "(d) SDW softens $K$, not $G$")
    bot = [(0.005, "STIFFEN", "(e) electron $\\to$ lattice: decomposition"), (0.505, "PROPAGATE", "(f) Rayleigh wave in TGS"),
           (0.755, "DETECT", "(g) (P)SAW vs pristine CASS")]
    for row, y in ((top, 0.99), (bot, 0.41)):
        for i, (x, verb, sub) in enumerate(row):
            fig.text(x + 0.012, y, verb, fontsize=10.5, fontweight="heavy", color="0.08", va="top")
            fig.text(x + 0.012, y - 0.031, sub, fontsize=6.6, color="0.3", va="top")
            nxt = row[i + 1][0] if i + 1 < len(row) else None
            if nxt is not None:
                fig.patches.append(matplotlib.patches.FancyArrow(nxt - 0.045, y - 0.014, 0.03, 0, width=0.005, head_width=0.016,
                                                                 head_length=0.010, color="#FF3B00", transform=fig.transFigure, figure=fig))
    fig.patches.append(matplotlib.patches.FancyArrow(0.97, 0.445, 0, -0.03, width=0.004, head_width=0.012, head_length=0.010,
                                                     color="#FF3B00", transform=fig.transFigure, figure=fig))
    # (a)
    panel_a(fig, [0.005, 0.445, 0.245, 0.49])
    # (b)
    if _NEST_V5:
        draw_bands_gh(fig.add_axes([0.30, 0.812, 0.185, 0.09]))
        draw_fs_change(fig.add_axes([0.268, 0.50, 0.22, 0.30]), legend_below=True)
    else:
        draw_fs_change(fig.add_axes([0.268, 0.47, 0.22, 0.34]))
    # (c): ETT sequence, spin-resolved Fe Fermi surfaces, note, Fermi-level DOS
    ett_strip(fig, 0.515, 0.845, 0.048, note_below=False)
    wsq = 0.105; hsq = wsq * W / H
    crit = draw_fs_spin(fig.add_axes([0.515, 0.665, wsq, hsq]), fig.add_axes([0.632, 0.665, wsq, hsq]))
    fig.text(0.515, 0.622, "Cr shifts the majority states (band filling alone moves\n$E_{\\rm F}$ away) $\\Rightarrow$ ETT at 5$-$10% Cr in spin $\\uparrow$ only\n"
             "(Korzhavyi 2009); unfolded SQS Fermi surfaces pending", fontsize=4.7, color="#8A5A00", va="top", linespacing=1.15)
    draw_nef(fig.add_axes([0.552, 0.465, 0.175, 0.088]))
    # (d)
    if _MOD_V6:
        cax = fig.add_axes([0.758, 0.735, 0.24, 0.17]); cdw_schematic(cax)
        cax.text(N_CDW - 1, -0.62, "schematic", fontsize=4.6, color="0.5", ha="right", va="top", style="italic")
        draw_kohn_v6(fig.add_axes([0.795, 0.585, 0.19, 0.095]))
        draw_af_vs_nm(fig.add_axes([0.83, 0.462, 0.155, 0.072]))
    else:
        cax = fig.add_axes([0.758, 0.665, 0.24, 0.235]); cdw_schematic(cax)
        cax.text(N_CDW - 1, -0.62, "schematic", fontsize=4.6, color="0.5", ha="right", va="top", style="italic")
        draw_kohn(fig.add_axes([0.80, 0.465, 0.185, 0.145]))
    # bottom row
    panel_c(fig, [0.03, 0.0, 0.46, 0.335], eq_at=(0.49, 0.355))
    sax = fig.add_axes([0.51, 0.148, 0.235, 0.205]); saw_schematic(sax)
    dec = mf.md_decomposition()[X0]
    f0, f1 = float(dec["f"][np.argmin(dec["phi"])]) * 1e3, float(dec["f"][np.argmax(dec["phi"])]) * 1e3
    fig.text(0.515, 0.10, r"$C_{ij}\!\uparrow\Rightarrow G\!\uparrow\Rightarrow v_R=\eta(\nu)\sqrt{G/\rho}\!\uparrow\Rightarrow f=v_R/\lambda\!\uparrow$" "\n"
             f"Fe$-$36Cr, MD: {f0:.0f} $\\to$ {f1:.0f} MHz (+{100 * (f1 / f0 - 1):.1f}%)\nelectronic lever rule: +7.7%",
             fontsize=5.6, va="top", color="0.1", linespacing=1.4)
    panel_d(fig, [0.755, 0.0, 0.245, 0.37])
    out = ROOT / "manuscript" / "figures" / "kdense_scripta" / out_name
    fig.savefig(f"{out}.png", dpi=600); fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.svg")
    print("wrote", out)


def draw_bands_gh(ax, y_fe=0.10):
    """NM Cr bands along Gamma-H near E_F (QE): the d band that bounds the electron jack (Gamma side) and the hole
    octahedron (H side) crosses E_F at k_F; with y_fe Fe (rigid band) the crossing moves, resizing both sheets."""
    import fermi_slices as fs
    k, ev, ef = fs.read_xml(ENG / "bandpath" / "Cr_nm" / "bands_data.xml")
    d = np.r_[0, np.cumsum(np.linalg.norm(np.diff(k, axis=0), axis=1))]
    m = d <= 1 + 1e-6
    kk, E = d[m], ev[0][m] - ef
    de = 2 * y_fe / N_EF_CR
    for b in range(E.shape[1]):
        if E[:, b].min() < 1.2 and E[:, b].max() > -1.5:
            ax.plot(kk, E[:, b], color="0.2", lw=0.8)
    ax.axhline(0, color="#7FDBFF", lw=1.0); ax.axhline(de, color="#FF8C1A", lw=1.0, ls="--")
    kf = {}
    for lev, col in ((0.0, "#2A9DF4"), (de, "#FF8C1A")):
        for b in range(E.shape[1]):
            e = E[:, b] - lev; s = np.where(np.sign(e[:-1]) != np.sign(e[1:]))[0]
            for i in s:
                x = kk[i] - e[i] * (kk[i + 1] - kk[i]) / (e[i + 1] - e[i])
                if 0.5 < x < 0.75:
                    ax.plot([x], [lev], "o", color=col, ms=3, zorder=5); kf[lev] = x
    ax.set(xlim=(0, 1), ylim=(-1.2, 1.0)); ax.tick_params(labelsize=5)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Γ", "H"]); ax.set_ylabel("$E-E_F$ (eV)", fontsize=5.5, labelpad=1)
    if 0.0 in kf and de in kf:
        ax.text(0.0, 1.03, f"Γ$-$H bands: $k_F$ {kf[0.0]:.2f} $\\to$ {kf[de]:.2f} $2\\pi/a$ with +{round(100 * y_fe)}% Fe", transform=ax.transAxes,
                fontsize=4.9, va="bottom", color="0.15")
    ax.text(0.13, 0.43, "jack (e)", transform=ax.transAxes, fontsize=4.8, color="#1B6FB5", ha="center")
    ax.text(0.97, 0.05, "octahedron (h)", transform=ax.transAxes, fontsize=4.8, color="#1B6FB5", ha="right")


def build_full_v5():
    """v4 with the NEST column filled: Gamma-H bands near E_F on top, the Fermi-surface section below."""
    import types
    src = build_full_v4.__code__
    # reuse v4 by monkey-patching the two NEST calls through a flag
    global _NEST_V5
    _NEST_V5 = True
    try:
        build_full_v4(out_name="fig4_mechanism_v5_full")
    finally:
        _NEST_V5 = False


_NEST_V5 = False


_MOD_V6 = False


def draw_kohn_v6(ax):
    """As draw_kohn, but labelled for what it is: a dip of the zone-boundary L phonons near Q. Elastic constants are the
    slopes at Gamma, which the dip does not touch; NM DFPT is steeper than the (AF) neutron data there."""
    draw_kohn(ax)
    for t in list(ax.texts):
        if "Kohn softening" in t.get_text():
            t.remove()
    ax.legend_.remove()
    ax.text(-0.2, 1.04, "Kohn dip near $\\mathbf{Q}$: zone-boundary L phonons only;\n"
            "$C_{ij}$ = slopes at $\\Gamma$ (NM too steep vs AF neutron data)", transform=ax.transAxes, fontsize=4.6,
            va="bottom", color="0.15", linespacing=1.1)
    ax.legend(fontsize=4.0, loc="lower right", frameon=False, handlelength=1.0, labelspacing=0.15)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Γ", "H"])


def draw_af_vs_nm(ax):
    """Pure Cr (QE, 2 atoms, k = 16^3): antiferromagnetic (commensurate SDW) vs non-magnetic. Magnetic order cuts the
    bulk modulus (C12 collapses) and leaves the shear moduli unchanged."""
    import json
    t = {}
    for tag in ("nm", "afm"):
        r = json.loads((ENG / "runs" / f"Cr_{tag}_k16x16x16_e60" / "result.json").read_text())
        c11, c12, c44 = r["C11"], r["C12"], r["C44"]
        t[tag] = {"$K$": (c11 + 2 * c12) / 3, "$C'$": (c11 - c12) / 2, "$C_{44}$": c44}
    labs = list(t["nm"]); y = np.arange(len(labs))[::-1]
    ax.barh(y + 0.19, [t["nm"][l] for l in labs], height=0.36, color="0.75", label="NM")
    ax.barh(y - 0.19, [t["afm"][l] for l in labs], height=0.36, color="#2A9DF4", label="AF (SDW)")
    for yi, l in zip(y, labs):
        d = t["afm"][l] / t["nm"][l] - 1
        ax.text(max(t["nm"][l], t["afm"][l]) + 6, yi, f"{100 * d:+.0f}%", va="center", fontsize=4.8,
                color="#C02F00" if abs(d) > 0.1 else "0.3", fontweight="bold" if abs(d) > 0.1 else "normal")
    ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=5.2)
    ax.set_xlim(0, 330); ax.tick_params(axis="x", labelsize=4.8); ax.set_xlabel("Cr, QE (GPa)", fontsize=5.2, labelpad=1)
    ax.legend(fontsize=4.2, loc="lower right", frameon=False, handlelength=0.9, borderaxespad=0.1)
    dk = 100 * (t["afm"]["$K$"] / t["nm"]["$K$"] - 1)
    ax.text(-0.42, 1.08, f"AF order: $K$ {dk:+.0f}%, shear unchanged".replace("-", "$-$"), transform=ax.transAxes, fontsize=4.8, color="0.15", va="bottom")


def build_full_v6():
    """v5 with an honest (d): the SDW softens the volume channel (K, via C12) and the zone-boundary L phonons (Kohn dip),
    not the shear moduli that decomposition stiffens."""
    global _NEST_V5, _MOD_V6
    _NEST_V5 = _MOD_V6 = True
    try:
        build_full_v4(out_name="fig4_mechanism_v6_full")
    finally:
        _NEST_V5 = _MOD_V6 = False


# ---------------------------------------------------------------------------------------------------------------------
# v7 (2026-10-15): panels redrawn from results that did not exist when v6 was made. Earlier versions are unchanged.
#  (a) <m_Cr>(x) from every relaxed QE SQS (5-95 % Cr) instead of five compositions
#  (c) the schematic Lifshitz strip is replaced by our unfolded SQS Fermi surfaces (spin up, exp. a, 0/5/10 % Cr)
#  (d) the SDW/CDW schematic is replaced by the computed Cr supercell (moments, site charges, forces) and the chain
#      m^2 -> magnetic pressure -> K -> f_SAW is given in numbers
#  (e) the exp.-a QE series is added, with its (absent) lever-rule gain stated
#  (g) a fourth violin: QE Fe-35Cr at the experimental lattice parameter
# ---------------------------------------------------------------------------------------------------------------------
_V7 = False
QE_LEVER = (float('nan'), float('nan'))
QE_AEXP = "#0072B2"


def qe_series_f(suf):
    """(x, f_SAW GHz) of the QE SQS + pure elements for suf '' (PBE-relaxed) or '_aexp' (experimental a)."""
    import json
    R = ENG / "runs"; out = []
    dirs = sorted(R.glob(f"sqs40_x*_k442_e60{suf}"))
    dirs = [d for d in dirs if suf or not d.name.split("e60")[1]]
    dirs += sorted(R.glob("sqs40_x0.80_k442_e60_aexp_crafm" if suf else "sqs40_x0.80_k442_e60_crafm2"))
    pure = {"": ("Fe_fm_k16x16x16_e60", "Cr_afm_k16x16x16_e60"), "_aexp": ("Fe_fm_k16x16x16_e60_aexp", "Cr_afm_k16x16x16_e60_aexp")}[suf]
    items = [(d / "result.json", None) for d in dirs] + [(R / n / "result.json", x) for x, n in zip((0.0, 1.0), pure)]
    for j, xfix in items:
        if j.exists():
            r = json.loads(j.read_text()); x = r["x_cr"] if xfix is None else xfix
            c11, c12, c44 = r["C11"], r["C12"], r["C44"]
            b, cp = (c11 + 2 * c12) / 3, (c11 - c12) / 2
            g = ((2 * cp + 3 * c44) / 5 + 5 * (2 * cp) * c44 / (4 * c44 + 3 * (2 * cp))) / 2
            e = 9 * b * g / (3 * b + g)
            out.append((x, float(mf.rayleigh_f(b, e, x * mf.RHO_CR + (1 - x) * mf.RHO_FE))))
    return np.array(sorted(set(out)))


def lever_gain(xs, fs, xa, xb, x0=X0):
    fa, fb, f0 = (float(np.interp(v, xs, fs)) for v in (xa, xb, x0))
    w = (x0 - xa) / (xb - xa)
    return 100 * (((1 - w) * fa + w * fb) / f0 - 1)


def panel_a_v7(fig, rect):
    ax = panel_a(fig, rect)
    import visualize_sqs as VS
    ins = [a for a in fig.axes if a.get_ylabel().startswith(r"$\langle m_{\rm Cr}")][-1]
    for ln in list(ins.lines):
        if ln.get_color() not in ("0.6",):
            ln.remove()
    for c in list(ins.collections):
        c.remove()
    pts = []
    for x in (0.05, 0.10, 0.20, 0.30, 0.35, 0.40, 0.45, 0.475, 0.50, 0.65, 0.70, 0.95):
        try:
            m = VS.mean_moments(x)[1]
        except Exception:
            m = None
        if m is not None:
            pts.append((x, m))
    xs, mcr = np.array(pts).T
    glow(ins, xs, mcr, CR, lw=1.1, marker="o", ms=2.6)
    ins.set_ylim(min(-0.9, mcr.min() - 0.05), 0.15)
    return ax


def draw_unfolded_row(fig, x0, y0, w):
    """Spin-up unfolded Fermi surfaces of the SQS at 0/5/10 % Cr (exp. a): the real version of the Lifshitz cartoon."""
    import unfold_fs as uf
    W, H = fig.get_size_inches()
    gap = 0.008; s = (w - 2 * gap) / 3; hs = s * W / H
    for i, x in enumerate((0.00, 0.05, 0.10)):
        ax = fig.add_axes([x0 + i * (s + gap), y0, s, hs])
        _, A, _ = uf.spectral(ENG / "unfold" / f"x{x:.2f}_aexp", "up")
        ax.imshow(A, extent=(-1, 1, -1, 1), origin="lower", cmap="magma", vmin=0, vmax=np.percentile(A, 99.3))
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.4)
        ax.text(0, -1.08, f"{round(100 * x)}% Cr", ha="center", va="top", fontsize=5.2)
        for (px, py, t) in ((-1, 1, "H"), (0, 1, "N")):
            ax.text(px * 0.86, py * 0.86, t, ha="center", va="center", fontsize=4.2, color="white")
    fig.text(x0, y0 + hs + 0.006, "our unfolded SQS, spin $\\uparrow$ at $E_{\\rm F}$ (exp. $a$): arms through N fade by 5% Cr",
             fontsize=4.9, color="0.15", va="bottom")


def sdw_numbers():
    """Numbers of the chain m^2 -> magnetic pressure -> K -> f_SAW from the Cr supercell and the 2-atom Cr tensors."""
    import json
    import sdw_cr
    res = {k: sdw_cr.sites(ENG / "sdw" / k / "scf.out") for k in ("nm", "afm", "sdw")}
    P = {}
    for k in ("nm", "afm", "sdw"):
        t = (ENG / "sdw" / k / "scf.out").read_text()
        P[k] = float(re.findall(r"P=\s*([-\d.]+)", t)[-1]) / 10           # kbar -> GPa
    m2 = {k: float(np.mean(res[k][1] ** 2)) for k in ("afm", "sdw")}
    tens = {}
    for tag in ("nm", "afm"):
        r = json.loads((ENG / "runs" / f"Cr_{tag}_k16x16x16_e60" / "result.json").read_text())
        c11, c12, c44 = r["C11"], r["C12"], r["C44"]
        b, cp = (c11 + 2 * c12) / 3, (c11 - c12) / 2
        g = ((2 * cp + 3 * c44) / 5 + 5 * (2 * cp) * c44 / (4 * c44 + 3 * (2 * cp))) / 2
        e = 9 * b * g / (3 * b + g)
        tens[tag] = (b, g, float(mf.rayleigh_f(b, e, mf.RHO_CR)) * 1e3)
    # strain wave: first harmonic of the forces at 2Q, static continuum response u = F / (C11 Omega q^2)
    a, NL = sdw_cr.A, sdw_cr.NL
    fz = res["sdw"][2] * 13.605693 / 0.529177                             # Ry/bohr -> eV/A
    F = np.abs(np.fft.rfft(fz) / NL * 2); k = int(np.argmax(F[1:]) + 1)
    q = 2 * np.pi * k / (NL * a / 2)
    rj = json.loads((ENG / "runs" / "Cr_afm_k16x16x16_e60" / "result.json").read_text())
    u = F[k] / (rj["C11"] / 160.21766 * a ** 3 / 2 * q ** 2)
    return dict(Pmag_afm=P["afm"] - P["nm"], Pmag_sdw=P["sdw"] - P["nm"], m2_afm=m2["afm"], m2_sdw=m2["sdw"],
                K_nm=tens["nm"][0], K_afm=tens["afm"][0], G_nm=tens["nm"][1], G_afm=tens["afm"][1],
                f_nm=tens["nm"][2], f_afm=tens["afm"][2], u_A=u, strain=q * u, period_A=NL * a / 2 / k)


def draw_sdw_computed(ax):
    """Computed Cr SDW supercell (QE): moments m(z) as stems, site-charge CDW and the force (strain-wave drive) at 2Q."""
    import sdw_cr
    ch, mg, fz, _ = sdw_cr.sites(ENG / "sdw" / "sdw" / "scf.out")
    z = np.arange(len(mg)) * sdw_cr.A / 2
    for zi, mi in zip(z, mg):
        ax.plot([zi, zi], [2.5, 2.5 + 0.45 * mi], color="0.15", lw=0.6)
    ax.scatter(z, 2.5 + 0.45 * mg, s=3, color=CR, zorder=3)
    glow(ax, z, 2.5 + 0.45 * np.abs(mg), "#FFB000", lw=0.6, ls="--")
    dq = (ch - ch.mean()) / np.abs(ch - ch.mean()).max()
    ax.fill_between(z, 1.3, 1.3 + 0.28 * dq, where=dq > 0, color="#7FDBFF", alpha=0.75, lw=0)
    ax.fill_between(z, 1.3, 1.3 + 0.28 * dq, where=dq < 0, color="#0B1026", alpha=0.55, lw=0)
    glow(ax, z, 1.3 + 0.28 * dq, "#2A9DF4", lw=0.7)
    fn = fz / np.abs(fz).max()
    glow(ax, z, 0.25 + 0.28 * fn, "#009E73", lw=0.7)
    ax.axhline(0.25, color="0.75", lw=0.3)
    for y, lab in ((3.1, "SDW moments $m(z)$ (computed, QE)"), (1.63, "CDW: site charge"), (0.58, "force at $2\\mathbf{Q}$ (drives the strain wave)")):
        ax.text(-1.0, y, lab, fontsize=5.0, va="bottom", ha="left", color="0.1")
    ax.set_xlim(-1.2, z[-1] + 1.0); ax.set_ylim(-0.42, 3.45); ax.axis("off")
    ax.text(z[-1] + 0.8, -0.40, r"Cr, $1\times1\times20$ cell, $\delta=1/20$", fontsize=4.4, color="0.45", ha="right", va="bottom")


def panel_c_v7(fig, rect, eq_at=None):
    panel_c(fig, rect, eq_at=eq_at)
    ax = fig.axes[-1]
    s = qe_series_f("_aexp")
    if len(s):
        ax.plot(s[:, 0], s[:, 1], "D", mfc="none", mec=QE_AEXP, ms=3.0, mew=0.8, label="QE, exp. $a$")
        ax.legend(fontsize=5.0, loc="upper left", bbox_to_anchor=(0.12, 1.0), frameon=False)
        pc = mf.read_csv("md_phase_compositions.csv")
        r = qe_series_f("")
        g_rel = lever_gain(r[:, 0], r[:, 1], XA, XB) if len(r) else float("nan")
        g_aex = lever_gain(s[:, 0], s[:, 1], XA, XB)
        global QE_LEVER
        QE_LEVER = (g_rel, g_aex)


def panel_d_v7(fig, rect):
    ax = fig.add_axes([rect[0] + 0.06, rect[1] + 0.09, rect[2] * 0.72, rect[3] * 0.71])
    ax.axhspan(609, 833, color="0.90", lw=0); ax.axhline(609, color="0.4", lw=0.8, ls="--")
    ax.text(1.5, 720, "unallowed in\npristine CASS", ha="center", va="center", fontsize=6, color="0.35")
    cases = ((0, "md", X0, MD), (1, "mddec", X0, "#0B5D3B"), (2, "qe", 0.35, QE_COL), (3, "qeaexp", 0.35, QE_AEXP))
    for pos, pref, xx, col in cases:
        p = mf.DATA / "saw" / f"{pref}_x{xx:.2f}.csv"
        if not p.exists():
            continue
        f = saw_top(pref, xx)
        frac = 100 * np.mean(grain_top(pref, xx) > 609)
        ax.text(pos, f.max() + 12, f"{f.max():.0f}" + (f"\n{frac:.0f}% in band" if frac >= 0.5 else ""), ha="center", va="bottom",
                fontsize=5.0, color="#FF3B00" if frac >= 0.5 else "0.3", fontweight="bold" if frac >= 0.5 else "normal")
        v = ax.violinplot(f, positions=[pos], widths=0.75, showextrema=False, bw_method=0.12)
        for b in v["bodies"]:
            b.set_facecolor(col); b.set_edgecolor(col); b.set_alpha(0.55)
        ax.plot([pos - 0.22, pos + 0.22], [f.max()] * 2, color="#FF3B00", lw=1.1)
    ax.set(xlim=(-0.6, 3.6), ylim=(330, 860)); ax.tick_params(labelsize=5.5)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(["MD\nrandom", "MD\ndecomp.", "QE\nrelaxed", "QE\nexp. $a$"], fontsize=5.0, linespacing=1.0)
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel("(P)SAW frequency (MHz), Fe$-$36Cr", fontsize=6)
    ax.text(0.99, 1.01, r"$10^4$ grains, $\lambda$ = 6.4 $\mu$m; QE at 35% Cr", transform=ax.transAxes, ha="right", va="bottom", fontsize=4.8, color="0.4")


def build_full_v7(out_name="fig4_mechanism_v7_full", pc=None, pd=None):
    """v6 with the panels the later results made possible (see the block comment above); earlier versions unchanged."""
    sys.path.insert(0, str(HERE / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    W, H = 190 / 25.4, 7.0
    fig = plt.figure(figsize=(W, H))
    bg = fig.add_axes([0, 0, 1, 1], zorder=-5); bg.axis("off")
    cmap = LinearSegmentedColormap.from_list("fecr", [(*to_rgb(FE), 1), (1, 1, 1, 1), (*to_rgb(CR), 1)])
    bg.imshow(np.linspace(0, 1, 512)[None, :], aspect="auto", cmap=cmap, alpha=0.10, extent=(0, 1, 0, 1))
    bg.axhline(0.425, xmin=0.01, xmax=0.99, color="0.75", lw=0.6, ls=(0, (4, 3)))
    top = [(0.005, "FRUSTRATE", "(a) atoms & moments"), (0.255, "NEST", "(b) Cr Fermi surface vs Fe content"),
           (0.505, "TOPOLOGY", "(c) Fe-rich ETT: spin $\\uparrow$ vs $\\downarrow$"), (0.755, "MODULATE", "(d) SDW couples to volume, not shear")]
    bot = [(0.005, "STIFFEN", "(e) electron $\\to$ lattice: decomposition"), (0.505, "PROPAGATE", "(f) Rayleigh wave in TGS"),
           (0.755, "DETECT", "(g) (P)SAW vs pristine CASS")]
    for row, y in ((top, 0.99), (bot, 0.41)):
        for i, (x, verb, sub) in enumerate(row):
            fig.text(x + 0.012, y, verb, fontsize=10.5, fontweight="heavy", color="0.08", va="top")
            fig.text(x + 0.012, y - 0.031, sub, fontsize=6.6, color="0.3", va="top")
            nxt = row[i + 1][0] if i + 1 < len(row) else None
            if nxt is not None:
                fig.patches.append(matplotlib.patches.FancyArrow(nxt - 0.045, y - 0.014, 0.03, 0, width=0.005, head_width=0.016,
                                                                 head_length=0.010, color="#FF3B00", transform=fig.transFigure, figure=fig))
    fig.patches.append(matplotlib.patches.FancyArrow(0.97, 0.445, 0, -0.03, width=0.004, head_width=0.012, head_length=0.010,
                                                     color="#FF3B00", transform=fig.transFigure, figure=fig))
    panel_a_v7(fig, [0.005, 0.445, 0.245, 0.49])
    draw_bands_gh(fig.add_axes([0.30, 0.812, 0.185, 0.09]))
    draw_fs_change(fig.add_axes([0.268, 0.50, 0.22, 0.30]), legend_below=True)
    # (c): real unfolded Fermi surfaces replace the Lifshitz cartoon
    draw_unfolded_row(fig, 0.515, 0.835, 0.225)
    wsq = 0.105; hsq = wsq * W / H
    draw_fs_spin(fig.add_axes([0.515, 0.640, wsq, hsq]), fig.add_axes([0.632, 0.640, wsq, hsq]))
    fig.text(0.515, 0.597, "Cr brings the majority H level to $E_{\\rm F}$\n(0.10 $\\to$ 0.005 eV at 5% Cr, exp. $a$) $\\Rightarrow$ ETT at\n"
             "$x_c\\approx0.08$, spin $\\uparrow$ only (S6; Korzhavyi: 5$-$10%)",
             fontsize=4.7, color="#8A5A00", va="top", linespacing=1.15)
    draw_nef(fig.add_axes([0.552, 0.465, 0.175, 0.075]))
    # (d): computed SDW + the numbers that tie it to the SAW
    n = sdw_numbers()
    draw_sdw_computed(fig.add_axes([0.758, 0.735, 0.24, 0.17]))
    fig.text(0.765, 0.728, f"$P_{{\\rm mag}}$ vs NM: AF {n['Pmag_afm']:.1f}, SDW {n['Pmag_sdw']:.1f} GPa\n"
             f"SDW/AF {n['Pmag_sdw'] / n['Pmag_afm']:.2f} $\\approx\\langle m^2\\rangle$ ratio {n['m2_sdw'] / n['m2_afm']:.2f} (Landau)\n"
             f"strain wave $u\\approx${n['u_A']:.2f} $\\AA$, $\\varepsilon\\approx${100 * n['strain']:.1f}%, period {n['period_A']:.0f} $\\AA$:\n"
             f"averages out over a 6.4 $\\mu$m SAW", fontsize=4.5, color="0.15", va="top", linespacing=1.2)
    kax = fig.add_axes([0.795, 0.583, 0.125, 0.062]); draw_kohn_v6(kax)
    kax.legend(fontsize=4.0, loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, handlelength=1.0, labelspacing=0.2)   # beside, not on, the data
    for t in kax.get_legend().get_texts():
        t.set_text(t.get_text().replace(" (Shaw 1971)", ""))      # cited in the caption; keeps the legend inside the figure
    draw_af_vs_nm(fig.add_axes([0.83, 0.462, 0.155, 0.072]))
    for t in fig.texts:
        if t.get_text().startswith("AF order"):
            t.remove()
    fig.axes[-1].text(-0.42, 1.08, f"AF order: $K$ {100 * (n['K_afm'] / n['K_nm'] - 1):+.0f}%, $f_{{\\rm SAW}}$ "
                      f"{100 * (n['f_afm'] / n['f_nm'] - 1):+.1f}% (shear unchanged)".replace("-", "$-$"),
                      transform=fig.axes[-1].transAxes, fontsize=4.8, color="0.15", va="bottom")
    for t in list(fig.axes[-1].texts):
        if t.get_text().startswith("AF order: $K$") and "SAW" not in t.get_text():
            t.remove()
    # bottom row
    (pc or panel_c_v7)(fig, [0.03, 0.0, 0.46, 0.335], eq_at=(0.49, 0.355))
    sax = fig.add_axes([0.51, 0.148, 0.235, 0.205]); saw_schematic(sax)
    dec = mf.md_decomposition()[X0]
    f0, f1 = float(dec["f"][np.argmin(dec["phi"])]) * 1e3, float(dec["f"][np.argmax(dec["phi"])]) * 1e3
    fig.text(0.515, 0.10, r"$C_{ij}\!\uparrow\Rightarrow G\!\uparrow\Rightarrow v_R=\eta(\nu)\sqrt{G/\rho}\!\uparrow\Rightarrow f=v_R/\lambda\!\uparrow$" "\n"
             f"Fe$-$36Cr, MD: {f0:.0f} $\\to$ {f1:.0f} MHz (+{100 * (f1 / f0 - 1):.1f}%, 5 alloys)\nelectronic lever rule (EMTO-CPA): +7.7%\n"
             f"QE lever rule: {QE_LEVER[0]:+.1f}% (relaxed), {QE_LEVER[1]:+.1f}% (exp. $a$): no gain".replace("-", "$-$").replace("$$-$$", "$-$"),
             fontsize=5.6, va="top", color="0.1", linespacing=1.4)
    (pd or panel_d_v7)(fig, [0.755, 0.0, 0.245, 0.37])
    out = ROOT / "manuscript" / "figures" / "kdense_scripta" / out_name
    fig.savefig(f"{out}.png", dpi=600); fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.svg")
    print("wrote", out)
    print(json.dumps({k: round(v, 3) for k, v in n.items()}))


# v8 (2026-10-15, user review of v7): (e) the Fe-36Cr star and the alpha/alpha' labels sat on the x axis (ylim
# started at 0.44 GHz, the markers at 0.443): the axis now starts lower and they sit above it. (g) the violins showed
# every mode of every grain pooled, while the in-band percentage uses each grain's highest mode, so the in-band tail
# looked thinner than the percentage; v8 draws each grain's highest mode, the quantity compared with the band.
def panel_c_v8(fig, rect, eq_at=None):
    panel_c_v7(fig, rect, eq_at=eq_at)
    ax = fig.axes[-1]
    for ln in ax.lines:
        if ln.get_marker() == "*":
            ln.set_ydata([0.4335])
    for t in ax.texts:
        if t.get_position()[1] == 0.443:
            t.set_position((t.get_position()[0], 0.4290))
    ax.set_ylim(0.425, 0.62)


def panel_d_v8(fig, rect):
    ax = fig.add_axes([rect[0] + 0.06, rect[1] + 0.09, rect[2] * 0.72, rect[3] * 0.71])
    ax.axhspan(609, 833, color="0.90", lw=0); ax.axhline(609, color="0.4", lw=0.8, ls="--")
    ax.text(1.5, 760, "unallowed in\npristine CASS", ha="center", va="center", fontsize=6, color="0.35")
    for pos, pref, xx, col in ((0, "md", X0, MD), (1, "mddec", X0, "#0B5D3B"), (2, "qe", 0.35, QE_COL), (3, "qeaexp", 0.35, QE_AEXP)):
        if not (mf.DATA / "saw" / f"{pref}_x{xx:.2f}.csv").exists():
            continue
        g = grain_top(pref, xx)                                   # highest (P)SAW mode of each of the 10^4 grains
        frac = 100 * np.mean(g > 609)
        ax.text(pos, g.max() + 12, f"{g.max():.0f}" + (f"\n{frac:.0f}% in band" if frac >= 0.5 else ""), ha="center", va="bottom",
                fontsize=5.0, color="#FF3B00" if frac >= 0.5 else "0.3", fontweight="bold" if frac >= 0.5 else "normal")
        v = ax.violinplot(g, positions=[pos], widths=0.75, showextrema=False, bw_method=0.12)
        for b in v["bodies"]:
            b.set_facecolor(col); b.set_edgecolor(col); b.set_alpha(0.55)
        ax.plot([pos - 0.22, pos + 0.22], [g.max()] * 2, color="#FF3B00", lw=1.1)
    ax.set(xlim=(-0.6, 3.6), ylim=(340, 860)); ax.tick_params(labelsize=5.5)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(["MD\nrandom", "MD\ndecomp.", "QE\nrelaxed", "QE\nexp. $a$"], fontsize=5.0, linespacing=1.0)
    ax.tick_params(axis="x", length=0)
    ax.set_ylabel("highest (P)SAW per grain (MHz), Fe$-$36Cr", fontsize=6)
    ax.text(0.99, 1.01, r"$10^4$ grains, $\lambda$ = 6.4 $\mu$m; QE at 35% Cr", transform=ax.transAxes, ha="right", va="bottom", fontsize=4.8, color="0.4")


def build_full_v8():
    """v7 with the (e) markers lifted off the axis and (g) drawn as each grain's highest mode."""
    build_full_v7(out_name="fig4_mechanism_v8_full", pc=panel_c_v8, pd=panel_d_v8)


def main():
    build_full()
    build_full_v2()
    build_full_v3()
    build_full_v4()
    build_full_v5()
    for cdw, saw, tag in ((False, False, ""), (True, False, "_cdw"), (False, True, "_saw"), (True, True, "_cdw_saw")):
        build(cdw, saw, tag)


def build(cdw, saw, tag):
    sys.path.insert(0, str(HERE / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    W, H = 190 / 25.4, 3.35
    fig = plt.figure(figsize=(W, H))
    bg = fig.add_axes([0, 0, 1, 1], zorder=-5); bg.axis("off")
    cmap = LinearSegmentedColormap.from_list("fecr", [(*to_rgb(FE), 1), (1, 1, 1, 1), (*to_rgb(CR), 1)])
    bg.imshow(np.linspace(0, 1, 512)[None, :], aspect="auto", cmap=cmap, alpha=0.10, extent=(0, 1, 0, 1))
    cols = [(0.005, "FRUSTRATE", "atoms & moments"), (0.255, "NEST", "electrons at $E_{\\rm F}$"),
            (0.505, "STIFFEN", "electron $\\to$ lattice"), (0.755, "DETECT", "surface acoustic waves")]
    for i, (x, verb, sub) in enumerate(cols):
        fig.text(x + 0.012, 0.955, verb, fontsize=11, fontweight="heavy", color="0.08", va="top")
        fig.text(x + 0.012, 0.885, f"({'abcd'[i]}) {sub}", fontsize=6.8, color="0.3", va="top")
        if i < 3:
            fig.patches.append(matplotlib.patches.FancyArrow(x + 0.215, 0.925, 0.03, 0, width=0.012, head_width=0.035,
                                                             head_length=0.012, color="#FF3B00", transform=fig.transFigure, figure=fig))
    panel_a(fig, [0.005, 0.0, 0.245, 0.86])
    panel_b(fig, [0.255, 0.0, 0.245, 0.86], cdw=cdw)
    panel_c(fig, [0.505, 0.0, 0.245, 0.86])
    panel_d(fig, [0.755, 0.0, 0.245, 0.86], saw=saw)
    out = ROOT / "manuscript" / "figures" / "kdense_scripta" / f"fig4_mechanism_v0{tag}"
    fig.savefig(f"{out}.png", dpi=600); fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.svg")
    print("wrote", out, "alpha/alpha' =", XA, XB)


if __name__ == "__main__":
    main()
