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
    motif(ax, 0.74, 0.62, dec, "decomposed\n$\\alpha$ FM | $\\alpha'$ AF")
    # inset: QE moments and demixing free energy
    ins = fig.add_axes([rect[0] + 0.04, rect[1] + 0.09, rect[2] * 0.62, rect[3] * 0.34])
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


def panel_b(fig, rect):
    import fermi_slices as fs
    ax = fig.add_axes([rect[0] + 0.01, rect[1] + rect[3] * 0.40, rect[2] * 0.95, rect[3] * 0.58])
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
    ax.text(0.0, 1.07, "Cr: jack ($\\Gamma$) nests onto octahedron (H)\n$\\Rightarrow$ spin- and charge-density waves",
            ha="center", va="bottom", fontsize=5.8, color="0.15", zorder=6)
    ax.set_xlim(-1.02, 1.02); ax.set_ylim(-1.02, 1.3); ax.set_aspect("equal"); ax.axis("off")
    ins = fig.add_axes([rect[0] + 0.05, rect[1] + 0.09, rect[2] * 0.75, rect[3] * 0.24])
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
    glow(ins, xs, nef, "#7FDBFF", lw=1.0, marker="o", ms=3)
    ins.axvspan(0.05, 0.10, color="#FFB000", alpha=0.35, lw=0)
    ins.text(0.075, 0.97, "ETT", rotation=90, fontsize=5, va="top", ha="center", color="#8A5A00", transform=ins.get_xaxis_transform())
    ins.set(xlim=(0, 1)); ins.tick_params(labelsize=5.5)
    ins.set_xlabel(r"$x_{\rm Cr}$", fontsize=6, labelpad=1); ins.set_ylabel(r"$N(E_{\rm F})$ (eV$^{-1}$)", fontsize=6)


def panel_c(fig, rect):
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
    ax.set(xlim=(0, 1), ylim=(0.44, 0.62)); ax.tick_params(labelsize=5.5)
    ax.set_xlabel(r"$x_{\rm Cr}$", fontsize=6, labelpad=1); ax.set_ylabel(r"$f_{\rm SAW}$ (GHz, 6.4 $\mu$m)", fontsize=6)
    ax.legend(fontsize=5.3, loc="upper left", frameon=False)
    ax.text(0.5, 1.03, r"$C = C_{\rm ion} + C_{\rm band}$,  $C_{\rm band}\propto-\!\sum_k(\partial\varepsilon_k/\partial\eta)^2\,\delta(\varepsilon_k-E_{\rm F})$",
            transform=ax.transAxes, ha="center", va="bottom", fontsize=5.6, color="0.2")


def saw_top(prefix, x):
    d = np.loadtxt(mf.DATA / "saw" / f"{prefix}_x{x:.2f}.csv", delimiter=",", skiprows=1)
    keep = d[:, 3] / d[:, 2] >= 2e-3
    return np.r_[d[:, 0], d[keep, 1]] / 6.4e-6 / 1e6


def panel_d(fig, rect):
    ax = fig.add_axes([rect[0] + 0.06, rect[1] + 0.09, rect[2] * 0.70, rect[3] * 0.71])
    ax.axhspan(609, 833, color="0.90", lw=0); ax.axhline(609, color="0.4", lw=0.8, ls="--")
    ax.text(0.5, 720, "unallowed in\npristine CASS", ha="center", va="center", fontsize=6, color="0.35")
    for pos, pref, col, lab in ((0, "md", MD, "random"), (1, "mddec", "#0B5D3B", r"decomposed")):
        f = saw_top(pref, X0)
        v = ax.violinplot(f, positions=[pos], widths=0.75, showextrema=False, bw_method=0.12)
        for b in v["bodies"]:
            b.set_facecolor(col); b.set_edgecolor(col); b.set_alpha(0.55)
        ax.text(pos, 350, lab, ha="center", fontsize=6.3)
        ax.plot([pos - 0.2, pos + 0.2], [f.max()] * 2, color="#FF3B00", lw=1.2)
    ax.annotate("", xy=(1, 600), xytext=(0, 560), arrowprops=dict(arrowstyle="-|>", color="#FF3B00", lw=1.2, connectionstyle="arc3,rad=-0.3"))
    ax.set(xlim=(-0.6, 1.6), ylim=(330, 860), xticks=[]); ax.tick_params(labelsize=5.5)
    ax.set_ylabel("(P)SAW frequency (MHz), Fe$-$36Cr", fontsize=6)


def main():
    sys.path.insert(0, str(HERE / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta")
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
    panel_b(fig, [0.255, 0.0, 0.245, 0.86])
    panel_c(fig, [0.505, 0.0, 0.245, 0.86])
    panel_d(fig, [0.755, 0.0, 0.245, 0.86])
    out = ROOT / "manuscript" / "figures" / "kdense_scripta" / "fig4_mechanism_v0"
    fig.savefig(f"{out}.png", dpi=600); fig.savefig(f"{out}.pdf")
    print("wrote", out, "alpha/alpha' =", XA, XB)


if __name__ == "__main__":
    main()
