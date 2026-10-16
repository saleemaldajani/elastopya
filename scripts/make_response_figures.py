#!/usr/bin/env python3
"""Figures for the response to M. P. Short that have no counterpart in the manuscript.

Same K-Dense 'scripta' style and checks as make_figures_kdense.py, at the Elsevier 1.5-column
width (140 mm). Output: manuscript/figures/kdense_scripta/response/R*.{pdf,png,export.json,audit.json}
  R_convergence    QE convergence study for bcc Fe (k-mesh, cutoff), vs the present VASP tensor
  R_md_validation  25^3 small-cell procedure vs the Acta 50^3 series, and MD random alloys
  R_decomposition  stiffening on decomposition at all 12 compositions, and the band check
  R_phonon_convergence  NM bcc Cr DFPT phonons: imaginary modes at k = 16^3 vanish with denser k / wider smearing
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_figures_kdense as mk  # noqa: E402  (sets the Scripta palette on make_figures)
from make_figures_kdense import mf, plt  # noqa: E402
from style_presets import apply_publication_style  # noqa: E402

OUT = mf.ROOT / "manuscript" / "figures" / "kdense_scripta" / "response"
W_IN = mk.PROFILE["widths_mm"]["one-and-half"] / 25.4
VASP_FE = dict(C11=324.43, C12=169.84, C44=127.72)      # present VASP tensor of bcc Fe (sheet_dft.csv)
X_WINDOWS = ((0.21, 0.75, "PWSCF–PAW"), (0.25, 0.99, "SIESTA–NC"))   # spinodal windows, Sec. 3.1


def finish(fig, stem, w_pad=1.2):
    fig.set_size_inches(W_IN, fig.get_size_inches()[1] * W_IN / fig.get_size_inches()[0])
    mk.normalize_text(fig)
    fig.tight_layout(w_pad=w_pad, pad=0.3)
    mk.export_figure(fig, OUT / stem, formats=["pdf", "png", "svg"], dpi=mk.DPI, bbox_inches=None, overwrite=True,
                     mkdir=True, write_manifest=True, provenance={"source": "scripts/make_response_figures.py"})
    found, sizes = mk.audit(fig)
    (OUT / f"{stem}.audit.json").write_text(json.dumps(dict(issues=found, font_sizes_pt=sizes), indent=1) + "\n")
    plt.close(fig)
    return stem, found, min(sizes)


def convergence():
    d = mf.read_csv("qe_fe_convergence.csv")
    k = np.array([str(v) for v in np.genfromtxt(mf.DATA / "qe_fe_convergence.csv", delimiter=",", names=True,
                                                 dtype=None, encoding="utf-8")["kmesh"]])
    order = ["2x2x5", "8x8x8", "12x12x12", "16x16x16", "20x20x20"]
    labels = ["2×2×5\n(Γ of SQS)", "8³", "12³", "16³", "20³"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(5.5, 2.6), gridspec_kw={"width_ratios": [1.6, 1]})
    cols = {"C11": mf.BLUE, "C12": mf.RED, "C44": mf.GREEN}
    for name, col in cols.items():
        y = [d[f"{name}_GPa"][(k == m) & (d["ecutwfc_Ry"] == 60.0)][0] for m in order]
        a.plot(range(5), y, "-o", color=col, ms=4, lw=1.0, label=f"$C_{{{name[1:]}}}$, 60 Ry")
        y80 = d[f"{name}_GPa"][(k == "16x16x16") & (d["ecutwfc_Ry"] == 80.0)][0]
        a.plot([3], [y80], "s", mfc="white", mec=col, mew=1.0, ms=5, zorder=4)
        a.plot([4.35, 4.75], [VASP_FE[name]] * 2, color=col, lw=2.0, solid_capstyle="butt")
    a.plot([], [], "s", mfc="white", mec="0.3", mew=1.0, ms=5, ls="none", label="80 Ry (16³)")
    a.plot([], [], color="0.3", lw=2.0, label="present VASP (40-atom SQS)")
    a.set_xticks(range(5)); a.set_xticklabels(labels); a.minorticks_on(); a.xaxis.set_minor_locator(plt.NullLocator())
    a.set(xlim=(-0.3, 4.9), ylim=(60, 410), xlabel="$k$-mesh (2-atom bcc cell)", ylabel="Elastic constant (GPa)")
    a.legend(loc="upper left", ncol=2, columnspacing=1.0)
    mk.style(a, "(a)")
    m = [d["moment_muB_per_atom"][(k == mm) & (d["ecutwfc_Ry"] == 60.0)][0] for mm in order]
    b.plot(range(5), m, "-o", color="0.15", ms=4, lw=1.0)
    b.set_xticks(range(5)); b.set_xticklabels(labels); b.xaxis.set_minor_locator(plt.NullLocator())
    b.set(xlim=(-0.3, 4.3), ylim=(1.5, 2.4), xlabel="$k$-mesh", ylabel=r"Fe moment ($\mu_{\rm B}$/atom)")
    mk.style(b, "(b)")
    return finish(fig, "R_convergence")


def md_validation():
    acta = mf.read_csv("sheet_md.csv")
    dd = mf.read_csv("md_decomposition.csv")
    m = np.abs(dd["x_cr"] - 0.5) < 1e-9
    o = np.argsort(dd["phi"][m])
    ph, c11, c12 = dd["phi"][m][o], dd["c11"][m][o], dd["c12"][m][o]
    fig, (a, b) = plt.subplots(1, 2, figsize=(5.5, 2.5))
    a.plot(acta["phi"], acta["bulk_gpa"], "o", color="0.15", ms=3.5, zorder=4, label=r"Acta, 50$^3$ (10 alloys)")
    a.plot(ph, (c11 + 2 * c12) / 3, "s", mfc="white", mec=mf.GREEN, mew=1.2, ms=6, label=r"present, 25$^3$ (1 alloy)")
    a.set(xlabel=r"Spinodal conversion fraction, $\phi$", ylabel="Bulk modulus (GPa), Fe$_{50}$Cr$_{50}$",
          xlim=(0.47, 1.03), ylim=(193.5, 199.5))
    a.legend(loc="upper left")
    a.text(0.98, 0.05, f"$C_{{11}}(\\phi=1)$: {acta['c11_phi1_gpa'][0]:.1f} vs {c11[-1]:.1f} GPa",
           transform=a.transAxes, ha="right", va="bottom")
    mk.style(a, "(a)")
    r = mf.read_csv("md_random_alloys.csv")
    b.errorbar(r["x_cr"], r["c11"], yerr=r["c11_sd"], fmt="s", color=mf.GREEN, ms=4, capsize=2, label="present, random alloys")
    b.plot([0.5], [252.7], "o", mfc="none", mec="0.15", mew=1.2, ms=8, label=r"Acta, 50/50 random ($\phi$ = 0.5)")
    b.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"MD $C_{11}$ (GPa)", xlim=(-0.03, 1.03), ylim=(230, 340))
    b.legend(loc="upper left")
    mk.style(b, "(b)")
    return finish(fig, "R_md_validation")


def decomposition():
    dec = mf.md_decomposition()
    xs = np.array(sorted(dec))
    df = np.array([100 * (dec[x]["f"][-1] / dec[x]["f"][0] - 1) for x in xs])
    fig, (a, b) = plt.subplots(1, 2, figsize=(5.5, 2.6))
    for (lo, hi, lab), alpha in zip(X_WINDOWS, (0.18, 0.09)):
        a.axvspan(lo, hi, color=mf.GREEN, alpha=alpha, lw=0)
    a.text(0.48, 0.3, "spinodal window\n(Fig. 1b)", ha="center", va="bottom", color="0.3")
    a.plot(xs, df, "-s", color="0.15", ms=4, lw=1.0)
    a.plot([mf.X_FE36CR], [np.interp(mf.X_FE36CR, xs, df)], "*", color=mf.BLUE, ms=11, zorder=5)
    a.annotate("Fe–36Cr", (mf.X_FE36CR, np.interp(mf.X_FE36CR, xs, df)), xytext=(-7, -2),
               textcoords="offset points", ha="right", va="top", color=mf.BLUE)
    a.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"$\Delta f_{\rm SAW}$, random $\to$ $\phi$ = 1 (%)", xlim=(0, 1), ylim=(0, 6.6))
    mk.style(a, "(a)")

    def top(prefix, x):
        p = mf.DATA / "saw" / f"{prefix}_x{x:.2f}.csv"
        if not p.exists():
            return np.nan
        d = np.loadtxt(p, delimiter=",", skiprows=1)
        keep = d[:, 3] / d[:, 2] >= 2e-3
        return np.where(keep, np.maximum(d[:, 0], d[:, 1]), d[:, 0]).max() / 6.4e-6 / 1e6
    b.axhspan(609, 640, color="0.92", lw=0)
    b.axhline(609, color="0.35", lw=0.9, ls="--")
    b.text(0.02, 612, "unallowed range begins (609 MHz)", va="bottom", color="0.3")
    xr = [0.0, 0.2, 0.36, 0.5, 0.95, 1.0]
    b.plot(xr, [top("md", x) for x in xr], "s", color=mf.GREEN, ms=4.5, label="MD, random alloy")
    b.plot(xs, [top("mddec", x) for x in xs], "s", mfc="white", mec="#0B5D3B", mew=1.2, ms=5.5,
           label=r"MD, fully decomposed ($\phi$ = 1)")
    b.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"Highest (P)SAW, 10$^4$ grains (MHz)", xlim=(-0.03, 1.03), ylim=(566, 625))
    b.legend(loc="lower right")
    mk.style(b, "(b)")
    return finish(fig, "R_decomposition")


def phonon_convergence():
    sys.path.insert(0, str(mf.ROOT / "dft_qe"))
    import bands_path as bp
    ph = mf.ROOT / "dft_qe" / "engaging" / "phonon"
    thz = lambda v: np.loadtxt(ph / v / "freq.gp")
    fig, (a, b) = plt.subplots(1, 2, figsize=(5.5, 2.6), gridspec_kw={"width_ratios": [1.35, 1]})
    for v, col, ls, lab in (("Cr_nm", mf.RED, "--", "$k$ = 16$^3$, 0.02 Ry"), ("Cr_nm_v5_k24", "0.1", "-", "$k$ = 24$^3$, 0.02 Ry")):
        g = thz(v); x, w = g[:, 0], g[:, 1:] * 0.0299792458
        for j in range(w.shape[1]):
            a.plot(x, w[:, j], color=col, ls=ls, lw=0.9, label=lab if j == 0 else None)
    _, _, ticks = bp.path_points(bp.PATH_BCC)
    sc = x[-1] / ticks[-1][0]
    a.set_xticks([t * sc for t, _ in ticks]); a.set_xticklabels([l for _, l in ticks]); a.xaxis.set_minor_locator(plt.NullLocator())
    for t, _ in ticks[1:-1]:
        a.axvline(t * sc, color="0.85", lw=0.5)
    a.axhline(0, color="0.5", lw=0.6)
    a.set(xlim=(0, x[-1]), ylim=(-5.5, 14.5), ylabel="Frequency (THz)")
    a.legend(loc="upper center", ncol=2, columnspacing=0.8, handlelength=1.6)
    mk.style(a, "(a)")
    runs = [("Cr_nm", "16$^3$, default mixing"), ("Cr_nm_v2_mix", "16$^3$, gentler mixing"),
            ("Cr_nm_v6_a2843", "16$^3$, $a$ − 0.4%"), ("Cr_nm_v4_deg04", "16$^3$, smearing 0.04 Ry"),
            ("Cr_nm_v3_k20", "20$^3$"), ("Cr_nm_v5_k24", "24$^3$")][::-1]
    lo = [float((thz(v)[:, 1:] * 0.0299792458).min()) for v, _ in runs]
    cols = [mf.RED if m < -0.05 else ("0.1" if m > -0.01 else mf.ORANGE) for m in lo]
    b.barh(range(len(runs)), lo, color=cols, height=0.6)
    for i, m in enumerate(lo):
        if m < -0.05:
            b.text(m - 0.12, i, f"{m:.2f}", va="center", ha="right")
        else:
            b.plot([0], [i], "o", color="0.1", ms=4, zorder=4)
            b.text(-0.25, i, "0.00, all real", va="center", ha="right")
    b.axvline(0, color="0.5", lw=0.6)
    b.set_yticks(range(len(runs))); b.set_yticklabels([l for _, l in runs]); b.yaxis.set_minor_locator(plt.NullLocator())
    b.set(xlim=(-6.2, 1.2), xlabel="Lowest phonon frequency (THz)")
    mk.style(b, "(b)")
    return finish(fig, "R_phonon_convergence")


def main():
    apply_publication_style("scripta", palette_name="okabe_ito_on_white")
    plt.rcParams["svg.fonttype"] = "none"
    for fn in (convergence, md_validation, decomposition, phonon_convergence):
        stem, found, mn = fn()
        print(stem, "min pt", mn, "issues", found)


if __name__ == "__main__":
    main()


def md_validation_v2():
    """R_md_validation with the five MD alloys of the 25^3 series shown as mean +- s.d. (v1 plotted every alloy
    under a '1 alloy' label, written before seeds 2-5 existed; 2026-10-16). v1 is kept."""
    acta = mf.read_csv("sheet_md.csv")
    dd = mf.read_csv("md_decomposition.csv")
    m = np.abs(dd["x_cr"] - 0.5) < 1e-9
    phs = np.unique(np.round(dd["phi"][m], 6))
    B = [(dd["c11"][m & np.isclose(dd["phi"], p)] + 2 * dd["c12"][m & np.isclose(dd["phi"], p)]) / 3 for p in phs]
    c11_1 = dd["c11"][m & np.isclose(dd["phi"], 1.0)]
    fig, (a, b) = plt.subplots(1, 2, figsize=(5.5, 2.5))
    a.plot(acta["phi"], acta["bulk_gpa"], "o", color="0.15", ms=3.5, zorder=4, label=r"Acta, 50$^3$ (10 alloys)")
    a.errorbar(phs, [v.mean() for v in B], yerr=[v.std(ddof=1) for v in B], fmt="s", mfc="white", mec=mf.GREEN,
               ecolor=mf.GREEN, mew=1.2, ms=6, capsize=2, label=r"present, 25$^3$ (5 alloys, mean $\pm$ s.d.)")
    a.set(xlabel=r"Spinodal conversion fraction, $\phi$", ylabel="Bulk modulus (GPa), Fe$_{50}$Cr$_{50}$",
          xlim=(0.47, 1.03), ylim=(193.5, 199.5))
    a.legend(loc="upper left")
    a.text(0.98, 0.05, f"$C_{{11}}(\\phi=1)$: {acta['c11_phi1_gpa'][0]:.1f} vs {c11_1.mean():.1f} GPa",
           transform=a.transAxes, ha="right", va="bottom")
    mk.style(a, "(a)")
    r = mf.read_csv("md_random_alloys.csv")
    b.errorbar(r["x_cr"], r["c11"], yerr=r["c11_sd"], fmt="s", color=mf.GREEN, ms=4, capsize=2, label="present, random alloys")
    b.plot([0.5], [252.7], "o", mfc="none", mec="0.15", mew=1.2, ms=8, label=r"Acta, 50/50 random ($\phi$ = 0.5)")
    b.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"MD $C_{11}$ (GPa)", xlim=(-0.03, 1.03), ylim=(230, 340))
    b.legend(loc="upper left")
    mk.style(b, "(b)")
    return finish(fig, "R_md_validation_v2")


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "v2":
    md_validation_v2()
