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


def structural_convergence():
    """R_structural_convergence (2026-10-16, Mike's convergence comment; also SM, with Table S1): (a) energy of every
    QE relaxation relative to its final value against ionic step, (b) residual total force at the last ionic step.
    Parsed from dft_qe/engaging/runs/*/relax*.out by make_convergence_table (resumed runs: parts concatenated)."""
    import re
    import make_convergence_table as ct
    sets = []
    for lat, pat, extra in (("relaxed", "sqs40_x*_k442_e60", "sqs40_x0.80_k442_e60_crafm2"),
                            ("exp", "sqs40_x*_k442_e60_aexp", "sqs40_x0.80_k442_e60_aexp_crafm")):
        for d in sorted(ct.RUNS.glob(pat)) + [ct.RUNS / extra]:
            if not ((d / "result.json").exists() or (d / "relax_result.json").exists() ):
                continue
            if ct.row(d) is None:
                continue
            parts = sorted(d.glob("relax_part*.out"), key=lambda q: int(q.stem[10:])) + [d / "relax.out"]
            e = [float(v) for p in parts for v in re.findall(r"^!\s+total energy\s+=\s+([-\d.]+)", ct.last_bfgs(p.read_text()), re.M)]
            x = float(re.match(r"sqs40_x([\d.]+)_", d.name).group(1))
            sets.append((lat, x, np.array(e), ct.row(d)[0]))
    fig, (a, b) = plt.subplots(1, 2, figsize=(5.5, 2.6), gridspec_kw={"width_ratios": [1.35, 1]})
    cmap = plt.get_cmap("viridis")
    for lat, x, e, f in sets:
        de = (e - e[-1]) / 40 * 13.605693e3          # meV/atom
        n = np.arange(1, len(e))
        a.semilogy(n, np.abs(de[:-1]) + 1e-6, "-" if lat == "relaxed" else "--", color=cmap(x), lw=0.8)
    a.axhline(1e-5 / 40 * 13.605693e3, color="0.4", lw=0.7, ls=":")
    a.text(0.97, 1e-5 / 40 * 13.605693e3 / 1.3, r"$10^{-5}$ Ry/cell", transform=a.get_yaxis_transform(), ha="right",
           va="top", color="0.35")
    a.plot([], [], "-", color="0.3", lw=0.8, label="PBE-relaxed cell")
    a.plot([], [], "--", color="0.3", lw=0.8, label="experimental $a$ (ions only)")
    a.set(xlabel="Ionic (BFGS) step", ylabel=r"$|E-E_{\rm final}|$ (meV/atom)", ylim=(1e-3, 30), xlim=(0.5, 24))
    a.legend(loc="upper right")
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
    cb = fig.colorbar(sm, ax=a, pad=0.02, fraction=0.05); cb.set_label(r"$x_{\rm Cr}$")
    mk.style(a, "(a)")
    for lat, mk_, col, lab in (("relaxed", "o", mf.RED, "PBE-relaxed cell"), ("exp", "D", mf.BLUE, "experimental $a$")):
        xs = [x for l, x, _, _ in sets if l == lat]; fs = [f for l, _, _, f in sets if l == lat]
        b.plot(xs, fs, mk_, mfc="white", mec=col, mew=1.0, ms=4.5, ls="none", label=lab)
    b.set(xlabel=r"$x_{\rm Cr}$", ylabel="Residual total force (Ry/bohr)", xlim=(-0.03, 1.03), ylim=(0, 0.004))
    b.legend(loc="upper left")
    mk.style(b, "(b)")
    return finish(fig, "R_structural_convergence")


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "structural":
    apply_publication_style("scripta", palette_name="okabe_ito_on_white")
    plt.rcParams["svg.fonttype"] = "none"
    print(structural_convergence())


def _trunc(name, lo, hi):
    """Colour map without its palest end (the bright yellow of viridis/plasma is hard to see on white)."""
    from matplotlib.colors import LinearSegmentedColormap
    base = plt.get_cmap(name)
    return LinearSegmentedColormap.from_list(f"{name}_{lo}_{hi}", base(np.linspace(lo, hi, 256)))


def fig3b_options():
    """R_fig3b_options (2026-10-16): three candidate Fig. 3b panels with all twelve MD compositions, for Mike's choice.
    (a) absolute f_SAW vs phi, x_Cr colour bar; (b) relative change vs phi, Fe-rich solid / Cr-rich dashed;
    (c) relative change vs x_Cr at every phi, phi colour bar. Exploration in scripts/explore_fig3b_options.py."""
    dec = mf.md_decomposition(); acta = mf.md(); xs = sorted(dec)
    cx = _trunc("viridis", 0.0, 0.82); nx = plt.Normalize(0, 1)
    cp = _trunc("plasma", 0.0, 0.78); npn = plt.Normalize(0.6, 1.0)
    rel = lambda s: 100 * (s["f"] / s["f"][np.argmin(s["phi"])] - 1)
    fig = plt.figure(figsize=(7.4, 2.7))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.12, 1.12], wspace=0.48)
    a, b, c = (fig.add_subplot(gs[0, i]) for i in range(3))
    ticks = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    for x in xs:
        s = dec[x]; o = np.argsort(s["phi"]); i0 = np.argmin(s["phi"]); col = cx(nx(x))
        sd = s["f_sd"][o] if s.get("f_sd") is not None else None
        a.errorbar(s["phi"][o], s["f"][o], yerr=sd, fmt="-s", color=col, ms=2.6, lw=0.8, elinewidth=0.5, capsize=1)
        fe = x <= 0.5 + 1e-9
        b.errorbar(s["phi"][o], rel(s)[o], yerr=None if sd is None else 100 * sd / s["f"][i0], fmt="-o" if fe else "--s",
                   color=col, mfc=col if fe else "white", ms=2.6, lw=0.8, elinewidth=0.5, capsize=1)
    a.plot(acta["phi"], acta["f"], "-o", color="black", ms=3.2, lw=0.8, zorder=5, label=r"Fe$_{50}$Cr$_{50}$, Acta 50$^3$")
    a.plot([], [], "-s", color="0.45", ms=2.6, label=r"present, 25$^3$ (5 alloys)")
    a.set(xlabel=r"Spinodal conversion fraction, $\phi$", ylabel=r"$f_{\rm SAW}$ (GHz), $\lambda$ = 6.4 $\mu$m",
          xlim=(0.47, 1.03), ylim=(0.440, 0.545))
    ra = 100 * (acta["f"] / acta["f"][np.argmin(acta["phi"])] - 1)
    b.plot(acta["phi"], ra, "o", mfc="none", mec="black", ms=4.5, mew=0.9, zorder=5, label=r"Fe$_{50}$Cr$_{50}$, Acta 50$^3$")
    b.plot([], [], "-o", color="0.45", ms=2.6, label=r"Fe-rich, 25$^3$")
    b.plot([], [], "--s", color="0.45", mfc="white", ms=2.6, label=r"Cr-rich, 25$^3$")
    b.set(xlabel=r"Spinodal conversion fraction, $\phi$", ylabel=r"$\Delta f_{\rm SAW}/f_{\rm SAW}(\phi=0.5)$ (%)",
          xlim=(0.47, 1.03), ylim=(-0.5, 8.0))
    for ax in (a, b):
        ax.set_xticks(ticks); ax.legend(loc="upper left")
    cbx = fig.colorbar(plt.cm.ScalarMappable(cmap=cx, norm=nx), ax=b, pad=0.03, fraction=0.07, aspect=22)
    cbx.set_label(r"$x_{\rm Cr}$")
    c.axvspan(0.30, 0.71, color="0.88", lw=0, zorder=0, label="spinodal, 748 K")
    for p in ticks[1:]:
        y, e = [], []
        for x in xs:
            s = dec[x]; i = int(np.argmin(abs(s["phi"] - p))); i0 = int(np.argmin(s["phi"]))
            y.append(100 * (s["f"][i] / s["f"][i0] - 1))
            e.append(100 * s["f_sd"][i] / s["f"][i0] if s.get("f_sd") is not None else 0)
        c.errorbar(xs, y, yerr=e, fmt="-o", color=cp(npn(p)), ms=2.6, lw=0.8, elinewidth=0.5, capsize=1)
        i = int(np.argmin(abs(acta["phi"] - p)))
        c.plot([0.5], [ra[i]], "o", mfc="none", mec=cp(npn(p)), ms=4.5, mew=0.9, zorder=5)
    c.plot([], [], "o", mfc="none", mec="0.3", ms=4.5, ls="none", label=r"Fe$_{50}$Cr$_{50}$, Acta 50$^3$")
    c.plot([], [], "-o", color="0.3", ms=2.6, label=r"present, 25$^3$ (5 alloys)")
    c.plot([0.36], [-0.25], "*", color=mf.BLUE, ms=6, zorder=6)
    c.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"$\Delta f_{\rm SAW}/f_{\rm SAW}(\phi=0.5)$ (%)", xlim=(0, 1), ylim=(-0.5, 8.0))
    c.legend(loc="upper right")
    cbp = fig.colorbar(plt.cm.ScalarMappable(cmap=cp, norm=npn), ax=c, pad=0.03, fraction=0.07, aspect=22, ticks=ticks[1:])
    cbp.set_label(r"$\phi$")
    for ax, t in ((a, "(a)"), (b, "(b)"), (c, "(c)")):
        mk.style(ax, t)
    fig.set_size_inches(mk.PROFILE["widths_mm"]["full"] / 25.4, 2.7 * mk.PROFILE["widths_mm"]["full"] / 25.4 / 7.4)
    mk.normalize_text(fig)
    stem = "R_fig3b_options"
    mk.export_figure(fig, OUT / stem, formats=["pdf", "png", "svg"], dpi=mk.DPI, bbox_inches="tight", overwrite=True,
                     mkdir=True, write_manifest=True, provenance={"source": "scripts/make_response_figures.py"})
    found, sizes = mk.audit(fig)
    (OUT / f"{stem}.audit.json").write_text(json.dumps(dict(issues=found, font_sizes_pt=sizes), indent=1) + "\n")
    plt.close(fig)
    return stem, found, min(sizes)


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "fig3b":
    apply_publication_style("scripta", palette_name="okabe_ito_on_white")
    plt.rcParams["svg.fonttype"] = "none"
    print(fig3b_options())


def fig3b_options_v2():
    """R_fig3b_options_v2 (2026-10-16): the current Fig. 3b as (a) next to the three candidates (b-d) of
    R_fig3b_options (kept), in a 2 x 2 layout so all four can be compared at once."""
    dec = mf.md_decomposition(); acta = mf.md(); xs = sorted(dec)
    cx = _trunc("viridis", 0.0, 0.82); nx = plt.Normalize(0, 1)
    cp = _trunc("plasma", 0.0, 0.78); npn = plt.Normalize(0.6, 1.0)
    rel = lambda s: 100 * (s["f"] / s["f"][np.argmin(s["phi"])] - 1)
    sdv = lambda s, o: s["f_sd"][o] if s.get("f_sd") is not None else None
    ticks = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    fig, axs = plt.subplots(2, 2, figsize=(7.4, 5.6))
    (a, b), (c, d) = axs
    # (a) current Fig. 3b
    a.plot(acta["phi"], acta["f"], "-", color="0.6", lw=0.8, zorder=1)
    a.plot(acta["phi"], acta["f"], "o", color="black", ms=4, ls="none", zorder=3, label=r"Fe$_{50}$Cr$_{50}$ (Acta, 50$^3$ cell)")
    cur = {0.05: mf.YELLOW, 0.2: mf.BLUE, 0.36: mf.PURPLE, 0.95: mf.ORANGE}
    for x in xs:
        s = dec[x]; o = np.argsort(s["phi"])
        if abs(x - 0.5) < 1e-9:
            a.plot(s["phi"][o], s["f"][o], "o", mfc="none", mec="0.55", ms=5, ls="none", zorder=2,
                   label=r"Fe$_{50}$Cr$_{50}$, 25$^3$ cell (validation)")
        elif round(x, 2) in cur:
            col = cur[round(x, 2)]
            a.plot(s["phi"][o], s["f"][o], "-", color=col, alpha=0.5, lw=0.8, zorder=1)
            a.errorbar(s["phi"][o], s["f"][o], yerr=sdv(s, o), fmt="s", color=col, ms=3.5, ls="none", zorder=3,
                       elinewidth=0.6, capsize=1.5, label=rf"Fe$_{{{100 - round(100 * x)}}}$Cr$_{{{round(100 * x)}}}$, 25$^3$ cell")
    h, l = a.get_legend_handles_labels()
    import re
    key = lambda t: (int(re.search(r"Cr\$_\{(\d+)\}\$", t).group(1)), "validation" in t)
    o = sorted(range(len(l)), key=lambda i: key(l[i]))
    a.legend([h[i] for i in o], [l[i] for i in o], loc="upper left")
    a.set(xlabel=r"Spinodal conversion fraction, $\phi$", ylabel=r"$f_{\rm SAW}$ (GHz), $\lambda$ = 6.4 $\mu$m",
          xlim=(0.47, 1.03), ylim=(0.445, 0.605))
    # (b) absolute, all twelve
    for x in xs:
        s = dec[x]; o = np.argsort(s["phi"])
        b.errorbar(s["phi"][o], s["f"][o], yerr=sdv(s, o), fmt="-s", color=cx(nx(x)), ms=2.6, lw=0.8, elinewidth=0.5, capsize=1)
    b.plot(acta["phi"], acta["f"], "-o", color="black", ms=3.2, lw=0.8, zorder=5, label=r"Fe$_{50}$Cr$_{50}$, Acta 50$^3$")
    b.plot([], [], "-s", color="0.45", ms=2.6, label=r"present, 25$^3$ (5 alloys)")
    b.set(xlabel=r"Spinodal conversion fraction, $\phi$", ylabel=r"$f_{\rm SAW}$ (GHz), $\lambda$ = 6.4 $\mu$m",
          xlim=(0.47, 1.03), ylim=(0.440, 0.545))
    # (c) relative vs phi
    for x in xs:
        s = dec[x]; o = np.argsort(s["phi"]); i0 = np.argmin(s["phi"]); fe = x <= 0.5 + 1e-9; col = cx(nx(x))
        sd = sdv(s, o)
        c.errorbar(s["phi"][o], rel(s)[o], yerr=None if sd is None else 100 * sd / s["f"][i0], fmt="-o" if fe else "--s",
                   color=col, mfc=col if fe else "white", ms=2.6, lw=0.8, elinewidth=0.5, capsize=1)
    ra = 100 * (acta["f"] / acta["f"][np.argmin(acta["phi"])] - 1)
    c.plot(acta["phi"], ra, "o", mfc="none", mec="black", ms=4.5, mew=0.9, zorder=5, label=r"Fe$_{50}$Cr$_{50}$, Acta 50$^3$")
    c.plot([], [], "-o", color="0.45", ms=2.6, label=r"Fe-rich, 25$^3$")
    c.plot([], [], "--s", color="0.45", mfc="white", ms=2.6, label=r"Cr-rich, 25$^3$")
    c.set(xlabel=r"Spinodal conversion fraction, $\phi$", ylabel=r"$\Delta f_{\rm SAW}/f_{\rm SAW}(\phi=0.5)$ (%)",
          xlim=(0.47, 1.03), ylim=(-0.5, 8.0))
    for ax in (a, b, c):
        ax.set_xticks(ticks)
    for ax in (b, c):
        ax.legend(loc="upper left")
    # (d) relative vs x_Cr at every phi
    d.axvspan(0.30, 0.71, color="0.88", lw=0, zorder=0, label="spinodal, 748 K")
    for p in ticks[1:]:
        y, e = [], []
        for x in xs:
            s = dec[x]; i = int(np.argmin(abs(s["phi"] - p))); i0 = int(np.argmin(s["phi"]))
            y.append(100 * (s["f"][i] / s["f"][i0] - 1))
            e.append(100 * s["f_sd"][i] / s["f"][i0] if s.get("f_sd") is not None else 0)
        d.errorbar(xs, y, yerr=e, fmt="-o", color=cp(npn(p)), ms=2.6, lw=0.8, elinewidth=0.5, capsize=1)
        i = int(np.argmin(abs(acta["phi"] - p)))
        d.plot([0.5], [ra[i]], "o", mfc="none", mec=cp(npn(p)), ms=4.5, mew=0.9, zorder=5)
    d.plot([], [], "o", mfc="none", mec="0.3", ms=4.5, ls="none", label=r"Fe$_{50}$Cr$_{50}$, Acta 50$^3$")
    d.plot([], [], "-o", color="0.3", ms=2.6, label=r"present, 25$^3$ (5 alloys)")
    d.plot([0.36], [-0.25], "*", color=mf.BLUE, ms=6, zorder=6)
    d.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"$\Delta f_{\rm SAW}/f_{\rm SAW}(\phi=0.5)$ (%)", xlim=(0, 1), ylim=(-0.5, 8.0))
    d.legend(loc="upper left")
    kw = dict(pad=0.03, fraction=0.06, aspect=22)
    fig.colorbar(plt.cm.ScalarMappable(cmap=cx, norm=nx), ax=a, **kw).ax.set_visible(False)   # keeps (a) the width of (b)
    fig.colorbar(plt.cm.ScalarMappable(cmap=cx, norm=nx), ax=b, **kw).set_label(r"$x_{\rm Cr}$")
    fig.colorbar(plt.cm.ScalarMappable(cmap=cx, norm=nx), ax=c, **kw).set_label(r"$x_{\rm Cr}$")
    fig.colorbar(plt.cm.ScalarMappable(cmap=cp, norm=npn), ax=d, ticks=ticks[1:], **kw).set_label(r"$\phi$")
    for ax, t in ((a, "(a) current"), (b, "(b)"), (c, "(c)"), (d, "(d)")):
        mk.style(ax, t)
    fig.subplots_adjust(wspace=0.38, hspace=0.42)
    fig.set_size_inches(mk.PROFILE["widths_mm"]["full"] / 25.4, 5.6 * mk.PROFILE["widths_mm"]["full"] / 25.4 / 7.4)
    mk.normalize_text(fig)
    stem = "R_fig3b_options_v2"
    mk.export_figure(fig, OUT / stem, formats=["pdf", "png", "svg"], dpi=mk.DPI, bbox_inches="tight", overwrite=True,
                     mkdir=True, write_manifest=True, provenance={"source": "scripts/make_response_figures.py"})
    found, sizes = mk.audit(fig)
    (OUT / f"{stem}.audit.json").write_text(json.dumps(dict(issues=found, font_sizes_pt=sizes), indent=1) + "\n")
    plt.close(fig)
    return stem, found, min(sizes)
