#!/usr/bin/env python3
"""SM Fig. S_sdw_cr_v2: the computed Cr SDW supercell (rows a-c, as S_sdw_cr) plus row (d), the chain from the
electrons to the SAW in numbers (added 2026-10-15; S_sdw_cr is kept unchanged):
  (d, left)  magnetic pressure P_mag = P(magnetic) - P(NM) of the AF and SDW cells against <m^2>: a line through
             the origin is the Landau magnetovolume coupling P_mag ~ <m^2>
  (d, right) what AF order does to the moduli and the SAW of pure Cr (QE, 2-atom cells): K drops, G and f_SAW barely move
The strain-wave amplitude (static continuum response to the 2Q force harmonic) is printed in the caption numbers.

  make_sdw_v2.py  -> manuscript/figures/kdense_scripta/supplementary/S_sdw_cr_v2.{png,pdf,svg}
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent / "dft_qe"))
import make_fig4_mechanism as m4  # noqa: E402  (sdw_numbers)
import sdw_cr  # noqa: E402

OUT = HERE.parent / "manuscript" / "figures" / "kdense_scripta" / "supplementary" / "S_sdw_cr_v2"
ENG = HERE.parent / "dft_qe" / "engaging" / "sdw"


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(HERE / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    n = m4.sdw_numbers()
    res = {k: sdw_cr.sites(ENG / k / "scf.out") for k in ("nm", "afm", "sdw")}
    z = np.arange(sdw_cr.NL) * sdw_cr.A / 2
    fig = plt.figure(figsize=(5.4, 6.6))
    gs = fig.add_gridspec(4, 2, height_ratios=[1, 0.9, 0.9, 1.05], hspace=0.42, wspace=0.42)
    a0 = fig.add_subplot(gs[0, :]); a1 = fig.add_subplot(gs[1, :], sharex=a0); a2 = fig.add_subplot(gs[2, :], sharex=a0)
    ch, mg, fz, e = res["sdw"]
    a0.plot(z, mg, "o-", color="#0072B2", ms=2.6, lw=0.8, label="SDW")
    a0.plot(z, np.abs(mg), "-", color="0.6", lw=0.6, label="|m|")
    a0.plot(z, res["afm"][1], ".", color="0.4", ms=2, label="commensurate AFM")
    a0.set_ylim(-1.35, 2.0); a0.set_ylabel(r"$m$ ($\mu_{\rm B}$)")
    a0.legend(fontsize=6, ncol=3, loc="upper center", frameon=False)
    de = lambda k: (res[k][3] - res["nm"][3]) / sdw_cr.NL * 13.605693e3
    a0.set_title(f"(a) SDW moments.  E(AFM) $-$ E(NM) = {de('afm'):.1f}, E(SDW) $-$ E(NM) = {de('sdw'):.1f} meV/atom",
                 fontsize=6.5, loc="left")
    a1.plot(z, 1e3 * (ch - ch.mean()), "s-", color="#D55E00", ms=2.6, lw=0.8)
    a1.set_ylim(-7.5, 6.5); a1.set_ylabel(r"$\Delta q$ (10$^{-3}$ e)")
    a1.set_title(r"(b) site charge, SDW $-$ mean: CDW with period $\lambda_{\rm SDW}/2$", fontsize=6.5, loc="left")
    a2.plot(z, 1e3 * fz * 13.605693 / 0.529177, "^-", color="#009E73", ms=2.6, lw=0.8)
    a2.set_ylim(-75, 95); a2.set_ylabel(r"$F_z$ (meV/$\AA$)"); a2.set_xlabel(r"$z$ ($\AA$) along [001]")
    a2.set_title(f"(c) force along Q: drives a strain wave, $u\\approx${n['u_A']:.2f} $\\AA$, "
                 f"$\\varepsilon\\approx${100 * n['strain']:.1f}% (period {n['period_A']:.0f} $\\AA$)", fontsize=6.5, loc="left")
    for ax in (a0, a1, a2):
        ax.axhline(0, color="0.7", lw=0.5)
    plt.setp(a0.get_xticklabels(), visible=False); plt.setp(a1.get_xticklabels(), visible=False)
    # (d, left): magnetovolume coupling
    b0 = fig.add_subplot(gs[3, 0])
    xs = np.array([0.0, n["m2_sdw"], n["m2_afm"]]); ys = np.array([0.0, n["Pmag_sdw"], n["Pmag_afm"]])
    slope = float(np.sum(xs * ys) / np.sum(xs * xs))
    xx = np.linspace(0, 1.45, 10)
    b0.plot(xx, slope * xx, "--", color="0.5", lw=0.7)
    for x, y, lab, col in zip(xs, ys, ("NM", "SDW", "AF"), ("0.5", "#0072B2", "0.15")):
        b0.plot(x, y, "o", color=col, ms=4)
        b0.annotate(lab, (x, y), xytext=(4, -9 if lab == "AF" else 3), textcoords="offset points", fontsize=6)
    b0.set(xlabel=r"$\langle m^2\rangle$ ($\mu_{\rm B}^2$)", ylabel=r"$P_{\rm mag}$ (GPa)", xlim=(-0.05, 1.5), ylim=(-0.4, 7))
    b0.set_title(r"(d) magnetovolume: $P_{\rm mag}\propto\langle m^2\rangle$", fontsize=6.5, loc="left")
    b0.text(0.03, 0.80, f"slope {slope:.1f} GPa/$\\mu_{{\\rm B}}^2$", transform=b0.transAxes, fontsize=5.5, color="0.3")
    # (d, right): AF order vs NM in pure Cr: moduli and SAW
    b1 = fig.add_subplot(gs[3, 1])
    labs = ["$K$", "$G$", "$f_{\\rm SAW}$"]
    vals = [100 * (n["K_afm"] / n["K_nm"] - 1), 100 * (n["G_afm"] / n["G_nm"] - 1), 100 * (n["f_afm"] / n["f_nm"] - 1)]
    cols = ["#C02F00" if abs(v) > 5 else "#0072B2" for v in vals]
    b1.bar(range(3), vals, color=cols, width=0.6)
    for i, v in enumerate(vals):
        b1.text(i, v + (-3.5 if v < 0 else 1.0), f"{v:+.1f}%".replace("-", "$-$"), ha="center", fontsize=6)
    b1.axhline(0, color="0.5", lw=0.6)
    b1.set_xticks(range(3)); b1.set_xticklabels(labs, fontsize=6.5)
    b1.set(ylabel="AF vs NM Cr (%)", ylim=(-36, 8))
    b1.set_title("AF order: volume channel only", fontsize=6.5, loc="left")
    fig.subplots_adjust(left=0.12, right=0.97, top=0.96, bottom=0.07)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT.with_suffix(f".{ext}"), dpi=600)
    print("wrote", OUT)
    print({k: round(v, 3) for k, v in n.items()}, "slope", round(slope, 2))


if __name__ == "__main__":
    main()
