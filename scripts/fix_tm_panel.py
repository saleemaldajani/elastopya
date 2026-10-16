#!/usr/bin/env python3
"""Corrected SM Fig. S5b (melting temperature vs x_Cr), drawn in the style of the original Excel panel
(Parts/Figures/TMvsX_v2.png); the original is kept.

The working sheet ('Analyzed data'!BP, BT) computed the Lindemann melting temperature
    T_m = 0.00321 f^2 M V^(2/3) Theta_D^2          (M g/mol, V cm^3/mol; 0.00321 = k_B m_u (1e-6/N_A)^(2/3) / (9 hbar^2) N_A...)
with a "Debye temperature" Theta_D = 2 hbar v_m / (h a0), where v_m was the ARITHMETIC mean of v_p, v_phi and v_s,
h stood in place of k_B and a0 was in angstrom, and with Lindemann constants f (Fe 0.0915, Cr 0.0767) of no stated source.
The prefactor 0.00321 is correct (Lindemann with r = V_at^(1/3)). Correction (2026-10-15):
  * Theta_D = (hbar / k_B)(6 pi^2 n)^(1/3) v_D with the Debye (Anderson) average v_D^-3 = (2 v_s^-3 + v_l^-3)/3, VRH moduli;
  * f calibrated on the measured melting points of bcc Fe (1811 K) and Cr (2180 K) with the same Theta_D construction,
    linear in x in between (standard practice for the Lindemann constant).
Same inputs as the original panel: the VASP tensors at x = 0, 0.2, 0.5, 0.95, 1 and the sheet densities, molar masses and
molar volumes. Error bars keep the relative size of the thermal panels (+-2.7 %).

  fix_tm_panel.py  -> manuscript/Parts/Figures/TMvsX_v3.png (+ .pdf); prints both versions
"""
from pathlib import Path

import numpy as np

from fix_thermal_panels import VASP, REL_ERR, OUT

HB, H, KB, NA = 1.054571817e-34, 6.62607e-34, 1.380649e-23, 6.02214076e23
TM_FE, TM_CR = 1811.0, 2180.0


def props(x, c11, c12, c44):
    rho = (x * 7.2 + (1 - x) * 7.87) * 1000; M = x * 52.0 + (1 - x) * 55.85
    V = x * 7.888423167 + (1 - x) * 7.043993577; a0 = x * 2.8846 + (1 - x) * 2.8665
    b = (c11 + 2 * c12) / 3; gv = (c11 - c12 + 3 * c44) / 5; gr = 5 * (c11 - c12) * c44 / (4 * c44 + 3 * (c11 - c12)); g = (gv + gr) / 2
    vs, vl, vphi = (np.sqrt(m * 1e9 / rho) for m in (g, b + 4 * g / 3, b))
    th_sheet = 2 * HB * (vl + vphi + vs) / 3 / (H * a0)
    vD = (1 / 3 * (2 / vs ** 3 + 1 / vl ** 3)) ** (-1 / 3)
    thD = HB / KB * (6 * np.pi ** 2 * rho / (M / 1000) * NA) ** (1 / 3) * vD
    return M, V, th_sheet, thD


def tm(f, M, V, th):
    return 0.00321 * f ** 2 * M * V ** (2 / 3) * th ** 2


def values():
    p = {x: props(x, *t) for x, t in VASP.items()}
    fFe = np.sqrt(TM_FE / tm(1, *p[0.0][:2], p[0.0][3])); fCr = np.sqrt(TM_CR / tm(1, *p[1.0][:2], p[1.0][3]))
    rows = []
    for x, (M, V, ths, thD) in p.items():
        f_old = x * 0.0767 + (1 - x) * 0.0915; f_new = x * fCr + (1 - x) * fFe
        rows.append((x, tm(f_old, M, V, ths), thD, tm(f_new, M, V, thD)))
    return np.array(rows), fFe, fCr


def panel(fname, x, y):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    plt.rcParams.update({"font.family": "Arial", "font.size": 11, "axes.linewidth": 0.9})
    fig = plt.figure(figsize=(1828 / 330, 1943 / 330), dpi=330)
    ax = fig.add_axes([0.162, 0.150, 0.81, 0.82])
    ax.errorbar(x, y, yerr=REL_ERR * y, fmt="o", ms=12, mfc="white", mec="black", mew=1.1,
                ecolor="0.45", elinewidth=0.8, capsize=3, capthick=0.8, zorder=3)
    ax.plot([], [], "o", ms=12, mfc="white", mec="black", mew=1.1, label="T_m, melting temperature [K]")
    ax.set_xlim(0, 1); ax.set_ylim(1400, 2500)
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0]); ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
    ax.xaxis.set_minor_locator(MultipleLocator(0.125))
    ax.yaxis.set_major_locator(MultipleLocator(300)); ax.yaxis.set_minor_locator(MultipleLocator(150))
    ax.tick_params(which="major", direction="out", length=5, width=0.8, top=False, right=False, labelsize=12, pad=5)
    ax.tick_params(which="minor", direction="out", length=3, width=0.8, top=False, right=False)
    ax.set_xlabel("X_Cr", fontsize=14, fontweight="bold", labelpad=4)
    ax.set_ylabel("Melting Temperature [K]", fontsize=14, fontweight="bold", labelpad=6)
    leg = ax.legend(loc="upper left", frameon=True, fancybox=False, edgecolor="black", handletextpad=0.1,
                    borderpad=1.3, markerscale=0.55, prop={"weight": "bold", "size": 12})
    leg.get_frame().set_linewidth(0.8)
    fig.savefig(OUT / f"{fname}.png", dpi=330); fig.savefig(OUT / f"{fname}.pdf")
    plt.close(fig)
    print("wrote", OUT / fname)


if __name__ == "__main__":
    v, fFe, fCr = values()
    print(f"calibrated Lindemann f: Fe {fFe:.4f}, Cr {fCr:.4f} (sheet: 0.0915, 0.0767)")
    print("x, T_m sheet (reproduced), Theta_D Debye [K], T_m corrected [K]")
    print(np.round(v, 0))
    panel("TMvsX_v3", v[:, 0], v[:, 3])
