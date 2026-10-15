#!/usr/bin/env python3
"""Corrected SM Fig. S7e,f (thermal diffusivity and thermal conductivity vs x_Cr), drawn to match the original Excel panels
(Parts/Figures/DvsX_v4.png, KvsX_v2.png) in size, font, markers, error bars, legend box and tick layout.

The original panels came from the working sheet ('Analyzed data'!CB, CF), which evaluated
    kappa = 0.87 h (2 N_A / M)^(2/3) rho^(1/6) E^(1/2) x 1e11
i.e. Planck's constant (x 1e11) in place of k_B and m = 2. Eq. (2) of SM Section S1 (Clarke 2003) is
    kappa_min = 0.87 k_B (N_A m / M)^(2/3) rho^(1/6) E^(1/2),   m = 1 atom per formula unit,
so the original kappa and D = kappa / (rho c_p) were 7.62 times too large; the composition trend is unchanged.
Same inputs as the original panels: the VASP tensors at x = 0, 0.2, 0.5, 0.95, 1 (sheet densities, molar masses, c_p).
The error bars keep the relative size of the original panels (+-2.7 %).

  fix_thermal_panels.py  -> Parts/Figures/DvsX_v5.png, KvsX_v3.png (+ .pdf)
"""
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE.parent / "manuscript" / "Parts" / "Figures"
NA, H, KB, CP = 6.02214e23, 6.62607e-34, 1.380649e-23, 460.548
VASP = {0.0: (324.434353, 169.838953, 127.7242), 0.2: (400.025537, 164.57738, 120.574293), 0.5: (347.46247, 186.14903, 123.354977),
        0.95: (499.051533, 145.051597, 99.27015), 1.0: (568.728963, 158.14884, 100.913353)}
REL_ERR = 0.027


def youngs(c11, c12, c44):
    b = (c11 + 2 * c12) / 3
    g = ((c11 - c12 + 3 * c44) / 5 + 5 * (c11 - c12) * c44 / (4 * c44 + 3 * (c11 - c12))) / 2
    return 9 * b * g / (3 * b + g)


def values():
    rows = []
    for x, t in VASP.items():
        E = youngs(*t) * 1e9
        rho = (x * 7.2 + (1 - x) * 7.87) * 1000
        M = (x * 52.0 + (1 - x) * 55.85) / 1000
        sheet = 0.87 * H * ((NA * 2) / M) ** (2 / 3) * rho ** (1 / 6) * np.sqrt(E) * 1e11
        s1 = 0.87 * KB * (NA / M) ** (2 / 3) * rho ** (1 / 6) * np.sqrt(E)
        rows.append((x, sheet, sheet / (rho * CP) * 1e6, s1, s1 / (rho * CP) * 1e6))
    return np.array(rows)


def panel(fname, x, y, ylim, ymajor, yminor, ylabel, legend, legend_loc):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator
    plt.rcParams.update({"font.family": "Arial", "font.size": 11, "axes.linewidth": 0.9})
    fig = plt.figure(figsize=(1828 / 330, 1943 / 330), dpi=330)
    ax = fig.add_axes([0.162, 0.150, 0.81, 0.82])
    ax.errorbar(x, y, yerr=REL_ERR * y, fmt="o", ms=12, mfc="white", mec="black", mew=1.1,
                ecolor="0.45", elinewidth=0.8, capsize=3, capthick=0.8, zorder=3)
    ax.plot([], [], "o", ms=12, mfc="white", mec="black", mew=1.1, label=legend)      # legend: bare circle, as in Excel
    ax.set_xlim(0, 1); ax.set_ylim(*ylim)
    ax.xaxis.set_major_locator(MultipleLocator(0.25)); ax.xaxis.set_minor_locator(MultipleLocator(0.125))
    ax.yaxis.set_major_locator(MultipleLocator(ymajor)); ax.yaxis.set_minor_locator(MultipleLocator(yminor))
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0]); ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
    ax.tick_params(which="major", direction="out", length=5, width=0.8, top=False, right=False, labelsize=12, pad=5)
    ax.tick_params(which="minor", direction="out", length=3, width=0.8, top=False, right=False)
    ax.set_xlabel("X_Cr", fontsize=14, fontweight="bold", labelpad=4)
    ax.set_ylabel(ylabel, fontsize=14, fontweight="bold", labelpad=6)
    leg = ax.legend(loc=legend_loc, frameon=True, fancybox=False, edgecolor="black", handletextpad=0.1,
                    borderpad=1.3, markerscale=0.55, prop={"weight": "bold", "size": 12})
    leg.get_frame().set_linewidth(0.8)
    fig.savefig(OUT / f"{fname}.png", dpi=330); fig.savefig(OUT / f"{fname}.pdf")
    plt.close(fig)
    print("wrote", OUT / fname)


if __name__ == "__main__":
    v = values()
    print("x, kappa_sheet, D_sheet, kappa_S1, D_S1")
    print(np.round(v, 3))
    print("ratio:", np.round(v[:, 1] / v[:, 3], 3))
    panel("DvsX_v5", v[:, 0], v[:, 4], (0.2, 0.6), 0.1, 0.05, r"Thermal Diffusivity [1e-6 square m/s]",
          "D, thermal diffusivity [1e-6 square m/s]", "upper left")
    panel("KvsX_v3", v[:, 0], v[:, 3], (1.0, 2.0), 0.5, 0.125, "Thermal Conductivity [W/mK]",
          "k, thermal conductivity [W/mK]", "upper right")
