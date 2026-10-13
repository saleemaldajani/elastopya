#!/usr/bin/env python3
"""DFT cells behind the first-principles tensors, for panel (d) of the figures: the 2-atom bcc cells of
pure Fe and Cr and the 40-atom SQS (2x2x5 conventional cells) of the Fe-Cr alloys, as relaxed in Quantum
ESPRESSO (dft_qe/engaging/runs), with arrows for the local magnetic moments (QE, integrated on atomic
spheres; up = parallel to the Fe moments).

Oblique view with the long (5a) axis horizontal: u = z + k y cos35, v = x + k y sin35 (y into the page).
All cells share one scale, so the unit cells appear at their true size next to the supercells.

  visualize_sqs.py      writes visualizations/dft_cells_strip.png (all compositions)
"""
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "dft_qe"))
import elastic_sqs_qe as eq  # noqa: E402
from visualize_decomposition import CR, FE  # noqa: E402

RUNS = ROOT / "dft_qe" / "engaging" / "runs"
SQS = ROOT / "dft_qe" / "sqs"
K_DEPTH, ANG = 0.5, np.radians(35)
R_ATOM = 0.2                  # sphere radius (units of a)
M_SCALE = 0.62 / 2.2          # arrow length per muB (units of a): Fe ~2.2 muB -> 0.62 a, centred on the atom
XLIM = (-0.3, 5 + K_DEPTH * 2 * np.cos(ANG) + 0.3)
YLIM = (-0.62, 2 + K_DEPTH * 2 * np.sin(ANG) + 0.42)    # room for the moment arrows of the edge atoms
ASPECT = (YLIM[1] - YLIM[0]) / (XLIM[1] - XLIM[0])      # block height / width, so atoms stay round
# (axes use set_aspect('auto'): with this block shape the data aspect is 1 to within a few %)
# Pure elements: QE 2-atom conventional cells (Fe FM; Cr as in the DFT row of the figures, no moment drawn)
PURE = {0.0: ("Fe", "Fe_fm_k16x16x16_e60"), 1.0: ("Cr", None)}


def proj(p):
    """p: (..., 3) in units of a, (x, y, z) -> (u, v)."""
    return np.stack([p[..., 2] + K_DEPTH * p[..., 1] * np.cos(ANG), p[..., 0] + K_DEPTH * p[..., 1] * np.sin(ANG)], -1)


def site_moments(out):
    blocks = re.split(r"Magnetic moment per site", Path(out).read_text())
    if len(blocks) < 2:
        return None
    return np.array([float(m) for m in re.findall(r"magn=\s*([-\d.]+)", blocks[-1].split("\n\n")[0])])


def cell_data(x):
    """(lengths in a, fractional positions, symbols, moments or None) of the DFT cell at composition x."""
    if x in PURE:
        el, run = PURE[x]
        m = site_moments(RUNS / run / "relax.out") if run else None
        return np.ones(3), np.array([[0, 0, 0], [0.5, 0.5, 0.5]]), [el, el], m
    run = RUNS / f"sqs40_x{x:.2f}_k442_e60"
    cell0, sym, frac0 = eq.read_poscar(SQS / f"x{x:.2f}.vasp")
    cell, frac = eq.final_structure((run / "relax.out").read_text(), cell0, frac0)
    a = abs(np.linalg.det(cell) / len(sym) * 2) ** (1 / 3)
    return np.diag(cell) / a, frac, sym, site_moments(run / "relax.out")


def draw(ax, x, moments=True):
    import matplotlib.patches as mp
    L, frac, sym, mag = cell_data(x)
    if len(sym) == 2:                                                  # bcc cube: add the 7 corner images
        frac = np.vstack([np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], float), frac[1:]])
        sym = [sym[0]] * 8 + sym[1:]
        mag = None if mag is None else np.r_[[mag[0]] * 8, mag[1:]]
    p = frac * L                                                       # (x, y, z) in units of a
    corners = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], float) * L
    for i in range(8):                                                 # cell edges (back ones dashed)
        for j in range(i + 1, 8):
            if np.count_nonzero(corners[i] != corners[j]) == 1:
                hidden = np.array([0, L[1], 0])                        # back-bottom-left corner
                back = any(np.allclose(c, hidden) for c in (corners[i], corners[j]))
                uv = proj(corners[[i, j]])
                ax.plot(uv[:, 0], uv[:, 1], color="0.35", lw=0.4, zorder=1, ls=(0, (2, 1.5)) if back else "-")
    order = np.argsort(-p[:, 1] - 1e-3 * p[:, 0])                     # far (large y) first
    uv = proj(p)
    shade = lambda y: 0.75 + 0.25 * (1 - y / max(L[1], 1e-9))         # deeper atoms slightly darker
    from matplotlib.colors import to_rgb
    for z, i in enumerate(order):
        base = np.array(to_rgb(CR if sym[i] == "Cr" else FE)) * shade(p[i, 1])
        ax.add_patch(mp.Circle(uv[i], R_ATOM, facecolor=base.clip(0, 1), edgecolor="0.15", lw=0.3, zorder=10 + 2 * z))
        if moments and mag is not None and abs(mag[i]) > 0.05:
            dv = M_SCALE * mag[i]
            ax.annotate("", xy=(uv[i, 0], uv[i, 1] + dv / 2), xytext=(uv[i, 0], uv[i, 1] - dv / 2),
                        arrowprops=dict(arrowstyle="-|>", lw=0.7, color="black", mutation_scale=4, shrinkA=0, shrinkB=0),
                        zorder=11 + 2 * z)
    ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    ax.set_aspect("auto"); ax.axis("off")
    ax._strip = True


def mean_moments(x):
    _, _, sym, m = cell_data(x)
    if m is None:
        return None, None
    sym = np.array(sym)
    f = lambda el: float(m[sym == el].mean()) if (sym == el).any() else None
    return f("Fe"), f("Cr")


QE_PURE = {0.0: "Fe_fm_k16x16x16_e60", 1.0: "Cr_nm_k16x16x16_e60"}   # 1.0: NM, as drawn


def qe_tensor(x):
    """QE (C11, C12, C44) in GPa at x (2-atom cells for the pure elements), or None if not computed."""
    import json
    run = RUNS / (QE_PURE[x] if x in QE_PURE else f"sqs40_x{x:.2f}_k442_e60") / "result.json"
    if not run.exists():
        return None
    d = json.loads(run.read_text())
    return d["C11"], d["C12"], d["C44"]


def vrh_f(x, c11, c12, c44):
    """Isotropic (VRH) Rayleigh frequency (GHz, 6.4 um), as for the DFT points of the figures."""
    import make_figures as mf
    b, cp = (c11 + 2 * c12) / 3, (c11 - c12) / 2
    gv = (2 * cp + 3 * c44) / 5
    gr = 5 * (2 * cp) * c44 / (4 * c44 + 3 * (2 * cp))
    g = (gv + gr) / 2
    e = 9 * b * g / (3 * b + g)
    return float(mf.rayleigh_f(b, e, x * mf.RHO_CR + (1 - x) * mf.RHO_FE))


def main():
    import matplotlib.pyplot as plt
    xs = [0.0, 0.2, 0.35, 0.5, 0.65, 0.95, 1.0]
    fig, axs = plt.subplots(1, len(xs), figsize=(1.4 * len(xs), 1.4 * ASPECT + 0.3))
    for ax, x in zip(axs, xs):
        draw(ax, x)
        ax.set_title(f"{round(100 * x)}% Cr", fontsize=8)
    out = ROOT / "visualizations" / "dft_cells_strip.png"
    fig.tight_layout(); fig.savefig(out, dpi=300)
    for x in xs:
        print(x, mean_moments(x))
    print("wrote", out)


if __name__ == "__main__":
    main()
