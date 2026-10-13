#!/usr/bin/env python3
"""Spin-polarisation and deformation-charge maps of the relaxed 40-atom SQS (scf of dos_sqs.py), on the (110)
plane through the cell origin, spanned by (2a, 2a, 0) and (0, 0, 5a); it contains 20 of the 40 sites
(10 corner and 10 body-centre). Shows the Fe moments and the Cr moments antiparallel to them at low Cr content
and vanishing between 35 and 50 % Cr (the mechanism Razumovskiy et al. give for the elastic anomaly).

  density_maps.py inputs <dos_dir>...       pp_spin.in (plot_num 6), pp_drho.in (plot_num 9: rho - atomic)
  density_maps.py plot <outfile> <dos_dir>...
"""
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import elastic_sqs_qe as e  # noqa: E402

NX, NY = 240, 480
KINDS = {"spin": 6, "drho": 9}


def cell_and_sites(d):
    s = Path(d, "scf.in").read_text()
    rows = s.split("CELL_PARAMETERS angstrom\n")[1].split("\n")[:3]
    cell = np.array([[float(v) for v in r.split()] for r in rows])
    pos = re.findall(r"^(Fe|Cr)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)", s.split("ATOMIC_POSITIONS crystal")[1], re.M)
    return cell, [(el, np.array([float(a), float(b), float(c)]) % 1.0) for el, a, b, c in pos]


def inputs(dirs):
    for d in dirs:
        cell, _ = cell_and_sites(d)
        alat = np.linalg.norm(cell[0])
        e1 = (cell[0] + cell[1]) / alat
        e2 = cell[2] / alat
        for k, num in KINDS.items():
            Path(d, f"pp_{k}.in").write_text(
                f"&inputpp\n  prefix = 'dos', outdir = './tmp', filplot = '{k}.pp', plot_num = {num}\n/\n"
                f"&plot\n  nfile = 1, filepp(1) = '{k}.pp', weight(1) = 1.0, iflag = 2, output_format = 7,\n"
                f"  fileout = '{k}_110.dat',\n"
                f"  e1(1) = {e1[0]:.6f}, e1(2) = {e1[1]:.6f}, e1(3) = {e1[2]:.6f},\n"
                f"  e2(1) = {e2[0]:.6f}, e2(2) = {e2[1]:.6f}, e2(3) = {e2[2]:.6f},\n"
                f"  x0(1) = 0.0, x0(2) = 0.0, x0(3) = 0.0, nx = {NX}, ny = {NY}\n/\n")
        print(d, "pp inputs")


def load(d, k):
    a = np.loadtxt(Path(d, f"{k}_110.dat"))
    return a[:, 2].reshape(NY, NX) if a.shape[0] == NX * NY else a[:, 2].reshape(NX, NY).T


def plot(outfile, dirs):
    import matplotlib.pyplot as plt
    from matplotlib.colors import TwoSlopeNorm
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    n = len(dirs)
    fig, axs = plt.subplots(2, n, figsize=(1.45 * n + 0.6, 5.2), squeeze=False)
    lim = {k: max(np.percentile(np.abs(load(d, k)), 99.5) for d in dirs) for k in KINDS}
    for j, d in enumerate(dirs):
        cell, sites = cell_and_sites(d)
        L1 = np.linalg.norm(cell[0] + cell[1]); L2 = np.linalg.norm(cell[2])
        x = float(re.match(r"x([\d.]+)", Path(d).name).group(1))
        for i, (k, cmap) in enumerate((("spin", "RdBu_r"), ("drho", "PuOr_r"))):
            ax = axs[i, j]
            im = ax.imshow(load(d, k), origin="lower", extent=(0, L1, 0, L2), cmap=cmap,
                           norm=TwoSlopeNorm(0, -lim[k], lim[k]), interpolation="bilinear")
            for el, f in sites:                    # sites on the plane: f1 = f2
                if abs(((f[0] - f[1]) + 0.5) % 1.0 - 0.5) < 0.04:
                    ax.plot(f[0] * L1, f[2] * L2, "o", ms=3.2, mfc="none", mec="black" if el == "Fe" else "#009E73", mew=0.8)
            ax.set_xticks([]); ax.set_yticks([])
            if i == 0:
                ax.set_title(f"{round(100 * x)}% Cr", fontsize=8)
            if j == n - 1:
                cb = fig.colorbar(im, ax=ax, fraction=0.08, pad=0.03)
                cb.set_label(r"$\rho_\uparrow-\rho_\downarrow$ (e/bohr$^3$)" if k == "spin" else r"$\rho-\rho_{\rm atoms}$ (e/bohr$^3$)", fontsize=7)
    axs[0, 0].set_ylabel("spin density, (110) plane", fontsize=7)
    axs[1, 0].set_ylabel("deformation charge", fontsize=7)
    fig.text(0.01, 0.005, "circles: Fe (black), Cr (green) sites in the plane; [001] vertical, [110] horizontal",
             fontsize=6, color="0.3")
    fig.tight_layout(w_pad=0.3, h_pad=0.4, rect=(0, 0.02, 1, 1))
    fig.savefig(outfile, dpi=600); fig.savefig(Path(outfile).with_suffix(".pdf")); fig.savefig(Path(outfile).with_suffix(".svg"))
    print("wrote", outfile)


if __name__ == "__main__":
    if sys.argv[1] == "inputs":
        inputs(sys.argv[2:])
    elif sys.argv[1] == "plot":
        plot(sys.argv[2], sys.argv[3:])
