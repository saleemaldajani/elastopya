#!/usr/bin/env python3
"""Fe-rich (alpha) and Cr-rich (alpha') compositions in the MD decomposition structures.

Local Cr fraction of every atom = Cr share among itself and its neighbours up to the 3rd bcc shell
(27 atoms, periodic). The distribution is split at the Otsu threshold into alpha and alpha';
the phase compositions are the mean local Cr fractions on each side, and the phase fraction is
the share of atoms on the Cr-rich side.
Outputs: data/md_phase_compositions.csv, visualizations/local_composition_histograms.png
"""
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import cKDTree

sys_path = Path(__file__).resolve().parent
import sys  # noqa: E402
sys.path.insert(0, str(sys_path))
from visualize_decomposition import read_data, series  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def local_cr(types, xyz, box):
    L = box[:, 1] - box[:, 0]
    a = L[0] / round(L[0] / 2.866)
    tree = cKDTree((xyz - box[:, 0]) % L, boxsize=L)
    nbrs = tree.query_ball_point((xyz - box[:, 0]) % L, r=1.05 * a * np.sqrt(2))
    cr = (types == 2).astype(float)
    return np.array([cr[n].mean() for n in nbrs])


def otsu(v, bins=54):
    h, e = np.histogram(v, bins=bins, range=(0, 1))
    c = 0.5 * (e[1:] + e[:-1]); w = np.cumsum(h); m = np.cumsum(h * c)
    with np.errstate(divide="ignore", invalid="ignore"):
        s = (m[-1] * w / w[-1] - m) ** 2 / (w * (w[-1] - w))
    return c[np.nanargmax(s[:-1])]


def main():
    rows = []
    ser = series()
    fig, axs = plt.subplots(1, max(1, len(ser)), figsize=(2.0 * max(1, len(ser)), 1.9), squeeze=False)
    for ax, (x, files) in zip(axs[0], ser.items()):
        for f in (files[0], files[-1]):
            phi = float(f.stem[3:])
            t, p, box = read_data(f)
            lc = local_cr(t, p, box)
            thr = otsu(lc)
            lo, hi = lc[lc < thr], lc[lc >= thr]
            rows.append(dict(x_cr=x, phi=phi, x_alpha=lo.mean() if lo.size else np.nan,
                             x_alpha_prime=hi.mean() if hi.size else np.nan,
                             f_alpha_prime=hi.size / lc.size, threshold=thr))
            ax.hist(lc, bins=27, range=(0, 1), histtype="stepfilled" if phi == 1 else "step",
                    alpha=0.45 if phi == 1 else 1, color="#2E5E9E" if phi == 1 else "0.4",
                    label=rf"$\phi$ = {phi:.1f}", density=True)
        ax.set_title(f"{round(100 * x)}% Cr", fontsize=8); ax.set_xlabel("local Cr fraction", fontsize=7)
        ax.tick_params(labelsize=6); ax.set_yticks([])
    axs[0, 0].legend(frameon=False, fontsize=6)
    fig.tight_layout()
    (ROOT / "visualizations").mkdir(exist_ok=True)
    fig.savefig(ROOT / "visualizations" / "local_composition_histograms.png", dpi=220, bbox_inches="tight")
    with open(ROOT / "data" / "md_phase_compositions.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0])); w.writeheader(); w.writerows(rows)
    for r in rows:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})


if __name__ == "__main__":
    main()
