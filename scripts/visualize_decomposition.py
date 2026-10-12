#!/usr/bin/env python3
"""Snapshots of the MD spinodal-decomposition structures (Lawrence's workflow, 25^3 bcc cells).

For every lammps/decomposition/runs/x*_L*_a*/phi*.data, one (001) atomic plane is drawn with
Fe and Cr coloured. Outputs in visualizations/:
  strip_x{X}.png   six panels phi = 0.5 ... 1.0 for one composition
  montage.png      rows = compositions, columns = phi
"""
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "lammps" / "decomposition" / "runs"
OUT = ROOT / "visualizations"
FE, CR = "#C0C0C8", "#2E5E9E"   # light grey Fe, blue Cr


def read_data(path):
    txt = path.read_text()
    lo_hi = [list(map(float, re.search(rf"([-\d.eE+]+) ([-\d.eE+]+) {a}lo {a}hi", txt).groups())) for a in "xyz"]
    body = txt.split("Atoms")[1].split("\n\n", 2)[1]
    arr = np.array([l.split()[:5] for l in body.strip().split("\n")], float)
    types, xyz = arr[:, 1].astype(int), arr[:, 2:5]
    return types, xyz, np.array(lo_hi)


def plane(types, xyz, box, layer=0):
    a = (box[0, 1] - box[0, 0]) / round((box[0, 1] - box[0, 0]) / 2.866)   # bcc lattice parameter
    z0 = box[2, 0] + layer * a
    m = np.abs(((xyz[:, 2] - z0) + a / 4) % a - a / 4) < a / 8          # atoms on the z = z0 (mod a) plane
    return types[m], xyz[m], a


def draw(ax, types, xyz, box):
    t, p, a = plane(types, xyz, box)
    ax.scatter(p[:, 0], p[:, 1], s=9, c=np.where(t == 2, CR, FE), marker="s", lw=0)
    ax.set_xlim(box[0]); ax.set_ylim(box[1]); ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    return np.mean(types == 2)


def series():
    out = {}
    for d in sorted(RUNS.glob("x*_L*_a*")):
        x = float(re.match(r"x([\d.]+)_", d.name).group(1))
        files = sorted(d.glob("phi*.data"), key=lambda p: float(p.stem[3:]))
        if len(files) == 6:
            out[x] = files
    return dict(sorted(out.items()))


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8})
    allser = series()
    for x, files in allser.items():
        fig, axs = plt.subplots(1, 6, figsize=(9, 1.9))
        for ax, f in zip(axs, files):
            types, xyz, box = read_data(f)
            draw(ax, types, xyz, box)
            ax.set_title(rf"$\phi$ = {float(f.stem[3:]):.1f}", fontsize=8)
        fig.suptitle(f"Fe$_{{{100 - round(100 * x)}}}$Cr$_{{{round(100 * x)}}}$: random (left) to fully "
                     f"decomposed (right); one (001) plane of a 25$^3$ bcc cell (grey Fe, blue Cr)",
                     fontsize=8, y=1.02)
        fig.savefig(OUT / f"strip_x{x:.2f}.png", dpi=220, bbox_inches="tight")
        plt.close(fig)
    if allser:
        fig, axs = plt.subplots(len(allser), 6, figsize=(9, 1.55 * len(allser)), squeeze=False)
        for r, (x, files) in enumerate(allser.items()):
            for c, f in enumerate(files):
                types, xyz, box = read_data(f)
                draw(axs[r, c], types, xyz, box)
                if r == 0:
                    axs[r, c].set_title(rf"$\phi$ = {float(f.stem[3:]):.1f}", fontsize=8)
            axs[r, 0].set_ylabel(f"{round(100 * x)}% Cr", fontsize=8)
        fig.suptitle("MD spinodal decomposition (Waldron/Acta workflow): grey Fe, blue Cr", fontsize=9, y=1.0)
        fig.savefig(OUT / "montage.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
    print("compositions:", list(allser))


if __name__ == "__main__":
    main()
