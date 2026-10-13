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
FE, CR = "#F06060", "#6060F0"   # Acta Mater. 246 (2023) Fig. 7 / OVITO defaults: red Fe, blue Cr


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


def md_grid3d(path):
    """Cr occupancy of the bcc corner sublattice of an MD cell, as an n x n x n array [ix, iy, iz]."""
    t, p, box = read_data(path)
    L = box[:, 1] - box[:, 0]
    n = round(L[0] / 2.866); a = L[0] / n
    rel = (p - box[:, 0]) / a
    on = np.all(np.abs(rel - np.round(rel)) < 0.1, axis=1)
    ijk = np.round(rel[on]).astype(int) % n
    g = np.zeros((n, n, n), int)
    g[ijk[:, 0], ijk[:, 1], ijk[:, 2]] = (t[on] == 2).astype(int)
    return g


_C30, _S30 = np.cos(np.pi / 6), 0.5
FACE_SHADE = {"top": 1.0, "right": 0.82, "left": 0.66}


def _cube_ax(ax, g):
    """Isometric view of the cell (u = (x - y) cos30, v = (x + y) sin30 + z): the visible faces are the
    top (z = max), left (x = 0) and right (y = 0), meeting at the front corner (0, 0, 1)."""
    from matplotlib.colors import to_rgb
    from matplotlib.transforms import Affine2D
    fe, cr = np.array(to_rgb(FE)), np.array(to_rgb(CR))
    rgb = lambda face, shade: (np.where(face[..., None] == 1, cr, fe) * shade).clip(0, 1)
    n = g.shape[0]
    faces = {   # (image with rows = second coordinate, affine map of the unit square onto the face)
        "top": (g[:, :, n - 1].T, [[_C30, -_C30, 0], [_S30, _S30, 1], [0, 0, 1]]),   # (x, y) at z = 1
        "left": (g[0, :, :].T, [[-_C30, 0, 0], [_S30, 1, 0], [0, 0, 1]]),            # (y, z) at x = 0
        "right": (g[:, 0, :].T, [[_C30, 0, 0], [_S30, 1, 0], [0, 0, 1]]),            # (x, z) at y = 0
    }
    for name, (img, m) in faces.items():
        ax.imshow(rgb(img, FACE_SHADE[name]), extent=(0, 1, 0, 1), origin="lower", interpolation="nearest",
                  transform=Affine2D(np.array(m, float)) + ax.transData)
    P = lambda x, y, z: ((x - y) * _C30, (x + y) * _S30 + z)
    hexagon = [P(0, 0, 0), P(1, 0, 0), P(1, 0, 1), P(1, 1, 1), P(0, 1, 1), P(0, 1, 0), P(0, 0, 0)]
    ax.plot(*zip(*hexagon), color="0.15", lw=0.5)
    for q in (P(1, 0, 1), P(0, 1, 1), P(0, 0, 0)):
        ax.plot(*zip(P(0, 0, 1), q), color="0.15", lw=0.5)
    ax.set_xlim(-_C30 - 0.02, _C30 + 0.02); ax.set_ylim(-0.02, 2.02)
    ax.set_aspect("equal"); ax.axis("off")
    ax._strip = True


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
                     f"decomposed (right); one (001) plane of a 25$^3$ bcc cell (red Fe, blue Cr)",
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
        fig.suptitle("MD spinodal decomposition (Waldron/Acta workflow): red Fe, blue Cr", fontsize=9, y=1.0)
        fig.savefig(OUT / "montage.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
    # cube versions: whole 25^3 cell, three outer faces (same rendering as the paper strips)
    for x, files in allser.items():
        fig, axs = plt.subplots(1, 6, figsize=(9, 2.0))
        for ax, f in zip(axs, files):
            _cube_ax(ax, md_grid3d(f))
            ax.set_title(rf"$\phi$ = {float(f.stem[3:]):.1f}", fontsize=8)
        fig.suptitle(f"Fe$_{{{100 - round(100 * x)}}}$Cr$_{{{round(100 * x)}}}$: random (left) to fully "
                     f"decomposed (right); whole 25$^3$ bcc cell, three outer faces (red Fe, blue Cr)",
                     fontsize=8, y=1.02)
        fig.savefig(OUT / f"cube_strip_x{x:.2f}.png", dpi=220, bbox_inches="tight")
        plt.close(fig)
    if allser:
        fig, axs = plt.subplots(len(allser), 6, figsize=(9, 1.6 * len(allser)), squeeze=False)
        for r, (x, files) in enumerate(allser.items()):
            for c, f in enumerate(files):
                _cube_ax(axs[r, c], md_grid3d(f))
                if r == 0:
                    axs[r, c].set_title(rf"$\phi$ = {float(f.stem[3:]):.1f}", fontsize=8)
            axs[r, 0].text(-0.15, 0.5, f"{round(100 * x)}% Cr", transform=axs[r, 0].transAxes, ha="right",
                           va="center", fontsize=8)
        fig.suptitle("MD spinodal decomposition, whole 25$^3$ cells (red Fe, blue Cr)", fontsize=9, y=1.0)
        fig.savefig(OUT / "cube_montage.png", dpi=200, bbox_inches="tight")
        plt.close(fig)
    print("compositions:", list(allser))


if __name__ == "__main__":
    main()
