#!/usr/bin/env python3
"""Read-only loaders for the project data (nothing outside this folder is written)."""
import csv
import json
import re
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DATA = ROOT / "data"
RUNS = ROOT / "dft_qe" / "engaging" / "runs"
MDRUNS = ROOT / "lammps" / "decomposition" / "runs"


def read_csv(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def md_random():
    """MD (EAM) random alloys: x -> (C11, C12, C44), duplicate x averaged. Returns sorted x, (n,3) array."""
    d = read_csv(DATA / "md_random_alloys.csv")
    xs = np.unique(np.round(d["x_cr"], 6))
    C = np.array([[d[k][np.isclose(d["x_cr"], x)].mean() for k in ("c11", "c12", "c44")] for x in xs])
    return xs, C


def md_decomposition():
    """{(x, phi): dict(mean=(C11,C12,C44), sd=(...), n)} averaged over the five alloy seeds."""
    d = read_csv(DATA / "md_decomposition.csv")
    out = {}
    for x in np.unique(d["x_cr"]):
        for phi in np.unique(d["phi"]):
            m = np.isclose(d["x_cr"], x) & np.isclose(d["phi"], phi)
            if m.sum():
                arr = np.column_stack([d["c11"][m], d["c12"][m], d["c44"][m]])
                out[(round(float(x), 4), round(float(phi), 3))] = dict(
                    mean=arr.mean(0), sd=arr.std(0, ddof=1) if m.sum() > 1 else np.zeros(3), n=int(m.sum()), all=arr)
    return out


def md_phases():
    d = read_csv(DATA / "md_phase_compositions.csv")
    return {round(float(x), 4): dict(xa=xa, xb=xb, f=f, thr=t)
            for x, phi, xa, xb, f, t in zip(d["x_cr"], d["phi"], d["x_alpha"], d["x_alpha_prime"],
                                            d["f_alpha_prime"], d["threshold"]) if np.isclose(phi, 1.0)}


def qe_aexp(include_pure=True):
    rows = {}
    for d in sorted(RUNS.glob("sqs40_x*_k442_e60_aexp")) + [RUNS / "sqs40_x0.80_k442_e60_aexp_crafm"]:
        j = d / "result.json"
        if j.exists():
            r = json.loads(j.read_text())
            rows[round(r["x_cr"], 3)] = (r["C11"], r["C12"], r["C44"])
    if include_pure:
        for x, name in ((0.0, "Fe_fm_k16x16x16_e60_aexp"), (1.0, "Cr_afm_k16x16x16_e60_aexp")):
            r = json.loads((RUNS / name / "result.json").read_text())
            rows[x] = (r["C11"], r["C12"], r["C44"])
    x = np.array(sorted(rows))
    return x, np.array([rows[v] for v in x])


def emto():
    """EMTO-CPA (Razumovskiy et al. 2011, digitized): C11 = B + 4C'/3, C12 = B - 2C'/3."""
    d = read_csv(DATA / "sheet_raw_razumovskiy.csv")
    x, cp, c44, b = d["x_cr"], d["cprime_gpa"], d["c44_gpa"], d["bulk_gpa"]
    o = np.argsort(x)
    return x[o], np.column_stack([b + 4 * cp / 3, b - 2 * cp / 3, c44])[o]


def interp(xs, C, x):
    return np.array([np.interp(x, xs, C[:, i]) for i in range(3)])


# ---------------------------------------------------------------- MD microstructures
def read_lammps(path):
    txt = Path(path).read_text()
    lo_hi = [list(map(float, re.search(rf"([-\d.eE+]+) ([-\d.eE+]+) {a}lo {a}hi", txt).groups())) for a in "xyz"]
    body = txt.split("Atoms")[1].split("\n\n", 2)[1]
    arr = np.array([l.split()[:5] for l in body.strip().split("\n")], float)
    return arr[:, 1].astype(int), arr[:, 2:5], np.array(lo_hi)


def md_voxel_field(path, shells=3):
    """Cr fraction per bcc unit cell (one voxel per cell, n^3 grid), after smoothing every atom's occupancy over itself
    and its neighbours up to the `shells`-th bcc shell (shells=3: 27 atoms, as scripts/phase_compositions.py)."""
    from scipy.spatial import cKDTree
    t, p, box = read_lammps(path)
    Lb = box[:, 1] - box[:, 0]
    n = int(round(Lb[0] / 2.866)); a = Lb[0] / n
    q = (p - box[:, 0]) % Lb
    cr = (t == 2).astype(float)
    if shells == 0:
        loc = cr
    else:
        rcut = {1: 1.05 * a * np.sqrt(3) / 2, 2: 1.05 * a, 3: 1.05 * a * np.sqrt(2)}[shells]
        tree = cKDTree(q, boxsize=Lb)
        nb = tree.query_ball_point(q, r=rcut)
        loc = np.array([cr[i].mean() for i in nb])
    idx = np.floor(q / a + 0.25).astype(int) % n
    grid = np.zeros((n, n, n)); cnt = np.zeros((n, n, n))
    np.add.at(grid, tuple(idx.T), loc); np.add.at(cnt, tuple(idx.T), 1)
    assert cnt.min() == 2 and cnt.max() == 2, "expected two atoms per bcc cell"
    return grid / cnt, cr.mean()
