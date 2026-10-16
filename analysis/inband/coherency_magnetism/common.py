#!/usr/bin/env python3
"""Shared helpers for the coherency / alpha'-magnetism analysis (analysis/inband/coherency_magnetism).

Only reads project data; writes only inside this folder.
"""
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCRIPTS = ROOT / "scripts"
RUNS = ROOT / "dft_qe" / "engaging" / "runs"
ANGLES = ROOT / "data" / "saw" / "angles_10k.csv"
OUT = HERE / "saw"
OUT.mkdir(exist_ok=True)
sys.path.insert(0, str(SCRIPTS))
import saw_du as sd  # noqa: E402

RHO_FE, RHO_CR = 7870.0, 7200.0
LAM = 6.4e-6
BAND = 609.0          # MHz, max (P)SAW of pristine CASS
XA, XB = 0.14, 0.88   # alpha, alpha' compositions
A_FE, A_CR = 2.8665, 2.8846


def rho(x):
    return x * RHO_CR + (1 - x) * RHO_FE


def a_vegard(x):
    return A_FE + x * (A_CR - A_FE)


def qe_aexp():
    """QE fixed-experimental-lattice SQS set: {x: (C11, C12, C44)} (+ x=0.80 crafm cell)."""
    rows = {}
    for d in sorted(RUNS.glob("sqs40_x*_k442_e60_aexp*")):
        if "INVALID" in d.name:
            continue
        j = d / "result.json"
        if j.exists():
            r = json.loads(j.read_text())
            rows[round(r["x_cr"], 3)] = (r["C11"], r["C12"], r["C44"])
    return rows


def qe_elem(name):
    r = json.loads((RUNS / name / "result.json").read_text())
    return np.array([r["C11"], r["C12"], r["C44"]])


def interp(rows, x, xs=None):
    xs = sorted(rows) if xs is None else xs
    C = np.array([rows[v] for v in xs])
    return np.array([np.interp(x, xs, C[:, i]) for i in range(3)])


def highest(path):
    """Highest (P)SAW per grain in MHz (PSAW kept if i_psaw/i_saw >= 2e-3)."""
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    keep = (d[:, 3] / d[:, 2] >= 2e-3) & np.isfinite(d[:, 1])
    v = np.where(keep, np.fmax(d[:, 0], d[:, 1]), d[:, 0])
    return v / LAM / 1e6


def summary(path):
    f = highest(path)
    return dict(fmax=np.nanmax(f), p999=np.nanpercentile(f, 99.9), p99=np.nanpercentile(f, 99),
                frac_in=float(np.mean(f > BAND)))


def run10k(tag, c11, c12, c44, rho_, workers=6):
    out = OUT / f"{tag}.csv"
    if not out.exists():
        subprocess.run([sys.executable, str(SCRIPTS / "saw_random_orientations.py"), "run", f"{c11}", f"{c12}",
                        f"{c44}", f"{rho_}", str(ANGLES), str(out), "10000", str(workers)], check=True)
    return out


def grain_fmax(C_ijkl, rho_, euler, sampling=4000, filt=True):
    """Highest mode (MHz) for one grain with a general C_ijkl (Pa). filt: apply the PSAW intensity rule (>= 2e-3);
    filt=False returns max(SAW, PSAW) regardless of intensity (the PSAW-branch velocity)."""
    v, inten = sd.get_saw(C_ijkl, rho_, euler, 0, sampling=sampling, return_intensity=True)
    keep = np.isfinite(v[1]) and (not filt or inten[1] / inten[0] >= 2e-3)
    return (max(v[0], v[1]) if keep else v[0]) / LAM / 1e6


def grain_fmax_dither(C_ijkl, rho_, euler, n=8, step=2.2e-3, filt=False):
    """Sub-grid estimate of the highest kept mode (MHz) at the standard sampling (4000).

    The solver's slowness grid quantises v in steps of ~0.1 %. Since v(C, s*rho) = v(C, rho)/sqrt(s), evaluating
    at n density scalings spread over one grid step and rescaling dithers the quantisation; the mean is then
    precise to ~0.05 MHz (needed for the sub-MHz coherency shifts). Default filt=False: the PSAW intensity ratio of
    the top grains sits at the 2e-3 threshold and flickers under the dither, so the branch velocity is tracked."""
    out = []
    for k in range(n):
        s = 1 + step * k / n
        out.append(grain_fmax(C_ijkl, rho_ * s, euler, 4000, filt) * np.sqrt(s))
    return float(np.mean(out))
