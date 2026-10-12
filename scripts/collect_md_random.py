#!/usr/bin/env python3
"""Collect elastic constants from lammps/random_alloys/runs/*/log.lammps.

Writes data/md_random_alloys_runs.csv (one row per alloy) and data/md_random_alloys.csv
(mean, SD, and SE per composition). Cubic constants are the averages printed by
stiffness.in: C11 = (C11+C22+C33)/3, C12 = (C12+C13+C23)/3, C44 = (C44+C55+C66)/3.
"""

import csv
import re
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "lammps" / "random_alloys" / "runs"
DATA = ROOT / "data"
AMU = 1.66053906660e-27  # kg
M_FE, M_CR = 55.85, 51.996  # masses used in init.mod


def parse(log):
    text = log.read_text()
    c = {k: float(v) for k, v in re.findall(r"^Elastic Constant (C\d\d)all = ([-\d.eE+]+) GPa", text, re.M)}
    box = re.search(r"^RELAXED_BOX ([\d.eE+-]+) ([\d.eE+-]+) ([\d.eE+-]+)", text, re.M)
    ncr = re.search(r"^\s*(\d+) settings made for type/fraction", text, re.M)
    natoms = re.search(r"Created (\d+) atoms", text)
    if len(c) < 9 or not box or not natoms:
        return None
    n = int(natoms.group(1))
    n_cr = int(ncr.group(1)) if ncr else 0
    vol = np.prod([float(v) for v in box.groups()]) * 1e-30
    rho = ((n - n_cr) * M_FE + n_cr * M_CR) * AMU / vol
    return dict(
        x_cr_actual=n_cr / n,
        c11=(c["C11"] + c["C22"] + c["C33"]) / 3,
        c12=(c["C12"] + c["C13"] + c["C23"]) / 3,
        c44=(c["C44"] + c["C55"] + c["C66"]) / 3,
        rho_relaxed=rho,
    )


def main():
    rows = []
    for d in sorted(RUNS.glob("x*_s*")):
        m = re.match(r"x([\d.]+)_s(\d+)", d.name)
        r = parse(d / "log.lammps") if (d / "log.lammps").exists() else None
        if m and r:
            rows.append(dict(x_cr=float(m.group(1)), seed=int(m.group(2)), **r))
    if not rows:
        raise SystemExit("no completed runs")
    keys = list(rows[0])
    with open(DATA / "md_random_alloys_runs.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, keys)
        w.writeheader()
        w.writerows(rows)
    out = []
    for x in sorted({r["x_cr"] for r in rows}):
        g = [r for r in rows if r["x_cr"] == x]
        rec = dict(x_cr=x, n=len(g))
        for k in ("c11", "c12", "c44", "rho_relaxed", "x_cr_actual"):
            v = np.array([r[k] for r in g])
            rec[k] = v.mean()
            rec[k + "_sd"] = v.std(ddof=1) if len(v) > 1 else 0.0
        out.append(rec)
    with open(DATA / "md_random_alloys.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, list(out[0]))
        w.writeheader()
        w.writerows(out)
    for r in out:
        print({k: round(v, 4) if isinstance(v, float) else v for k, v in r.items()})


if __name__ == "__main__":
    main()
