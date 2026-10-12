#!/usr/bin/env python3
"""Collect elastic constants from lammps/decomposition/runs/x*_L*_a*/stiff_phi*/log.lammps.

Writes data/md_decomposition.csv: one row per (x_cr, cell size L, alloy, phi), with the
cubic averages printed by stiffness.in and the relaxed box volume.
"""

import csv
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "lammps" / "decomposition" / "runs"
OUT = ROOT / "data" / "md_decomposition.csv"


def main():
    rows = []
    for log in sorted(RUNS.glob("x*_L*_a*/stiff_phi*/log.lammps")):
        m = re.match(r"x([\d.]+)_L(\d+)_a(\d+)", log.parents[1].name)
        phi = float(log.parent.name.replace("stiff_phi", ""))
        text = log.read_text()
        c = {k: float(v) for k, v in re.findall(r"^Elastic Constant (C\d\d)all = ([-\d.eE+]+) GPa", text, re.M)}
        if len(c) < 9 or not m:
            continue
        box = re.search(r"^RELAXED_BOX ([\d.eE+-]+) ([\d.eE+-]+) ([\d.eE+-]+)", text, re.M)
        rows.append(dict(
            x_cr=float(m.group(1)), L=int(m.group(2)), alloy=int(m.group(3)), phi=phi,
            c11=(c["C11"] + c["C22"] + c["C33"]) / 3,
            c12=(c["C12"] + c["C13"] + c["C23"]) / 3,
            c44=(c["C44"] + c["C55"] + c["C66"]) / 3,
            volume_A3=(float(box.group(1)) * float(box.group(2)) * float(box.group(3))) if box else "",
        ))
    if not rows:
        raise SystemExit("no completed decomposition runs")
    with open(OUT, "w", newline="") as fh:
        w = csv.DictWriter(fh, list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    for r in rows:
        print(f"x={r['x_cr']} L={r['L']} phi={r['phi']:.1f} C11={r['c11']:.2f}")


if __name__ == "__main__":
    main()
