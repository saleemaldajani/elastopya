#!/usr/bin/env python3
"""Collect runs/*.csv into results.csv: per-grain highest (P)SAW frequency at 6.4 um for several PSAW intensity cuts.
f_shear = sqrt(C44/rho)/Lambda is listed because the band top with the 2e-3 cut sits at it."""
import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LAM = 6.4e-6
T = {r["tag"]: r for r in csv.DictReader(open(HERE / "tensors.csv"))}
T["ss304"] = dict(c11=209, c12=133, c44=121, rho=7950, A=2 * 121 / 76, note="Acta Fig. 6 austenite tensor")
T["x036_random_eq_voigt"] = dict(c11=274.096, c12=110.616, c44=110.784, rho=7628.8, A=1.355,
                                 note="Vegard random 0.36 (= Vegard Voigt decomposed); Fe C12 134.7")
out = []
for p in sorted((HERE / "runs").glob("*.csv")):
    if p.stem.endswith(".tmp"):
        continue
    d = np.genfromtxt(p, delimiter=",", names=True)
    r = d["i_psaw"] / d["i_saw"]; f = d["psaw_m_s"] / LAM / 1e6; s = d["saw_m_s"] / LAM / 1e6
    t = T[p.stem]; row = dict(tag=p.stem, c11=t["c11"], c12=t["c12"], c44=t["c44"], rho=t["rho"], A=round(float(t["A"]), 3),
                              f_shear=round(np.sqrt(float(t["c44"]) * 1e9 / float(t["rho"])) / LAM / 1e6, 1),
                              saw_min=round(np.nanmin(s), 1), saw_max=round(np.nanmax(s), 1))
    for cut in (2e-3, 5e-4, 0.0):
        top = np.where(np.isfinite(f) & (r >= cut), np.fmax(s, f), s)
        k = f"{cut:g}"
        row[f"fmed_{k}"] = round(np.median(top), 1); row[f"fmax_{k}"] = round(top.max(), 1)
        row[f"pct609_{k}"] = round(100 * np.mean(top > 609), 2); row[f"pct627_{k}"] = round(100 * np.mean(top > 627), 2)
    out.append(row)
with open(HERE / "results.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(out[0])); w.writeheader(); w.writerows(out)
for row in out:
    print(f"{row['tag']:24s} A={row['A']:<6} fsh={row['f_shear']:6.1f} SAW {row['saw_min']:5.0f}-{row['saw_max']:4.0f} | "
          f"2e-3: med {row['fmed_0.002']:4.0f} max {row['fmax_0.002']:4.0f} >609 {row['pct609_0.002']:5.2f}% | "
          f"5e-4: max {row['fmax_0.0005']:4.0f} >609 {row['pct609_0.0005']:5.2f}% | 0: >609 {row['pct609_0']:5.2f}%")
