#!/usr/bin/env python3
"""Summary of coherency strains, stresses and tensor changes (ranges over all TOEC choices); no SAW."""
import sys
import numpy as np
import coherency as co, common as cm
T = {r[0]: np.array([float(v) for v in r[2:6]]) for r in
     (l.split(",") for l in (cm.HERE / "tensors_magnetic.csv").read_text().split("\n")[1:] if l)}
for eta in (0.0039, 0.0047):
    for apk in ("aprime_AF", "aprime_MEAS"):
        Ca, Cb = T["alpha_qe"][:3], T[apk][:3]
        rows = co.sweep(Ca, Cb, eta)
        r0 = rows[0]
        print(f"\n== eta={eta}, alpha'={apk}: strains alpha {np.diag(r0['ea']).round(5)}, alpha' {np.diag(r0['eb']).round(5)};"
              f" stress (GPa) alpha {np.diag(r0['sa']).round(3)}, alpha' {np.diag(r0['sb']).round(3)}")
        D = [co.describe(r, Ca, Cb) for r in rows]
        for ph in ("alpha", "alpha'"):
            for k in ("dC11", "dC33", "dC12", "dC13", "dC44", "dC66"):
                v = np.array([d[ph][k] for d in D]); print(f"  {ph:7s} {k}: {v.min():6.2f} .. {v.max():6.2f}")
            v = np.array([d[ph]["cub"] for d in D]); print(f"  {ph:7s} 3-variant cubic dC11,dC12,dC44: {v.min(0).round(2)} .. {v.max(0).round(2)}")
        v = np.array([d["dec_cubic"] for d in D]); print(f"  dec36 3-variant cubic dC11,dC12,dC44: {v.min(0).round(2)} .. {v.max(0).round(2)}")
        v = np.array([list(d["dec_tetra"].values()) for d in D]); print(f"  dec36 single variant (11,33,12,13,44,66): {v.min(0).round(2)} .. {v.max(0).round(2)}")
        print(f"  density change dec36: {np.mean([r['r_dec'] - r['r0'] for r in rows]):.3f} kg/m3")
