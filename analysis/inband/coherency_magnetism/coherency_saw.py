#!/usr/bin/env python3
"""Explicit (P)SAW check of the coherency estimate on the top grains of dec36_AF (QE alpha + QE alpha' AF).

For every TOEC combination of coherency.sweep (Fe source x Cr source x deviatoric splits), the decomposed-medium
change is computed; the extreme cases (largest |change| of C44-type constants, cubic 3-variant average and single
tetragonal variant) are evaluated on the top-K grains (sub-grid dithering, PSAW branch tracked) at 1x and 10x the
coherency strain. Writes coherency_table.csv and coherency_saw.csv.
Usage: coherency_saw.py [eta]   (default 0.0039; Vegard misfit 0.0047 as sensitivity)
"""
import sys

import numpy as np

import coherency as co
import common as cm
import topgrains as tg


def main():
    eta = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0039
    T = {r[0]: np.array([float(v) for v in r[2:6]]) for r in
         (l.split(",") for l in (cm.HERE / "tensors_magnetic.csv").read_text().split("\n")[1:] if l)}
    Ca, Cb = T["alpha_qe"][:3], T["aprime_AF"][:3]
    rows = co.sweep(Ca, Cb, eta)
    tab = ["fe_source,cr_source,d3_fe,d3_cr,d1_cr,"
           "a_dC11,a_dC33,a_dC12,a_dC13,a_dC44,a_dC66,ap_dC11,ap_dC33,ap_dC12,ap_dC13,ap_dC44,ap_dC66,"
           "dec_cub_dC11,dec_cub_dC12,dec_cub_dC44,dec_tet_dC11,dec_tet_dC33,dec_tet_dC12,dec_tet_dC13,"
           "dec_tet_dC44,dec_tet_dC66"]
    D = []
    for r in rows:
        d = co.describe(r, Ca, Cb)
        D.append(d)
        v = [d["alpha"][k] for k in ("dC11", "dC33", "dC12", "dC13", "dC44", "dC66")]
        v += [d["alpha'"][k] for k in ("dC11", "dC33", "dC12", "dC13", "dC44", "dC66")]
        v += list(d["dec_cubic"]) + list(d["dec_tetra"].values())
        tab.append(f"{r['fe']},{r['cr']},{r['d3a']:.0f},{r['d3b']:.0f},{r['d1b']:.0f}," + ",".join(f"{x:.3f}" for x in v))
    (cm.HERE / f"coherency_table_eta{eta:g}.csv").write_text("\n".join(tab) + "\n")

    cub44 = np.array([d["dec_cubic"][2] for d in D])
    tet = np.array([[d["dec_tetra"]["dC44"], d["dec_tetra"]["dC66"]] for d in D])
    pick = {"cubic, most negative dC44": int(np.argmin(cub44)), "cubic, most positive dC44": int(np.argmax(cub44)),
            "tetra, max dC44": int(np.argmax(tet[:, 0])), "tetra, max dC66": int(np.argmax(tet[:, 1])),
            "tetra, min dC66": int(np.argmin(tet[:, 1]))}
    ids = tg.top_ids(cm.OUT / "dec36_AF.csv", 40)
    a0, r0 = rows[0]["a0"], rows[0]["r0"]
    f0 = tg.eval_top(ids, a0 * 1e9, r0)
    out = [f"case,amp,fmax_top_MHz,dfmax_MHz,dfmax_per_unit_strain_MHz,base_fmax_top_MHz"]
    print(f"base top-{len(ids)} max {f0.max():.2f} MHz", flush=True)
    for name, i in pick.items():
        r = rows[i]
        for amp in (1, 10):
            da = r["a_dec"] - a0
            if name.startswith("cubic"):
                da = co.rot_variants(da)
            a = a0 + amp * da
            rho = r0 + amp * (r["r_dec"] - r0)
            f = tg.eval_top(ids, a * 1e9, rho)
            out.append(f"{name},{amp},{f.max():.3f},{f.max() - f0.max():.3f},{(f.max() - f0.max()) / amp:.3f},"
                       f"{f0.max():.3f}")
            print(out[-1], "| Fe", r["fe"], "Cr", r["cr"], "d3", r["d3a"], r["d3b"], flush=True)
    (cm.HERE / f"coherency_saw_eta{eta:g}.csv").write_text("\n".join(out) + "\n")


if __name__ == "__main__":
    main()
