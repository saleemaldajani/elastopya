#!/usr/bin/env python3
"""Task 3 (tensors): effective tensors of decomposed Fe-36Cr (alpha 0.14 / alpha' 0.88, CALPHAD 748 K) from the
QE (experimental lattice) and EMTO-CPA random-alloy tensors, for every homogenization scheme.
Writes fecr36_tensors.csv (all schemes, cubic constants or full tetragonal constants, isotropic Rayleigh f) and
saw_runs/<model>_<scheme>.npy (6x6 Voigt, GPa) for the SAW runs.
"""
import csv
import json

import numpy as np

import data_io as io
import geometry as geo
import homog as h

X0, XA, XB = 0.36, 0.14, 0.88
W = (X0 - XA) / (XB - XA)


def models():
    xq, Cq = io.qe_aexp()
    xe, Ce = io.emto()
    return {
        "qe": dict(Ca=h.cubic(*io.interp(xq, Cq, XA)), Cb=h.cubic(*io.interp(xq, Cq, XB)),
                   Crand=h.cubic(*Cq[np.isclose(xq, 0.35)][0]), x_rand=0.35,
                   note="random = SQS x=0.35 (rho at 0.35); alpha/alpha' linear interpolation of the aexp QE set"),
        "emto": dict(Ca=h.cubic(*io.interp(xe, Ce, XA)), Cb=h.cubic(*io.interp(xe, Ce, XB)),
                     Crand=h.cubic(*io.interp(xe, Ce, X0)), x_rand=X0,
                     note="random = linear interpolation at 0.36; alpha/alpha' interpolated"),
    }


def fft_geom(Ca, Cb, kind, seeds=(0, 1, 2), n=45, **kw):
    out = []
    for s in seeds:
        m = geo.spinodal_mask(n, W, kind=kind, seed=s, **kw)
        fields = [np.where(m, h.cubic_constants(Cb)[i], h.cubic_constants(Ca)[i]) for i in range(3)]
        out.append(h.fft_homogenize(*fields, tol=1e-6))
    Cm = np.mean(out, 0)
    return Cm, np.std([h.cubic_constants(c) for c in out], 0, ddof=1)


def all_schemes(Ca, Cb, w, with_fft=True):
    L = h.laminate(Ca, Cb, w, axis=2)
    hsu, hsl = h.hs_bounds_cubic(Ca, Cb, w)
    h1u, h1l = h.hs_bounds_cubic(Ca, Cb, w, pfun=h.p_100)
    S = {
        "voigt": h.voigt(Ca, Cb, w), "reuss": h.reuss(Ca, Cb, w), "hill": h.hill(Ca, Cb, w),
        "laminate_tet": L, "laminate_V3": h.cubic_project(L), "laminate_R3": h.cubic_project_reuss(L),
        "hs_sphere_upper": hsu, "hs_sphere_lower": hsl, "sc_sphere": h.self_consistent(Ca, Cb, w),
        "hs100_upper": h1u, "hs100_lower": h1l, "sc_100": h.self_consistent(Ca, Cb, w, pfun=h.p_100),
    }
    sd = {}
    if with_fft:
        for kind, kw in (("axes", {}), ("cone", dict(spread_deg=15.0)), ("iso", {})):
            C, s = fft_geom(Ca, Cb, kind, **kw)
            S[f"fft_{kind}_raw"] = C
            S[f"fft_{kind}"] = h.cubic_project(C)
            sd[f"fft_{kind}"] = s
    return S, sd


def main():
    rows = []
    meta = {}
    for name, m in models().items():
        Ca, Cb = m["Ca"], m["Cb"]
        meta[name] = dict(alpha=h.cubic_constants(Ca), alpha_prime=h.cubic_constants(Cb), w=W, note=m["note"],
                          A_alpha=h.zener(*h.cubic_constants(Ca)), A_alpha_prime=h.zener(*h.cubic_constants(Cb)))
        S, sd = all_schemes(Ca, Cb, W)
        entries = [("random", m["Crand"], m["x_rand"]), ("alpha", Ca, XA), ("alpha_prime", Cb, XB)] + \
                  [(k, v, X0) for k, v in S.items()]
        for sch, C, x in entries:
            c11, c12, c44 = h.cubic_constants(C)
            r = dict(model=name, scheme=sch, rho=h.rho_of(x), C11=C[0, 0], C22=C[1, 1], C33=C[2, 2], C12=C[0, 1],
                     C13=C[0, 2], C23=C[1, 2], C44=C[3, 3], C55=C[4, 4], C66=C[5, 5],
                     c11_cub=c11, c12_cub=c12, c44_cub=c44, A_cub=h.zener(c11, c12, c44),
                     A_44=2 * C[3, 3] / (C[0, 0] - C[0, 1]), A_66=2 * C[5, 5] / (C[0, 0] - C[0, 1]),
                     f_iso=h.rayleigh_f(C, h.rho_of(x)),
                     max_offcubic=float(np.max(np.abs(C - h.cubic_project(C)))))
            if sch in sd:
                r["sd_c11"], r["sd_c12"], r["sd_c44"] = sd[sch]
            rows.append(r)
            if not sch.endswith("_raw"):
                np.save(io.HERE / "saw_runs" / f"{name}_{sch}.npy", C)
    keys = list(rows[0].keys()) + ["sd_c11", "sd_c12", "sd_c44"]
    with open(io.HERE / "fecr36_tensors.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, keys); wr.writeheader()
        for r in rows:
            wr.writerow({k: (f"{r[k]:.4f}" if isinstance(r.get(k), (float, np.floating)) else r.get(k, "")) for k in keys})
    (io.HERE / "fecr36_meta.json").write_text(json.dumps(meta, indent=1, default=float))
    for r in rows:
        print(f"{r['model']:5s} {r['scheme']:16s} C11 {r['C11']:6.1f} C33 {r['C33']:6.1f} C12 {r['C12']:6.1f} C13 {r['C13']:6.1f} "
              f"C44 {r['C44']:6.1f} C66 {r['C66']:6.1f} | cub {r['c11_cub']:6.1f} {r['c12_cub']:6.1f} {r['c44_cub']:6.1f} "
              f"A {r['A_cub']:.3f} f_iso {r['f_iso']:6.1f} offcub {r['max_offcubic']:.2f}")
    print(json.dumps(meta, indent=1, default=float))


if __name__ == "__main__":
    main()
