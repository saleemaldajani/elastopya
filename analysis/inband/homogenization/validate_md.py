#!/usr/bin/env python3
"""Task 2: validate homogenization schemes against the explicit MD (EAM) decomposed cells.

For each MD composition x: phase tensors = MD random-alloy tensors interpolated at the MD phase compositions
(x_alpha, x_alpha'; data/md_phase_compositions.csv, phi = 1), w = (x - x_alpha)/(x_alpha' - x_alpha) (lever rule),
homogenize, take the cubic projection (MD reports (C11+C22+C33)/3 etc., which is the same projection), and compare
with the explicit MD phi = 1 tensor (mean of five alloy seeds).
FFT schemes use the actual MD microstructures (one voxel per bcc cell, 25^3), per seed, then averaged.
Writes md_validation.csv and md_validation_summary.txt.
"""
import csv
import sys

import numpy as np

import data_io as io
import homog as h

xs_r, C_r = io.md_random()
dec = io.md_decomposition()
phases = io.md_phases()
XS = sorted(phases)


def schemes(Ca, Cb, w):
    L = h.laminate(Ca, Cb, w, axis=2)
    hsu, hsl = h.hs_bounds_cubic(Ca, Cb, w)
    h1u, h1l = h.hs_bounds_cubic(Ca, Cb, w, pfun=h.p_100)
    return {
        "voigt": h.voigt(Ca, Cb, w), "reuss": h.reuss(Ca, Cb, w), "hill": h.hill(Ca, Cb, w),
        "lam_V3": h.cubic_project(L), "lam_R3": h.cubic_project(h.cubic_project_reuss(L)),
        "hs_upper": hsu, "hs_lower": hsl, "sc_sphere": h.self_consistent(Ca, Cb, w),
        "hs100_upper": h1u, "hs100_lower": h1l, "sc_100": h.self_consistent(Ca, Cb, w, pfun=h.p_100),
    }


def fft_md(x, seeds=range(1, 6), shells=3):
    """FFT homogenization on the real MD cells: (i) two-phase map thresholded at the Otsu threshold with the phase
    tensors; (ii) continuous map C(x_local) from the MD random-alloy curve."""
    ph = phases[x]
    Ca, Cb = io.interp(xs_r, C_r, ph["xa"]), io.interp(xs_r, C_r, ph["xb"])
    two, cont, fr = [], [], []
    for s in seeds:
        f = io.MDRUNS / f"x{x:g}_L25_a{s}" / "phi1.000.data"
        if not f.exists():
            continue
        g, _ = io.md_voxel_field(f, shells=shells)
        m = g > ph["thr"]
        fr.append(m.mean())
        fields = [np.where(m, Cb[i], Ca[i]) for i in range(3)]
        two.append(h.cubic_project(h.fft_homogenize(*fields, tol=1e-6)))
        Cl = np.stack([np.interp(g, xs_r, C_r[:, i]) for i in range(3)])
        cont.append(h.cubic_project(h.fft_homogenize(*Cl, tol=1e-6)))
    return np.mean(two, 0), np.mean(cont, 0), float(np.mean(fr)), len(fr)


def metrics(C, x):
    c11, c12, c44 = h.cubic_constants(C)
    return dict(c11=c11, c12=c12, c44=c44, A=h.zener(c11, c12, c44), f=h.rayleigh_f(h.cubic(c11, c12, c44), h.rho_of(x)))


def main():
    do_fft = "--nofft" not in sys.argv
    rows = []
    for x in XS:
        ph = phases[x]
        Ca, Cb = h.cubic(*io.interp(xs_r, C_r, ph["xa"])), h.cubic(*io.interp(xs_r, C_r, ph["xb"]))
        w = (x - ph["xa"]) / (ph["xb"] - ph["xa"])
        md1 = dec[(x, 1.0)]; md05 = dec[(x, 0.5)]
        ref = metrics(h.cubic(*md1["mean"]), x)
        # seed scatter of the MD reference
        fs = [h.rayleigh_f(h.cubic(*c), h.rho_of(x)) for c in md1["all"]]
        As = [h.zener(*c) for c in md1["all"]]
        base = dict(x=x, xa=ph["xa"], xb=ph["xb"], w=w)
        rows.append(dict(base, scheme="MD_phi1", **ref, sd_c11=md1["sd"][0], sd_c12=md1["sd"][1], sd_c44=md1["sd"][2],
                         sd_A=np.std(As, ddof=1), sd_f=np.std(fs, ddof=1)))
        rows.append(dict(base, scheme="MD_phi0.5", **metrics(h.cubic(*md05["mean"]), x)))
        rows.append(dict(base, scheme="MDrandom_interp", **metrics(h.cubic(*io.interp(xs_r, C_r, x)), x)))
        for name, C in schemes(Ca, Cb, w).items():
            rows.append(dict(base, scheme=name, **metrics(C, x)))
        if do_fft:
            two, cont, fr, n = fft_md(x)
            rows.append(dict(base, scheme="fft_md_2phase", **metrics(two, x), w_fft=fr))
            rows.append(dict(base, scheme="fft_md_continuous", **metrics(cont, x)))
            print(f"x={x}: FFT on {n} seeds, voxel alpha' fraction {fr:.3f} (lever {w:.3f})", flush=True)
    keys = ["x", "scheme", "xa", "xb", "w", "c11", "c12", "c44", "A", "f", "sd_c11", "sd_c12", "sd_c44", "sd_A", "sd_f", "w_fft"]
    with open(io.HERE / "md_validation.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, keys); wr.writeheader()
        for r in rows:
            wr.writerow({k: (f"{r[k]:.5g}" if isinstance(r.get(k), (float, np.floating)) else r.get(k, "")) for k in keys})
    summarize(rows)


def summarize(rows):
    names = [r["scheme"] for r in rows if r["x"] == XS[0]][1:]
    out = []
    out.append("RMS / max |error| vs MD phi=1 over all MD compositions (and over the miscibility-gap subset 0.2<=x<=0.8)")
    out.append(f"{'scheme':20s} {'dC11':>6s} {'dC12':>6s} {'dC44':>6s} {'dA':>7s} {'df_MHz':>7s} | {'df gap':>7s} {'max df':>7s} {'dA gap':>7s}  shift_f%")
    ref = {r["x"]: r for r in rows if r["scheme"] == "MD_phi1"}
    rnd = {r["x"]: r for r in rows if r["scheme"] == "MD_phi0.5"}
    for nm in names:
        sel = [r for r in rows if r["scheme"] == nm]
        d = {k: np.array([r[k] - ref[r["x"]][k] for r in sel]) for k in ("c11", "c12", "c44", "A", "f")}
        gap = np.array([0.2 <= r["x"] <= 0.8 for r in sel])
        rms = {k: np.sqrt(np.mean(v ** 2)) for k, v in d.items()}
        # fraction of the MD random->decomposed frequency shift reproduced, gap compositions
        sh = np.array([(r["f"] - rnd[r["x"]]["f"]) / (ref[r["x"]]["f"] - rnd[r["x"]]["f"]) for r in sel])[gap]
        out.append(f"{nm:20s} {rms['c11']:6.2f} {rms['c12']:6.2f} {rms['c44']:6.2f} {rms['A']:7.3f} {rms['f']:7.2f} | "
                   f"{np.sqrt(np.mean(d['f'][gap]**2)):7.2f} {np.max(np.abs(d['f'])):7.2f} {np.sqrt(np.mean(d['A'][gap]**2)):7.3f}  "
                   f"{100*np.median(sh):5.0f}")
    sdf = np.array([ref[x]["sd_f"] for x in XS]); sdA = np.array([ref[x]["sd_A"] for x in XS])
    out.append(f"MD seed s.d. (phi=1): f median {np.median(sdf):.2f} MHz, A median {np.median(sdA):.4f}")
    out.append("shift_f% = median over 0.2<=x<=0.8 of (f_scheme - f_MD,phi0.5)/(f_MD,phi1 - f_MD,phi0.5)")
    out.append("")
    out.append("Fe-36Cr detail (x=0.36):")
    for r in rows:
        if r["x"] == 0.36:
            out.append(f"  {r['scheme']:20s} C11 {r['c11']:7.2f} C12 {r['c12']:7.2f} C44 {r['c44']:7.2f} A {r['A']:.3f} f {r['f']:6.1f} MHz")
    txt = "\n".join(out)
    print(txt)
    (io.HERE / "md_validation_summary.txt").write_text(txt + "\n")


if __name__ == "__main__":
    main()
