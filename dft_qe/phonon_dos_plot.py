#!/usr/bin/env python3
"""Phonon DOS and vibrational thermodynamics of bcc Fe-Cr, compared with Fultz et al., PRB 52, 3280 (1995).

  phonon_dos_plot.py <outfile> [--root dft_qe/engaging]

(a) Phonon DOS of Fe, Fe70Cr30, Fe53Cr47 (our x = 0.475), Fe30Cr70 and Cr, stacked as in Fultz Fig. 1:
    40-atom SQS by finite displacements (phonopy/<x>/phonons.npz), with the DFPT DOS of the 1-atom cells
    (phonon/<case>/dos.dat) dashed for the pure metals.
(b) Vibrational entropy of mixing dS_vib = S(alloy) - (1-x) S(Fe) - x S(Cr) at 150 K, 300 K and T -> infinity
    (classical limit, 3 k_B [(1-x)<ln nu>_Fe + x <ln nu>_Cr - <ln nu>_alloy]), against Fultz Table II.
(c) Free energy of mixing vs T: dG = dH_static + dF_vib(T) - T S_conf (relaxed SQS energies; FM Fe and NM Cr
    references, as for the phonons), with dF_vib = 0 for comparison (configurational entropy only).
Panels without data yet are left empty with a note.
"""
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
FULTZ = {0.30: (0.098, 0.128, 0.141), 0.47: (0.142, 0.184, 0.201), 0.70: (0.148, 0.194, 0.214)}  # Table II, k_B/atom
ALLOYS = {0.30: "x0.30", 0.475: "x0.475", 0.70: "x0.70"}
CM1_THZ = 0.0299792458
KB = 8.617333262e-5


def npz(root, name):
    p = root / "phonopy" / name / "phonons.npz"
    return dict(np.load(p)) if p.exists() else None


def ln_nu(d):
    f, g = d["freq_thz"], d["dos_per_atom"]
    m = f > 0.05
    return float(np.trapezoid(g[m] * np.log(f[m]), f[m]) / np.trapezoid(g[m], f[m]))


def e_last(out):
    return float(re.findall(r"^!\s+total energy\s+=\s+([-\d.]+)", Path(out).read_text(), re.M)[-1])


def main(outfile, root):
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(HERE.parent / "scripts" / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    root = Path(root)
    fe, cr = npz(root, "Fe_ref"), npz(root, "Cr_ref")
    al = {x: npz(root, n) for x, n in ALLOYS.items()}
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(7.4, 2.9), gridspec_kw={"width_ratios": [1.1, 1, 1]})

    # (a) stacked DOS
    rows = [(0.0, "Fe", fe, "phonon/Fe_fm/dos.dat")] + [(x, f"Fe$_{{{round(100*(1-x))}}}$Cr$_{{{round(100*x)}}}$", al[x], None)
                                                    for x in ALLOYS] + [(1.0, "Cr", cr, "phonon/Cr_nm_v5_k24/dos.dat")]
    for k, (x, lab, d, dfpt) in enumerate(rows):
        off = 0.5 * k
        if d is not None:
            a.plot(d["freq_thz"], d["dos_per_atom"] / 3 + off, color="0.1", lw=0.9)   # unit area, as Fultz Fig. 1
        if dfpt and (root / dfpt).exists():
            q = np.loadtxt(root / dfpt)
            f, g = q[:, 0] * CM1_THZ, q[:, 1]
            g = g / np.trapezoid(g, f)                         # unit area
            a.plot(f, g + off, "--", color="#0072B2", lw=0.8)
        a.text(0.3, off + 0.32, lab, fontsize=7)
    a.set(xlim=(0, 11), ylim=(0, 0.5 * len(rows) + 0.4), xlabel=r"$\nu$ (THz)", ylabel=r"$g(\nu)$ (THz$^{-1}$), offset")
    a.plot([], [], color="0.1", lw=0.9, label="SQS, finite displ."); a.plot([], [], "--", color="#0072B2", lw=0.8, label="DFPT, 1-atom")
    a.legend(loc="upper right", fontsize=6)

    # (b) dS_vib vs x
    for xf, vals in FULTZ.items():
        for v, mk, lab in zip(vals, ("v", "s", "o"), ("150 K", "300 K", r"$T\to\infty$")):
            b.plot([xf], [v], mk, mfc="white", mec="0.2", ms=5, label=f"Fultz, {lab}" if xf == 0.30 else None)
    if fe is not None and cr is not None:
        xs, s150, s300, sinf = [], [], [], []
        for x, d in al.items():
            if d is None:
                continue
            s = lambda dd, T: float(np.interp(T, dd["T"], dd["S_kB_per_atom"]))
            xs.append(x)
            s150.append(s(d, 150) - (1 - x) * s(fe, 150) - x * s(cr, 150))
            s300.append(s(d, 300) - (1 - x) * s(fe, 300) - x * s(cr, 300))
            sinf.append(3 * ((1 - x) * ln_nu(fe) + x * ln_nu(cr) - ln_nu(d)))
        for ys, mk in ((s150, "v"), (s300, "s"), (sinf, "o")):
            b.plot(xs, ys, mk, color="#D55E00", ms=5)
        b.plot([], [], "o", color="#D55E00", label="QE (this work)")
    else:
        b.text(0.5, 0.5, "pending: Fe, Cr reference cells", transform=b.transAxes, ha="center", fontsize=7)
    b.axhline(0, color="0.6", lw=0.5)
    b.set(xlim=(0, 1), xlabel=r"$x_{\rm Cr}$", ylabel=r"$\Delta S_{\rm vib}$ ($k_{\rm B}$/atom)")
    b.legend(loc="upper left", fontsize=6)

    # (c) dG_mix(T)
    runs = root / "runs"
    try:
        e_fe = e_last(runs / "Fe_fm_k8x8x10_e60" / "relax.out") / 2
        e_cr = e_last(runs / "Cr_nm_k16x16x16_e60" / "relax.out") / 2
    except (FileNotFoundError, IndexError):
        e_fe = e_cr = None
    T = np.linspace(0, 1200, 121)
    for x, col in zip(ALLOYS, ("#E69F00", "#D55E00", "#0072B2")):
        rel = runs / f"sqs40_x{x:.2f}_k442_e60" / "relax.out"
        if x == 0.475:
            rel = runs / "sqs40_x0.475_k442_e60" / "relax.out"
        if e_fe is None or not rel.exists() or "End final coordinates" not in rel.read_text():
            continue
        dh = (e_last(rel) / 40 - (1 - x) * e_fe - x * e_cr) * 13.605693123
        sconf = -KB * (x * np.log(x) + (1 - x) * np.log(1 - x))
        c.plot(T, dh - T * sconf, ":", color=col, lw=0.9)
        d = al[x]
        if d is not None and fe is not None and cr is not None:
            df = np.interp(T, d["T"], d["F_eV_per_atom"]) - (1 - x) * np.interp(T, fe["T"], fe["F_eV_per_atom"]) \
                 - x * np.interp(T, cr["T"], cr["F_eV_per_atom"])
            c.plot(T, dh + df - T * sconf, "-", color=col, lw=1.1, label=f"{round(100*x)}% Cr")
        else:
            c.plot([], [], ":", color=col, label=f"{round(100*x)}% Cr")
    c.axhline(0, color="0.6", lw=0.5)
    c.plot([], [], "-", color="0.3", label="with vibrations"); c.plot([], [], ":", color="0.3", label="configurational only")
    c.set(xlim=(0, 1200), xlabel="T (K)", ylabel=r"$\Delta G_{\rm mix}$ (eV/atom)")
    c.legend(loc="lower left", fontsize=6)
    for ax, l in zip((a, b, c), "abc"):
        ax.text(-0.18, 1.02, f"({l})", transform=ax.transAxes, fontweight="bold", fontsize=9)
    fig.tight_layout(w_pad=0.8)
    fig.savefig(outfile, dpi=600); fig.savefig(Path(outfile).with_suffix(".pdf")); fig.savefig(Path(outfile).with_suffix(".svg"))
    print("wrote", outfile)


if __name__ == "__main__":
    args = sys.argv[1:]
    root = args[args.index("--root") + 1] if "--root" in args else str(HERE / "engaging")
    main(args[0], root)
