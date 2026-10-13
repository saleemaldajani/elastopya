#!/usr/bin/env python3
"""Spin-resolved Fermi-surface cross-sections of bcc Fe and Cr (Quantum ESPRESSO, PBE-PAW, 60 Ry).

  fermi_slices.py inputs <case> <dir>     write scf.in and bands.in for one case
  fermi_slices.py plot <outfile> <dir>... plot the k_z = 0 cross-sections from the runs' XML files
  fermi_slices.py korzhavyi <outfile> <dir>...  same data as a spectral function in the layout of
                                          Korzhavyi et al., PRB 79, 054202 (2009), Fig. 6

Cases
  Fe_fm   bcc Fe, 1-atom cell, ferromagnetic (a = 2.827 A, the converged QE lattice parameter)
  Cr_nm   bcc Cr, 1-atom cell, non-magnetic (a = 2.85 A)
  Cr_afm  Cr, 2-atom simple-cubic cell, commensurate antiferromagnet (a = 2.87 A); folded zone

The bands are computed on a uniform grid in the k_z = 0 plane (units 2 pi / a) and the Fermi
surface is drawn as the E = E_F contour of every band, majority and minority spin separately.
"""
import os
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

PSEUDO = os.environ.get("QE_PSEUDO", "/private/tmp/claude-501/-Users-saleem-MIT-research-Scripta/"
                        "37f8fe43-8e5f-4d11-883b-ebb44fa03560/scratchpad/pseudo")
HA_EV = 27.211386
N_GRID = 121
CASES = {
    "Fe_fm": dict(a=2.827, ibrav=3, nat=1, species=[("Fe", 55.845, "Fe.pbe-spn-kjpaw_psl.0.2.1.UPF")],
                  atoms=[("Fe", 0, 0, 0)], mag=["starting_magnetization(1) = 0.4"], nspin=2, kscf=16,
                  kmax=1.0, nbnd=16),
    "Cr_nm": dict(a=2.85, ibrav=3, nat=1, species=[("Cr", 51.996, "Cr.pbe-spn-kjpaw_psl.1.0.0.UPF")],
                  atoms=[("Cr", 0, 0, 0)], mag=[], nspin=1, kscf=16, kmax=1.0, nbnd=16),
    "Cr_afm": dict(a=2.87, ibrav=1, nat=2, species=[("Cr1", 51.996, "Cr.pbe-spn-kjpaw_psl.1.0.0.UPF"),
                                                    ("Cr2", 51.996, "Cr.pbe-spn-kjpaw_psl.1.0.0.UPF")],
                   atoms=[("Cr1", 0, 0, 0), ("Cr2", 0.5, 0.5, 0.5)],
                   mag=["starting_magnetization(1) = 0.5", "starting_magnetization(2) = -0.5"], nspin=2,
                   kscf=12, kmax=0.5, nbnd=28),
}
# Same cases at the room-temperature experimental lattice parameters (Fe 2.8665 A, Cr 2.8846 A), the volume
# convention of the EMTO-CPA literature; NM Cr with k = 24^3 (16^3 under-samples its nested Fermi surface in DFPT).
A_EXP = {"Fe_fm": 2.8665, "Cr_nm": 2.8846, "Cr_afm": 2.8846}
for _c, _a in A_EXP.items():
    CASES[_c + "_aexp"] = {**CASES[_c], "a": _a, **({"kscf": 24} if _c == "Cr_nm" else {})}


def base(case):
    """Fe_fm_aexp -> Fe_fm (structure, path and title do not depend on the lattice parameter)."""
    return case.removesuffix("_aexp")


def write_inputs(case, d):
    c = CASES[case]
    d = Path(d); d.mkdir(parents=True, exist_ok=True)
    bohr = c["a"] / 0.529177210903
    head = lambda calc, extra="": (
        f"&control\n  calculation = '{calc}', prefix = 'fs', outdir = './tmp', pseudo_dir = '{PSEUDO}'\n/\n"
        f"&system\n  ibrav = {c['ibrav']}, celldm(1) = {bohr:.6f}, nat = {c['nat']}, ntyp = {len(c['species'])},\n"
        f"  ecutwfc = 60, ecutrho = 480, occupations = 'smearing', smearing = 'mv', degauss = 0.02,\n"
        f"  nspin = {c['nspin']}, nbnd = {c['nbnd']}{extra}\n" + "".join(f"  {m}\n" for m in c["mag"]) + "/\n"
        "&electrons\n  conv_thr = 1e-9, mixing_beta = 0.3\n/\n"
        "ATOMIC_SPECIES\n" + "".join(f"{s} {m} {p}\n" for s, m, p in c["species"]) +
        "ATOMIC_POSITIONS alat\n" + "".join(f"{s} {x} {y} {z}\n" for s, x, y, z in c["atoms"]))
    k = c["kscf"]
    (d / "scf.in").write_text(head("scf") + f"K_POINTS automatic\n{k} {k} {k} 0 0 0\n")
    g = np.linspace(-c["kmax"], c["kmax"], N_GRID)
    kx, ky = np.meshgrid(g, g, indexing="ij")
    pts = "".join(f"{x:.6f} {y:.6f} 0.000000 1\n" for x, y in zip(kx.ravel(), ky.ravel()))
    (d / "bands.in").write_text(head("bands", ", nosym = .true.") + f"K_POINTS tpiba\n{kx.size}\n" + pts)


def read_xml(path):
    """k-points (2 pi/a), eigenvalues (eV) per spin [nspin, nk, nbnd], Fermi energy (eV), from QE XML."""
    root = ET.parse(path).getroot()
    bs = root.find("output/band_structure")
    lsda = bs.find("lsda").text.strip() == "true"
    ef = float(bs.find("fermi_energy").text) * HA_EV
    ks = bs.findall("ks_energies")
    alat = float(root.find("output/atomic_structure").get("alat"))
    k = np.array([[float(v) for v in e.find("k_point").text.split()] for e in ks])   # units 2 pi / alat
    ev = np.array([[float(v) for v in e.find("eigenvalues").text.split()] for e in ks]) * HA_EV
    if lsda:
        nb = ev.shape[1] // 2
        ev = np.stack([ev[:, :nb], ev[:, nb:]])
    else:
        ev = ev[None]
    return k, ev, ef


def plot(outfile, dirs):
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta")
    titles = {"Fe_fm": "bcc Fe (FM)", "Cr_nm": "bcc Cr (NM)", "Cr_afm": "Cr (AFM, folded zone)"}
    cols = {0: "#D55E00", 1: "#0072B2"}           # majority, minority
    fig, axs = plt.subplots(1, len(dirs), figsize=(2.5 * len(dirs), 2.6))
    axs = np.atleast_1d(axs)
    for ax, d in zip(axs, dirs):
        d = Path(d); case = d.name
        k, ev, ef = read_xml(d / "fs_data.xml" if (d / "fs_data.xml").exists() else d / "tmp" / "fs.save" / "data-file-schema.xml")
        n = int(round(np.sqrt(len(k))))
        kx, ky = k[:, 0].reshape(n, n), k[:, 1].reshape(n, n)
        for s in range(ev.shape[0]):
            for b in range(ev.shape[2]):
                e = ev[s, :, b].reshape(n, n) - ef
                if e.min() < 0 < e.max():
                    ax.contour(kx, ky, e, levels=[0], colors=cols[s], linewidths=1.0 if s == 0 else 0.9,
                               linestyles="-" if s == 0 else "--")
        km = CASES[case]["kmax"]
        if CASES[case]["ibrav"] == 3:   # bcc zone in the k_z = 0 plane: square with corners at H
            ax.plot([1, 0, -1, 0, 1], [0, 1, 0, -1, 0], color="0.5", lw=0.6)
            for lab, (x, y) in {"Γ": (0, 0), "H": (1, 0), "N": (0.5, 0.5)}.items():
                ax.text(x + 0.03, y + 0.03, lab, fontsize=7, color="0.3")
        ax.set(xlim=(-km, km), ylim=(-km, km), xlabel=r"$k_x$ ($2\pi/a$)", ylabel=r"$k_y$ ($2\pi/a$)")
        ax.set_aspect("equal")
        ax.set_title(titles.get(base(case), case) + (" (exp. a)" if case.endswith("_aexp") else ""), fontsize=8)
    axs[0].plot([], [], color=cols[0], label="majority"); axs[0].plot([], [], "--", color=cols[1], label="minority")
    axs[0].legend(loc="lower left", fontsize=7)
    fig.tight_layout()
    fig.savefig(outfile, dpi=600); fig.savefig(Path(outfile).with_suffix(".pdf"))
    print("wrote", outfile)


def spectral_map(k, ev, ef, eta=0.08, n=241):
    """Bloch spectral function at E_F, A(k) = sum_n (eta/pi) / ((E_n - E_F)^2 + eta^2) (eV), on the bcc
    (001) zone section drawn as Korzhavyi et al., PRB 79, 054202 (2009), Fig. 6: square with H at the
    corners, N at the edge midpoints, Gamma in the centre. Square coordinates (p, q) in [-1, 1] map to
    k = p (1/2, -1/2, 0) + q (1/2, 1/2, 0) in units of 2 pi / a."""
    from scipy.interpolate import RegularGridInterpolator
    m = int(round(np.sqrt(len(k))))
    g = k[:, 0].reshape(m, m)[:, 0]
    pq = np.linspace(-1, 1, n)
    P, Q = np.meshgrid(pq, pq, indexing="xy")
    KX, KY = (P + Q) / 2, (Q - P) / 2
    out = []
    for s in range(ev.shape[0]):
        A = np.zeros_like(P)
        for b in range(ev.shape[2]):
            e = ev[s, :, b].reshape(m, m) - ef
            if e.min() > 3 or e.max() < -3:
                continue
            f = RegularGridInterpolator((g, g), e)
            eb = f(np.stack([KX.ravel(), KY.ravel()], -1)).reshape(P.shape)
            A += (eta / np.pi) / (eb ** 2 + eta ** 2)
        out.append(A)
    return pq, out


def plot_korzhavyi(outfile, dirs):
    """Spin-up (top) and spin-down (bottom) Fermi-surface sections, one column per run directory."""
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta")
    labels = {"Fe_fm": "Fe (this work, QE)"}
    fig, axs = plt.subplots(2, len(dirs), figsize=(2.2 * len(dirs), 4.4), squeeze=False)
    for j, d in enumerate(dirs):
        d = Path(d)
        k, ev, ef = read_xml(d / "fs_data.xml" if (d / "fs_data.xml").exists() else d / "tmp" / "fs.save" / "data-file-schema.xml")
        pq, maps = spectral_map(k, ev, ef)
        for i, (A, spin) in enumerate(zip(maps, ("spin up", "spin down"))):
            ax = axs[i, j]
            ax.imshow(A, extent=(-1, 1, -1, 1), origin="lower", cmap="Greys", vmin=0, vmax=np.percentile(A, 99.5))
            ax.set_xticks([]); ax.set_yticks([]); ax.minorticks_off()
            for (x, y, t) in ((-1, -1, "H"), (1, -1, "H"), (-1, 1, "H"), (1, 1, "H"),
                              (0, -1, "N"), (0, 1, "N"), (-1, 0, "N"), (1, 0, "N")):
                ax.text(x * 1.08, y * 1.08, t, ha="center", va="center", fontsize=7, clip_on=False)
            ax.text(0.03, 0.03, "Γ", fontsize=7)
            ax.set_title(f"{labels.get(d.name, d.name)}, {spin}", fontsize=7, pad=14)
    fig.tight_layout(h_pad=1.5)
    fig.savefig(outfile, dpi=600); fig.savefig(Path(outfile).with_suffix(".pdf"))
    print("wrote", outfile)


if __name__ == "__main__":
    if sys.argv[1] == "inputs":
        write_inputs(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "korzhavyi":
        plot_korzhavyi(sys.argv[2], sys.argv[3:])
    elif sys.argv[1] == "plot":
        plot(sys.argv[2], sys.argv[3:])
