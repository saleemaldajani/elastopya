#!/usr/bin/env python3
"""Fermi-surface nesting of non-magnetic bcc Cr: the nesting function (Fermi-surface autocorrelation)
    xi(q) = (1/N) sum_k sum_{n,m} delta(e_nk - E_F) delta(e_m,k+q - E_F)
on a full n^3 grid of the primitive reciprocal cell (bands without symmetry), computed by FFT, along
Gamma-H ([001]) where the electron 'jack' at Gamma nests onto the hole octahedron at H, giving the SDW wavevector
Q = (2 pi/a)(1 - delta). Rigid-band filling with y extra Fe (2 more valence electrons per Fe) shows how Fe
detunes the nesting that drives the SDW of Cr (Fawcett, RMP 60, 209 (1988)).

  nesting.py inputs <dir> [n]     scf.in (k 24^3) and bands.in (full n^3 grid, nosym), Cr_nm of fermi_slices
  nesting.py analyse <dir>        xi(q) along Gamma-H for y = 0, 0.05, 0.10, 0.20 -> nesting.npz
"""
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fermi_slices as fs  # noqa: E402

SIGMA = 0.05          # eV, Gaussian width of the Fermi-surface delta functions
Y_FE = (0.0, 0.05, 0.10, 0.20)


def inputs(d, n=32):   # 32^3 = 32768 k-points: QE pw.x allows at most 40000 (npk); 40^3 aborted (2026-10-15)
    d = Path(d)
    fs.write_inputs("Cr_nm", d)
    scf = (d / "scf.in").read_text().replace("16 16 16 0 0 0", "24 24 24 0 0 0")
    (d / "scf.in").write_text(scf)
    head = scf.split("K_POINTS")[0].replace("calculation = 'scf'", "calculation = 'bands'")
    head = head.replace("nbnd = 16", "nbnd = 12, nosym = .true., noinv = .true.")
    # Davidson from random wavefunctions aborted in cdiaghg (MPI error 37) on the 32^3 nosym grid (2026-10-15). CG then
    # finished but left 47% of the k-points unconverged (3s level anywhere from -96 to -10 eV instead of -52.4 eV; the
    # grid E_F came out 1.6 eV above the SCF value) -- found 2026-10-16. ParO, as in the unfolded SQS bands, with full
    # accuracy for the empty states.
    head = head.replace("conv_thr = 1e-9, mixing_beta = 0.3",
                        "conv_thr = 1e-9, mixing_beta = 0.3, diagonalization = 'paro', diago_full_acc = .true.")
    g = np.arange(n) / n
    kk = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).reshape(-1, 3)
    (d / "bands.in").write_text(head + f"K_POINTS crystal\n{len(kk)}\n" +
                                "".join(f"{a:.6f} {b:.6f} {c:.6f} 1\n" for a, b, c in kk))
    (d / "grid.txt").write_text(f"{n}\n")
    print(d, "bands k-points:", len(kk))


def fermi_level(ev, nel, sigma=0.02):
    from scipy.special import erfc
    occ = lambda mu: 2 * np.mean(np.sum(0.5 * erfc((ev - mu) / (np.sqrt(2) * sigma)), axis=-1))
    lo, hi = ev.min(), ev.max()
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if occ(mid) < nel else (lo, mid)
    return 0.5 * (lo + hi)


def analyse(d):
    d = Path(d)
    n = int((d / "grid.txt").read_text())
    xml = d / "tmp" / "fs.save" / "data-file-schema.xml"
    if not xml.exists():
        xml = d / "bands_data.xml"
    _, ev, _ = fs.read_xml(xml)
    ev = ev[0]                                            # NM: one spin channel, [nk, nb]
    root = ET.parse(xml).getroot()
    B = np.array([[float(v) for v in root.find(f"output/basis_set/reciprocal_lattice/b{i}").text.split()]
                  for i in (1, 2, 3)])                    # rows b_i, units 2 pi/alat (alat = a)
    nel = 14.0                                            # Cr spn PAW: 3s2 3p6 3d5 4s1
    # guards (2026-10-16): every k-point must have the flat 3s level, and the grid E_F must match the SCF one
    if np.ptp(ev[:, 0]) > 1.0:
        sys.exit(f"bands not converged: lowest band spans {ev[:, 0].min():.1f} to {ev[:, 0].max():.1f} eV")
    ef_scf = float(root.find("output/band_structure/fermi_energy").text) * 27.211386
    if abs(fermi_level(ev, nel) - ef_scf) > 0.1:
        sys.exit(f"grid E_F {fermi_level(ev, nel):.3f} eV differs from the SCF {ef_scf:.3f} eV")
    grid = ev.reshape(n, n, n, -1)
    s = np.linspace(0, 1, 201)                            # Gamma -> H = (0, 0, 1) 2 pi/a
    cryst = np.outer(s, [0, 0, 1.0]) @ np.linalg.inv(B)   # crystal coordinates of q
    out = {"s": s}
    for y in Y_FE:
        ef = fermi_level(grid, nel + 2 * y)
        w = np.sum(np.exp(-0.5 * ((grid - ef) / SIGMA) ** 2), axis=-1) / (np.sqrt(2 * np.pi) * SIGMA)
        F = np.fft.fftn(w)
        xi = np.real(np.fft.ifftn(np.conj(F) * F)) / w.size ** 2      # autocorrelation, periodic on the grid
        idx = (cryst * n) % n                                          # trilinear interpolation
        i0 = np.floor(idx).astype(int); t = idx - i0
        val = 0.0
        for dx in (0, 1):
            for dy in (0, 1):
                for dz in (0, 1):
                    wgt = np.prod(np.where([dx, dy, dz], t, 1 - t), axis=1)
                    val = val + wgt * xi[(i0[:, 0] + dx) % n, (i0[:, 1] + dy) % n, (i0[:, 2] + dz) % n]
        out[f"xi_y{y:.2f}"] = val
        out[f"ef_y{y:.2f}"] = ef
        peak = s[50:][np.argmax(val[50:])]
        print(f"y_Fe = {y:.2f}: E_F shift {ef - out['ef_y0.00'] if y else 0:+.3f} eV, nesting peak at q = {peak:.3f} (2pi/a),"
              f" height {val[50:].max():.4g}")
    np.savez(d / "nesting.npz", **out)


def plot(npz, out):
    """SM Fig. S_nesting: xi(q) along Gamma-H for pure Cr and rigid-band Fe additions, normalised to the pure-Cr peak."""
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    z = np.load(npz); s = z["s"]; m = s >= 0.5
    ref = z["xi_y0.00"][m].max()
    fig, ax = plt.subplots(figsize=(3.4, 2.4))
    cols = {"0.00": "#0072B2", "0.05": "#009E73", "0.10": "#E69F00", "0.20": "#D55E00"}
    for y, col in cols.items():
        de = float(z[f"ef_y{y}"] - z["ef_y0.00"])
        lab = "pure Cr" if y == "0.00" else f"+{round(100 * float(y))}% Fe ($\\Delta E_F$ = +{de:.2f} eV)"
        ax.plot(s[m], z[f"xi_y{y}"][m] / ref, color=col, lw=1.0, label=lab)
    ax.axvline(0.952, color="0.4", lw=0.7, ls=":")
    ax.text(0.950, 0.06, "measured\n$Q$ = 0.952", fontsize=5.5, color="0.35", ha="right", transform=ax.get_xaxis_transform())
    ax.set(xlim=(0.5, 1.0), xlabel=r"$q$ along $\Gamma$H ($2\pi/a$)", ylabel=r"nesting function $\xi(q)$ (rel.)")
    ax.legend(fontsize=5.5, loc="upper left", frameon=False)
    fig.tight_layout()
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{out}.{ext}", dpi=600)
    print("wrote", out)


def full_grid_from_ibz(xml, n):
    """Band energies on the full n^3 grid (crystal coordinates, i slowest) from an nscf run on the irreducible wedge
    of the same n^3 mesh, unfolded with the 48 operations of O_h (signed permutations of Cartesian k). Added 2026-10-16:
    the full-grid nosym bands runs failed (CG left 47% of k unconverged; ParO aborted), whereas the pure_dos nscf of
    non-magnetic Cr (36^3, tetrahedra, symmetry on) is converged."""
    import itertools
    k, ev, ef = fs.read_xml(xml)
    ev = ev[0]
    root = ET.parse(xml).getroot()
    B = np.array([[float(v) for v in root.find(f"output/basis_set/reciprocal_lattice/b{i}").text.split()] for i in (1, 2, 3)])
    Binv = np.linalg.inv(B)
    ops = [np.diag(sg)[list(pm)] for pm in itertools.permutations(range(3)) for sg in itertools.product((1, -1), repeat=3)]
    owner = -np.ones(n ** 3, int)
    for j, kc in enumerate(k):
        for R in ops:
            c = np.rint((R @ kc) @ Binv * n).astype(int) % n
            if not np.allclose((R @ kc) @ Binv * n, np.rint((R @ kc) @ Binv * n), atol=1e-4):
                sys.exit("k-point not on the n^3 mesh: wrong n?")
            owner[(c[0] * n + c[1]) * n + c[2]] = j
    if (owner < 0).any():
        sys.exit(f"{(owner < 0).sum()} grid points not reached by symmetry")
    return ev[owner], ef


def analyse_ibz(xml, outdir, n=36):
    """As analyse(), on the full grid rebuilt from the irreducible wedge (full_grid_from_ibz)."""
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    ev, ef_nscf = full_grid_from_ibz(xml, n)
    root = ET.parse(xml).getroot()
    B = np.array([[float(v) for v in root.find(f"output/basis_set/reciprocal_lattice/b{i}").text.split()] for i in (1, 2, 3)])
    nel = 14.0
    ef0 = fermi_level(ev, nel)
    print(f"grid E_F {ef0:.3f} eV, nscf E_F {ef_nscf:.3f} eV")
    if abs(ef0 - ef_nscf) > 0.1:
        sys.exit("grid E_F differs from the nscf value")
    grid = ev.reshape(n, n, n, -1)
    s = np.linspace(0, 1, 201)
    cryst = np.outer(s, [0, 0, 1.0]) @ np.linalg.inv(B)
    out = {"s": s, "n": n}
    for y in Y_FE:
        ef = fermi_level(grid, nel + 2 * y)
        w = np.sum(np.exp(-0.5 * ((grid - ef) / SIGMA) ** 2), axis=-1) / (np.sqrt(2 * np.pi) * SIGMA)
        F = np.fft.fftn(w)
        xi = np.real(np.fft.ifftn(np.conj(F) * F)) / w.size ** 2
        idx = (cryst * n) % n
        i0 = np.floor(idx).astype(int); t = idx - i0
        val = 0.0
        for dx in (0, 1):
            for dy in (0, 1):
                for dz in (0, 1):
                    wgt = np.prod(np.where([dx, dy, dz], t, 1 - t), axis=1)
                    val = val + wgt * xi[(i0[:, 0] + dx) % n, (i0[:, 1] + dy) % n, (i0[:, 2] + dz) % n]
        out[f"xi_y{y:.2f}"] = val
        out[f"ef_y{y:.2f}"] = ef
        peak = s[50:][np.argmax(val[50:])]
        print(f"y_Fe = {y:.2f}: E_F shift {ef - out['ef_y0.00'] if y else 0:+.3f} eV, nesting peak at q = {peak:.3f} (2pi/a),"
              f" height {val[50:].max():.4g}")
    np.savez(outdir / "nesting.npz", **out)


if __name__ == "__main__":
    if sys.argv[1] == "inputs":
        inputs(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 32)
    elif sys.argv[1] == "analyse":
        analyse(sys.argv[2])
    elif sys.argv[1] == "plot":              # nesting.py plot <nesting.npz> <out stem>
        plot(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "analyse-ibz":       # nesting.py analyse-ibz <data-file-schema.xml> <outdir> [n]
        analyse_ibz(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 36)
