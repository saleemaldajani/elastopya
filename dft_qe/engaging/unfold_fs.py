#!/usr/bin/env python3
"""Unfolded (Bloch spectral function) Fermi surfaces of Fe-rich SQS supercells, for comparison with
Korzhavyi et al., PRB 79, 054202 (2009), Fig. 6.

The 40-atom SQS (2x2x5 conventional bcc cells) is a 40-fold supercell of the 1-atom bcc cell. Its bands
are computed at the supercell K-points that unfold onto a grid covering the (001) zone section through
Gamma (square with H at the corners, as in Korzhavyi Fig. 6), and unfolded with BandUPpy for each spin.

  unfold_fs.py prepare <sqs.vasp> <dir> [--a 2.83] [--n 61]   write scf.in, bands.in, kmap.npz
  unfold_fs.py unfold <dir>                                      run BandUPpy per spin -> unfolded_<spin>.npz
  unfold_fs.py plot <outfile> <dir>...                           Korzhavyi-style spectral-function maps

The lattice parameter is fixed (default 2.83 A) and ions are not relaxed: the figure is about the Fermi
surface topology, not total energies.
"""
import os
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
PSEUDO = os.environ.get("QE_PSEUDO", "/orcd/pool/004/sdajani/fecr_qe/pseudo")
UPF = {"Fe": "Fe.pbe-spn-kjpaw_psl.0.2.1.UPF", "Cr": "Cr.pbe-spn-kjpaw_psl.1.0.0.UPF"}
MASS = {"Fe": 55.845, "Cr": 51.996}
VAL = {"Fe": 16, "Cr": 14}                    # valence electrons of the spn PAW potentials
# supercell vectors in terms of the bcc primitive vectors p1 = a/2(-1,1,1), p2 = a/2(1,-1,1), p3 = a/2(1,1,-1):
# 2a x = 2(p2 + p3), 2a y = 2(p1 + p3), 5a z = 5(p1 + p2)
M_SC = np.array([[0, 2, 2], [2, 0, 2], [5, 5, 0]])
B_PC = np.array([[0, 1, 1], [1, 0, 1], [1, 1, 0]], float)     # PC reciprocal vectors, units 2 pi / a


def read_poscar(path):
    L = Path(path).read_text().split("\n")
    scale = float(L[1]); cell = np.array([[float(v) for v in L[i].split()] for i in (2, 3, 4)]) * scale
    names, counts = L[5].split(), [int(v) for v in L[6].split()]
    syms = [n for n, c in zip(names, counts) for _ in range(c)]
    frac = np.array([[float(v) for v in L[8 + i].split()[:3]] for i in range(len(syms))])
    return cell, syms, frac


def grid_pq(n):
    pq = np.linspace(-1, 1, n)
    P, Q = np.meshgrid(pq, pq, indexing="xy")
    kx, ky = (P + Q) / 2, (Q - P) / 2                          # Cartesian, units 2 pi / a
    return pq, np.stack([kx.ravel(), ky.ravel(), np.zeros(kx.size)], -1)


def qe_input(calc, cell, syms, frac, kblock, nbnd):
    sp = sorted(set(syms))
    mag = {"Fe": 0.6, "Cr": -0.3}
    return (f"&control\n  calculation = '{calc}', prefix = 'uf', outdir = './tmp', pseudo_dir = '{PSEUDO}'\n/\n"
            f"&system\n  ibrav = 0, nat = {len(syms)}, ntyp = {len(sp)}, ecutwfc = 60, ecutrho = 480,\n"
            f"  occupations = 'smearing', smearing = 'mv', degauss = 0.02, nspin = 2, nbnd = {nbnd}\n"
            + "".join(f"  starting_magnetization({i + 1}) = {mag[s]}\n" for i, s in enumerate(sp))
            + ("  nosym = .true., noinv = .true.\n" if calc == "bands" else "")
            + "/\n&electrons\n  conv_thr = 1e-8, mixing_beta = 0.2, mixing_mode = 'local-TF', electron_maxstep = 300"
            + (", diago_full_acc = .true., diagonalization = 'cg'" if calc == "bands" else "") + "\n/\n"
            "ATOMIC_SPECIES\n" + "".join(f"{s} {MASS[s]} {UPF[s]}\n" for s in sp)
            + "CELL_PARAMETERS angstrom\n" + "".join(" ".join(f"{v:.10f}" for v in r) + "\n" for r in cell)
            + "ATOMIC_POSITIONS crystal\n"
            + "".join(f"{s} {f[0]:.10f} {f[1]:.10f} {f[2]:.10f}\n" for s, f in zip(syms, frac)) + kblock)


def prepare(sqs, d, a=2.83, n=61):
    import banduppy
    d = Path(d); d.mkdir(parents=True, exist_ok=True)
    cell0, syms, frac = read_poscar(sqs)
    cell = np.diag([2 * a, 2 * a, 5 * a])
    assert np.allclose(cell0 / np.abs(cell0).max(), cell / cell.max(), atol=1e-6), "expected a 2x2x5 cell"
    nel = sum(VAL[s] for s in syms)
    nbnd = int(nel / 2 * 1.2) + 8
    pq, kcart = grid_pq(n)
    kpc = kcart @ np.linalg.inv(B_PC)                           # reduced coordinates of the PC reciprocal basis
    uf = banduppy.Unfolding(supercell=M_SC, print_log="low")
    k_full, k_uniq, K_sc, kmap, _ = uf.generate_SC_Kpts_from_pc_kpts(kpointsPBZ=kpc)
    K_sc = np.asarray(K_sc)
    (d / "scf.in").write_text(qe_input("scf", cell, syms, frac, "K_POINTS automatic\n4 4 2 0 0 0\n", nbnd))
    kb = f"K_POINTS crystal\n{len(K_sc)}\n" + "".join(f"{k[0]:.10f} {k[1]:.10f} {k[2]:.10f} 1\n" for k in K_sc)
    (d / "bands.in").write_text(qe_input("bands", cell, syms, frac, kb, nbnd))
    np.savez(d / "kmap.npz", pq=pq, k_full=np.asarray(k_full), K_sc=K_sc, kmap=np.array([kmap], dtype=object),
             n=n, x_cr=syms.count("Cr") / len(syms))
    print(f"{d}: {len(syms)} atoms, x_Cr = {syms.count('Cr') / len(syms):.3f}, nbnd = {nbnd}, "
          f"PC k = {len(k_full)}, SC K = {len(K_sc)}")


def unfold(d):
    import banduppy
    d = Path(d)
    z = np.load(d / "kmap.npz", allow_pickle=True)
    uf = banduppy.Unfolding(supercell=M_SC, print_log="low")
    for spin in ("up", "dw"):
        bands, kpts = uf.Unfold(PBZ_kpts_list_full=z["k_full"], SBZ_kpts_list=z["K_sc"],
                                SBZ_PBZ_kpts_map=z["kmap"][0], ab_initio_code="qe",
                                qe_keywards={"output_file_dir": str(d / "tmp"), "save_file_prefix": "uf",
                                             "unfold_spin_channel": spin},
                                unfold_kpts_in_batch=True, kpt_batch_size=20)
        np.savez(d / f"unfolded_{spin}.npz", kpts=np.asarray(kpts), bands=np.asarray(bands))
        print(spin, np.asarray(bands).shape)


def fermi_energy(d):
    """E_F (eV) of the SCF run in <d> (from scf.out)."""
    import re
    t = (Path(d) / "scf.out").read_text()
    m = re.findall(r"the Fermi energy is\s+([-\d.]+) ev", t)
    return float(m[-1])


def spectral(d, spin, eta=0.08):
    """A(k, E_F) on the (p, q) grid from the unfolded weights: sum_n w_n (eta/pi)/((E_n - E_F)^2 + eta^2).
    BandUPpy returns energies relative to E_F (it reads the Fermi energy from the QE output)."""
    z = np.load(Path(d) / "kmap.npz", allow_pickle=True)
    u = np.load(Path(d) / f"unfolded_{spin}.npz")["bands"]   # BandUPpy: k index, k on path, energy (eV), weight
    n = int(z["n"])
    A = np.zeros(n * n)
    ik, e, w = u[:, 0].astype(int), u[:, 2], u[:, 3]
    e = e - fermi_energy(d) if np.median(e) > 5 else e      # absolute energies -> relative to E_F
    np.add.at(A, ik, w * (eta / np.pi) / (e ** 2 + eta ** 2))
    return z["pq"], A.reshape(n, n), float(z["x_cr"])


def plot(outfile, dirs):
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(HERE.parent / "scripts" / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta")
    fig, axs = plt.subplots(2, len(dirs), figsize=(2.2 * len(dirs), 4.4), squeeze=False)
    for j, d in enumerate(dirs):
        for i, (spin, lab) in enumerate((("up", "spin up"), ("dw", "spin down"))):
            pq, A, x = spectral(d, spin)
            ax = axs[i, j]
            ax.imshow(A, extent=(-1, 1, -1, 1), origin="lower", cmap="Greys", vmin=0, vmax=np.percentile(A, 99.5))
            ax.set_xticks([]); ax.set_yticks([]); ax.minorticks_off()
            for (px, py, t) in ((-1, -1, "H"), (1, -1, "H"), (-1, 1, "H"), (1, 1, "H"),
                                (0, -1, "N"), (0, 1, "N"), (-1, 0, "N"), (1, 0, "N")):
                ax.text(px * 1.08, py * 1.08, t, ha="center", va="center", fontsize=7, clip_on=False)
            ax.text(0.03, 0.03, "Γ", fontsize=7)
            ax.set_title(f"{round(100 * x)}% Cr, {lab}", fontsize=7, pad=14)
    fig.tight_layout(h_pad=1.5)
    fig.savefig(outfile, dpi=600); fig.savefig(Path(outfile).with_suffix(".pdf"))
    print("wrote", outfile)


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "prepare":
        args = sys.argv[2:]
        a = float(args[args.index("--a") + 1]) if "--a" in args else 2.83
        n = int(args[args.index("--n") + 1]) if "--n" in args else 61
        prepare(args[0], args[1], a, n)
    elif cmd == "unfold":
        unfold(sys.argv[2])
    elif cmd == "plot":
        plot(sys.argv[2], sys.argv[3:])
