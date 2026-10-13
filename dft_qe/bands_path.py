#!/usr/bin/env python3
"""Electronic band structures along Gamma-H-N-Gamma-P-H (bcc), spin-resolved.

Pure elements (1- or 2-atom cells, settings of fermi_slices.py):
  bands_path.py pure-inputs <case> <dir>          scf.in, bands.in (Fe_fm, Cr_nm, Cr_afm)
  bands_path.py pure-plot <outfile> <dir>...      E - E_F along the path, majority/minority

Fe-rich SQS (40 atoms), unfolded onto the bcc zone with BandUPpy (effective band structure):
  bands_path.py sqs-prepare <unfold_dir> <dir>    bands.in along the path, reusing the SCF of <unfold_dir>
  bands_path.py sqs-unfold <dir>                  BandUPpy per spin -> ebs_<spin>.npz
  bands_path.py sqs-plot <outfile> <dir>...       spectral weight A(k, E), spin up / spin down
"""
import shutil
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import fermi_slices as fs  # noqa: E402

PATH_BCC = [("Γ", (0, 0, 0)), ("H", (1, 0, 0)), ("N", (0.5, 0.5, 0)), ("Γ", (0, 0, 0)),
            ("P", (0.5, 0.5, 0.5)), ("H", (1, 0, 0))]                      # Cartesian, units 2 pi / a
PATH_SC = [("Γ", (0, 0, 0)), ("X", (0.5, 0, 0)), ("M", (0.5, 0.5, 0)), ("Γ", (0, 0, 0)),
           ("R", (0.5, 0.5, 0.5)), ("X", (0.5, 0, 0))]                     # 2-atom simple-cubic cell (AFM Cr)
NSEG = 60


def path_points(path):
    pts, dist, ticks = [], [], [(0.0, path[0][0])]
    d = 0.0
    for (l0, a), (l1, b) in zip(path[:-1], path[1:]):
        a, b = np.array(a, float), np.array(b, float)
        for t in np.linspace(0, 1, NSEG, endpoint=False):
            p = a + t * (b - a)
            if pts:
                d += np.linalg.norm(p - pts[-1])
            pts.append(p); dist.append(d)
        ticks.append((d + np.linalg.norm(b - pts[-1]), l1))
    pts.append(np.array(path[-1][1], float)); d += np.linalg.norm(pts[-1] - pts[-2]); dist.append(d)
    return np.array(pts), np.array(dist), ticks


def pure_inputs(case, d):
    fs.write_inputs(case, d)                                   # scf.in (and a grid bands.in, replaced below)
    d = Path(d)
    pts, dist, _ = path_points(PATH_SC if case.startswith("Cr_afm") else PATH_BCC)
    txt = (d / "bands.in").read_text().split("K_POINTS")[0]
    (d / "bands.in").write_text(txt + f"K_POINTS tpiba\n{len(pts)}\n"
                                + "".join(f"{p[0]:.6f} {p[1]:.6f} {p[2]:.6f} 1\n" for p in pts))


def _style():
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(HERE.parent / "scripts" / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta")
    return plt


def pure_plot(outfile, dirs):
    plt = _style()
    titles = {"Fe_fm": "bcc Fe (FM)", "Cr_nm": "bcc Cr (NM)", "Cr_afm": "Cr (AFM)"}
    fig, axs = plt.subplots(1, len(dirs), figsize=(2.4 * len(dirs), 2.8), squeeze=False)
    for ax, d in zip(axs[0], dirs):
        d = Path(d); case = d.name
        k, ev, ef = fs.read_xml(d / "bands_data.xml" if (d / "bands_data.xml").exists()
                                else d / "tmp" / "fs.save" / "data-file-schema.xml")
        pts, dist, ticks = path_points(PATH_SC if case.startswith("Cr_afm") else PATH_BCC)
        for s, (col, ls) in enumerate((("#D55E00", "-"), ("#0072B2", "--"))[:ev.shape[0]]):
            for b in range(ev.shape[2]):
                ax.plot(dist, ev[s, :, b] - ef, ls, color=col if ev.shape[0] == 2 else "0.15", lw=0.8)
        ax.axhline(0, color="0.5", lw=0.5, ls=":")
        ax.set_xticks([t for t, _ in ticks]); ax.set_xticklabels([l for _, l in ticks]); ax.minorticks_off()
        for t, _ in ticks[1:-1]:
            ax.axvline(t, color="0.8", lw=0.5)
        ax.set(xlim=(0, dist[-1]), ylim=(-8, 4), ylabel=r"$E - E_{\rm F}$ (eV)")
        ax.set_title(titles.get(fs.base(case), case) + (" (exp. a)" if case.endswith("_aexp") else ""), fontsize=8)
    axs[0, 0].plot([], [], color="#D55E00", label="majority"); axs[0, 0].plot([], [], "--", color="#0072B2", label="minority")
    axs[0, 0].legend(loc="lower left", fontsize=7)
    fig.tight_layout()
    fig.savefig(outfile, dpi=600); fig.savefig(Path(outfile).with_suffix(".pdf"))
    print("wrote", outfile)


def sqs_prepare(unfold_dir, d):
    """Reuse the converged SCF of an unfold_fs.py run: copy its .save (charge density, no wavefunctions)."""
    import banduppy
    import unfold_fs as uf
    src, d = Path(unfold_dir), Path(d)
    (d / "tmp" / "uf.save").mkdir(parents=True, exist_ok=True)
    for f in (src / "tmp" / "uf.save").iterdir():
        if not f.name.startswith("wfc"):
            shutil.copy2(f, d / "tmp" / "uf.save" / f.name)
    pts, dist, ticks = path_points(PATH_BCC)
    kpc = pts @ np.linalg.inv(uf.B_PC)
    u = banduppy.Unfolding(supercell=uf.M_SC, print_log="low")
    k_full, k_uniq, K_sc, kmap, _ = u.generate_SC_Kpts_from_pc_kpts(kpointsPBZ=kpc)
    K_sc = np.asarray(K_sc)
    txt = (src / "bands.in").read_text().split("K_POINTS")[0]
    (d / "bands.in").write_text(txt + f"K_POINTS crystal\n{len(K_sc)}\n"
                                + "".join(f"{k[0]:.10f} {k[1]:.10f} {k[2]:.10f} 1\n" for k in K_sc))
    (d / "scf_source.txt").write_text(str(src.resolve()))
    z = np.load(src / "kmap.npz", allow_pickle=True)
    np.savez(d / "pathmap.npz", k_full=np.asarray(k_full), K_sc=K_sc, kmap=np.array([kmap], dtype=object),
             dist=dist, ticks=np.array([t for t, _ in ticks]), labels=np.array([l for _, l in ticks]),
             x_cr=float(z["x_cr"]))
    print(d, "PC k", len(k_full), "SC K", len(K_sc))


def sqs_unfold(d):
    import banduppy
    import unfold_fs as uf
    d = Path(d)
    z = np.load(d / "pathmap.npz", allow_pickle=True)
    u = banduppy.Unfolding(supercell=uf.M_SC, print_log="low")
    for spin in ("up", "dw"):
        bands, kpts = u.Unfold(PBZ_kpts_list_full=z["k_full"], SBZ_kpts_list=z["K_sc"], SBZ_PBZ_kpts_map=z["kmap"][0],
                               ab_initio_code="qe", qe_keywards={"output_file_dir": str(d / "tmp"),
                               "save_file_prefix": "uf", "unfold_spin_channel": spin},
                               unfold_kpts_in_batch=True, kpt_batch_size=20)
        np.savez(d / f"ebs_{spin}.npz", kpts=np.asarray(kpts), bands=np.asarray(bands))


def uf_fermi(d):
    """E_F (eV) of the SCF this run reuses (the source unfold directory's scf.out), if known."""
    import re
    src = Path(d) / "scf_source.txt"
    if not src.exists():
        return None
    t = Path(src.read_text().strip(), "scf.out").read_text()
    return float(re.findall(r"the Fermi energy is\s+([-\d.]+) ev", t)[-1])


def sqs_plot(outfile, dirs, eta=0.06):
    plt = _style()
    egrid = np.linspace(-8, 4, 481)
    fig, axs = plt.subplots(2, len(dirs), figsize=(2.4 * len(dirs), 4.8), squeeze=False)
    for j, d in enumerate(dirs):
        z = np.load(Path(d) / "pathmap.npz", allow_pickle=True)
        dist = z["dist"]
        for i, (spin, lab) in enumerate((("up", "spin up"), ("dw", "spin down"))):
            u = np.load(Path(d) / f"ebs_{spin}.npz")["bands"]   # BandUPpy: k index, k on path, energy (eV), weight
            ef = uf_fermi(d)
            A = np.zeros((len(egrid), len(dist)))
            for ik, _, e, w in u[:, :4]:
                e = e - ef if ef is not None else e
                A[:, int(ik)] += w * (eta / np.pi) / ((egrid - e) ** 2 + eta ** 2)
            ax = axs[i, j]
            ax.imshow(A, extent=(0, dist[-1], egrid[0], egrid[-1]), origin="lower", aspect="auto", cmap="Greys",
                      vmin=0, vmax=np.percentile(A, 99.5))
            ax.axhline(0, color="#D55E00", lw=0.6, ls=":")
            ax.set_xticks(z["ticks"]); ax.set_xticklabels(z["labels"]); ax.minorticks_off()
            ax.set_ylabel(r"$E - E_{\rm F}$ (eV)")
            ax.set_title(f"{round(100 * float(z['x_cr']))}% Cr, {lab}", fontsize=7)
    fig.tight_layout()
    fig.savefig(outfile, dpi=600); fig.savefig(Path(outfile).with_suffix(".pdf"))
    print("wrote", outfile)


if __name__ == "__main__":
    c = sys.argv[1]
    if c == "pure-inputs":
        pure_inputs(sys.argv[2], sys.argv[3])
    elif c == "pure-plot":
        pure_plot(sys.argv[2], sys.argv[3:])
    elif c == "sqs-prepare":
        sqs_prepare(sys.argv[2], sys.argv[3])
    elif c == "sqs-unfold":
        sqs_unfold(sys.argv[2])
    elif c == "sqs-plot":
        sqs_plot(sys.argv[2], sys.argv[3:])
