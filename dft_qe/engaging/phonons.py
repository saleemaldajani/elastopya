#!/usr/bin/env python3
"""Phonon dispersions of bcc Fe (FM) and Cr (NM, commensurate AFM) from DFPT (QE ph.x, q2r.x, matdyn.x).

  phonons.py inputs <case> <dir> [nq]     scf.in (conv 1e-12), ph.in (nq^3 grid), q2r.in, matdyn.in (path)
  phonons.py plot <outfile> <dir>...      frequencies (THz) along Gamma-H-N-Gamma-P-H (Cr AFM: Gamma-X-M-Gamma-R-X)
Same pseudopotentials, cutoff and smearing as the elastic-constant runs (fermi_slices.CASES).
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bands_path as bp  # noqa: E402
import fermi_slices as fs  # noqa: E402

CM1_THZ = 0.0299792458


def inputs(case, d, nq=6):
    fs.write_inputs(case, d)
    d = Path(d)
    scf = (d / "scf.in").read_text().replace("conv_thr = 1e-9", "conv_thr = 1e-12")
    (d / "scf.in").write_text(scf)
    (d / "bands.in").unlink()
    c = fs.CASES[case]
    masses = "".join(f"  amass({i + 1}) = {m}\n" for i, (_, m, _) in enumerate(c["species"]))
    (d / "ph.in").write_text(f"phonons\n&inputph\n  prefix = 'fs', outdir = './tmp', fildyn = 'dyn',\n"
                             f"  tr2_ph = 1e-15, ldisp = .true., nq1 = {nq}, nq2 = {nq}, nq3 = {nq},\n{masses}/\n")
    (d / "q2r.in").write_text("&input\n  fildyn = 'dyn', zasr = 'crystal', flfrc = 'fc'\n/\n")
    path = bp.PATH_SC if case == "Cr_afm" else bp.PATH_BCC
    lines = "".join(f"{q[0]} {q[1]} {q[2]} {bp.NSEG}\n" for _, q in path[:-1]) + f"{path[-1][1][0]} {path[-1][1][1]} {path[-1][1][2]} 1\n"
    (d / "matdyn.in").write_text(f"&input\n  asr = 'crystal', flfrc = 'fc', flfrq = 'freq', q_in_band_form = .true.,\n"
                                 f"{masses}/\n{len(path)}\n{lines}")


def plot(outfile, dirs):
    plt = bp._style()
    titles = {"Fe_fm": "bcc Fe (FM)", "Cr_nm": "bcc Cr (NM)", "Cr_afm": "Cr (AFM)"}
    fig, axs = plt.subplots(1, len(dirs), figsize=(2.4 * len(dirs), 2.6), squeeze=False)
    for ax, d in zip(axs[0], dirs):
        d = Path(d); case = d.name
        gp = np.loadtxt(d / "freq.gp")
        x, w = gp[:, 0], gp[:, 1:] * CM1_THZ
        for b in range(w.shape[1]):
            ax.plot(x, w[:, b], color="0.15", lw=0.9)
        _, _, ticks = bp.path_points(bp.PATH_SC if case == "Cr_afm" else bp.PATH_BCC)
        scale = x[-1] / ticks[-1][0]
        ax.set_xticks([t * scale for t, _ in ticks]); ax.set_xticklabels([l for _, l in ticks]); ax.minorticks_off()
        for t, _ in ticks[1:-1]:
            ax.axvline(t * scale, color="0.8", lw=0.5)
        ax.axhline(0, color="0.6", lw=0.5)
        ax.set(xlim=(0, x[-1]), ylabel="Frequency (THz)")
        ax.set_title(titles.get(case, case), fontsize=8)
    fig.tight_layout()
    fig.savefig(outfile, dpi=600); fig.savefig(Path(outfile).with_suffix(".pdf"))
    print("wrote", outfile)


if __name__ == "__main__":
    if sys.argv[1] == "inputs":
        inputs(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 6)
    elif sys.argv[1] == "plot":
        plot(sys.argv[2], sys.argv[3:])
