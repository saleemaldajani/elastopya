#!/usr/bin/env python3
"""Phonon dispersions of bcc Fe (FM) and Cr (NM, commensurate AFM) from DFPT (QE ph.x, q2r.x, matdyn.x).

  phonons.py inputs <case> <dir> [nq]     scf.in (conv 1e-12), ph.in (nq^3 grid), q2r.in, matdyn.in (path)
  phonons.py check <dir>... [--allow-imaginary]   expected-stable systems must have no imaginary mode
  phonons.py plot <outfile> <dir>... [--allow-imaginary]   (checks first) frequencies (THz) along Gamma-H-N-Gamma-P-H (Cr AFM: Gamma-X-M-Gamma-R-X)
  phonons.py dos-input <dir> [nq]         matdyn_dos.in: phonon DOS on an nq^3 mesh (default 40) -> dos.dat
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
                             f"  tr2_ph = 1e-15, ldisp = .true., nq1 = {nq}, nq2 = {nq}, nq3 = {nq},\n"
                             f"  alpha_mix(1) = 0.2, nmix_ph = 12,\n{masses}/\n")   # gentler mixing: default 0.7 diverged for Cr
    (d / "q2r.in").write_text("&input\n  fildyn = 'dyn', zasr = 'crystal', flfrc = 'fc'\n/\n")
    path = bp.PATH_SC if case.startswith("Cr_afm") else bp.PATH_BCC
    lines = "".join(f"{q[0]} {q[1]} {q[2]} {bp.NSEG}\n" for _, q in path[:-1]) + f"{path[-1][1][0]} {path[-1][1][1]} {path[-1][1][2]} 1\n"
    (d / "matdyn.in").write_text(f"&input\n  asr = 'crystal', flfrc = 'fc', flfrq = 'freq', q_in_band_form = .true.,\n"
                                 f"{masses}/\n{len(path)}\n{lines}")


NEG_TOL_THZ = -0.05   # numerical-noise tolerance near Gamma
# Systems expected to be dynamically stable: imaginary modes there signal a numerical problem.
# A deliberately unstable phase (e.g. a hypothetical structure) is plotted with --allow-imaginary,
# and its imaginary branches are then a result, shown below zero.
EXPECT_STABLE = {"Fe_fm", "Cr_nm", "Cr_afm"}   # also their _aexp versions (fermi_slices.base)


def check(d, allow_imaginary=False):
    """Lowest frequency (THz) along the path. For an expected-stable system, imaginary modes beyond
    NEG_TOL_THZ stop the plot unless allow_imaginary is set."""
    w = np.loadtxt(Path(d) / "freq.gp")[:, 1:] * CM1_THZ
    wmin = float(w.min())
    name = Path(d).name
    print(f"{name}: lowest frequency {wmin:.3f} THz")
    if wmin < NEG_TOL_THZ and fs.base(name.split("_v")[0]) in EXPECT_STABLE and not allow_imaginary:
        raise SystemExit(f"{name}: imaginary modes down to {wmin:.2f} THz in a phase expected to be stable; "
                         "fix the calculation (or pass --allow-imaginary if instability is the result)")
    return wmin


def plot(outfile, dirs, allow_imaginary=False):
    for d in dirs:
        check(d, allow_imaginary)
    plt = bp._style()
    titles = {"Fe_fm": "bcc Fe (FM)", "Cr_nm": "bcc Cr (NM)", "Cr_afm": "Cr (AFM)",
              "Fe_cub": "bcc Fe (FM, finite displ.)"}   # Fe_cub: phonopy_cubic.py (DFPT is too stiff for FM Fe)
    fig, axs = plt.subplots(1, len(dirs), figsize=(2.4 * len(dirs), 2.6), squeeze=False)
    for ax, d in zip(axs[0], dirs):
        d = Path(d); case = d.name.split("_v")[0]
        gp = np.loadtxt(d / "freq.gp")
        x, w = gp[:, 0], gp[:, 1:] * CM1_THZ
        for b in range(w.shape[1]):
            ax.plot(x, w[:, b], color="0.15", lw=0.9)
        _, _, ticks = bp.path_points(bp.PATH_SC if case.startswith("Cr_afm") else bp.PATH_BCC)
        scale = x[-1] / ticks[-1][0]
        ax.set_xticks([t * scale for t, _ in ticks]); ax.set_xticklabels([l for _, l in ticks]); ax.minorticks_off()
        for t, _ in ticks[1:-1]:
            ax.axvline(t * scale, color="0.8", lw=0.5)
        ax.axhline(0, color="0.6", lw=0.5)
        ax.set(xlim=(0, x[-1]), ylabel="Frequency (THz)")
        suffix = d.name[len(case) + 1:] if d.name != case else ""
        suffix = {"v5_k24": "($k=24^3$)", "v3_k20": "($k=20^3$)", "v4_deg04": "(smearing 0.04 Ry)"}.get(suffix, suffix)
        ax.set_title(titles.get(fs.base(case), case) + (" (exp. a)" if case.endswith("_aexp") else "") + (f" {suffix}" if suffix else ""), fontsize=8)
    fig.tight_layout()
    fig.savefig(outfile, dpi=600); fig.savefig(Path(outfile).with_suffix(".pdf")); fig.savefig(Path(outfile).with_suffix(".svg"))
    print("wrote", outfile)


def dos_input(d, nq=40):
    """matdyn input for the phonon DOS from the same force constants (fc) as the dispersion."""
    d = Path(d)
    masses = "".join(l + "\n" for l in (d / "matdyn.in").read_text().splitlines() if "amass" in l)
    (d / "matdyn_dos.in").write_text(f"&input\n  asr = 'crystal', flfrc = 'fc', dos = .true., fldos = 'dos.dat',\n"
                                     f"  nk1 = {nq}, nk2 = {nq}, nk3 = {nq}, deltaE = 1.0,\n{masses}/\n")


if __name__ == "__main__":
    if sys.argv[1] == "dos-input":
        dos_input(sys.argv[2], int(sys.argv[3]) if len(sys.argv) > 3 else 40)
    elif sys.argv[1] == "inputs":
        inputs(sys.argv[2], sys.argv[3], int(sys.argv[4]) if len(sys.argv) > 4 else 6)
    elif sys.argv[1] in ("check", "plot"):
        allow = "--allow-imaginary" in sys.argv
        args = [a for a in sys.argv[2:] if a != "--allow-imaginary"]
        if sys.argv[1] == "check":
            for d in args:
                check(d, allow)
        else:
            plot(args[0], args[1:], allow)
