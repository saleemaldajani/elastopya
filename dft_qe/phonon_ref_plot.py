#!/usr/bin/env python3
"""Our DFPT phonons of bcc Fe and Cr against experiment.

  phonon_ref_plot.py <outfile>              Fe from DFPT (the original version of the figure)
  phonon_ref_plot.py <outfile> --fe-fd      Fe from finite displacements (phonopy_cubic.py, 3x3x3 cubic supercell)
                                            at the PBE-relaxed and the experimental lattice parameter. Used since
                                            2026-10-14: DFPT (ph.x) gives FM Fe frequencies ~60% too high, while a
                                            frozen phonon with identical settings agrees with neutron data (WORKFLOW_NOTES.md).

(a) Fe (FM, DFPT) dispersion along Gamma-H-N-Gamma-P-H with the room-temperature neutron data of Minkiewicz,
    Shirane & Nathans, Phys. Rev. 162, 528 (1967) (tabulated).
(b) Cr (NM, k = 24^3) with the 300 K neutron data of Shaw & Muhlestein, PRB 4, 969 (1971) (digitized).
(c) Phonon DOS of Fe: DFPT (40^3 q-mesh, Gaussian-smoothed 0.15 THz) vs NRIXS at 30 K and 300 K (Mauger et al.,
    PRB 90, 064303 (2014)) and the force-constant DOS of Fultz et al., PRB 52, 3280 (1995).
(d) Phonon DOS of Cr: DFPT (NM and AFM) vs Fultz et al. (1995).
All DOS normalised to unit area.
"""
import csv
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ENG, REF = HERE / "engaging" / "phonon", HERE / "refdata"
CM1_THZ = 0.0299792458
S2, S3 = np.sqrt(2), np.sqrt(3)
X = dict(G=0.0, H=1.0, N=1.0 + S2 / 2, G2=1.0 + S2, P=1.0 + S2 + S3 / 2, H2=1.0 + S2 + S3)


def path_x(direction, zeta):
    """Position on the Gamma-H-N-Gamma-P-H path (units 2 pi/a) of a point given in the experimental convention."""
    d = direction.replace(" ", "")
    if d in ("[z00]", "[00z]"):
        return X["G"] + zeta
    if d in ("[zz1]",):                                   # H (0,0,1) -> N' (1/2,1/2,1)
        return X["H"] + zeta * S2
    if d in ("[zz0]", "[0zz]"):                           # Gamma -> N, drawn on the N-Gamma segment
        return X["G2"] - zeta * S2
    if d == "[zzz]":                                      # Gamma -> P -> H
        return X["G2"] + zeta * S3
    return None


def rows(name):
    with open(REF / name) as fh:
        r = list(csv.DictReader(l for l in fh if not l.startswith("#")))
    return r


def ph_dir(case):
    """DFPT runs in engaging/phonon/<case>; finite-displacement runs in phonopy_derived/<case> (same file formats)."""
    d = HERE / "phonopy_derived" / case
    return d if d.exists() else ENG / case


def dispersion(ax, case, refname, label, variants=None):
    """variants: list of (case, colour, linestyle, label); default: the one DFPT case in blue."""
    for c, col, ls, lab in (variants or [(case, "#2A9DF4", "-", "QE DFPT")]):
        gp = np.loadtxt(ph_dir(c) / "freq.gp")
        x, w = gp[:, 0], gp[:, 1:] * CM1_THZ
        for b in range(w.shape[1]):
            ax.plot(x, w[:, b], ls, color=col, lw=1.0, label=lab if b == 0 else None)
    pts = [(path_x(r["direction"], float(r["zeta"])), float(r["nu_THz"])) for r in rows(refname)]
    pts = np.array([p for p in pts if p[0] is not None])
    ax.plot(pts[:, 0], pts[:, 1], "o", mfc="white", mec="0.1", ms=2.6, mew=0.6, label=label)
    for k, v in X.items():
        ax.axvline(v, color="0.85", lw=0.5)
    ax.set_xticks(list(X.values())); ax.set_xticklabels(["Γ", "H", "N", "Γ", "P", "H"]); ax.minorticks_off()
    ax.set(xlim=(0, X["H2"]), ylim=(0, 11.5), ylabel="frequency (THz)")
    ax.legend(fontsize=6, loc="lower right", frameon=True, framealpha=0.9)


def dfpt_dos(case, sigma=0.15):
    a = np.loadtxt(ph_dir(case) / "dos.dat")
    f, g = a[:, 0] * CM1_THZ, a[:, 1]
    grid = np.linspace(0, 12, 1201)
    gi = np.interp(grid, f, g, left=0, right=0)
    k = np.exp(-0.5 * (np.arange(-60, 61) * 0.01 / sigma) ** 2)
    gs = np.convolve(gi, k / k.sum(), mode="same")
    return grid, gs / np.trapezoid(gs, grid)


def main(outfile, fe_fd=False):
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(HERE.parent / "scripts" / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.0), gridspec_kw={"width_ratios": [1.4, 1]})
    fe_var = [("Fe_cub_aexp", "#2A9DF4", "-", "QE finite displ., exp. $a$"), ("Fe_cub", "#0B5D3B", "--", "QE finite displ., relaxed $a$")] if fe_fd else None
    dispersion(axs[0, 0], "Fe_fm", "fe_phonon_minkiewicz1967.csv", "neutron, 295 K (Minkiewicz 1967)", variants=fe_var)
    axs[0, 0].set_title("(a) bcc Fe (FM" + (", finite displacements)" if fe_fd else ")"), fontsize=8, loc="left")
    dispersion(axs[1, 0], "Cr_nm_v5_k24", "cr_phonon_shaw1971.csv", "neutron, 300 K (Shaw & Muhlestein 1971)")
    axs[1, 0].set_title("(b) bcc Cr (NM, $k=24^3$)", fontsize=8, loc="left")
    ax = axs[0, 1]
    for c, col, ls, lab in (fe_var or [("Fe_fm", "#2A9DF4", "-", "QE DFPT")]):
        f, g = dfpt_dos(c); ax.plot(f, g, ls, color=col, lw=1.2, label=lab)
    mau = rows("fe_phonon_dos_mauger2014.csv")
    for T, ls in (("30", "-"), ("300", "--")):
        d = np.array([[float(r["nu_THz"]), float(r["g_perTHz"])] for r in mau if r["T_K"] == T])
        if len(d):
            ax.plot(d[:, 0], d[:, 1], ls, color="0.1", lw=0.9, label=f"NRIXS {T} K (Mauger 2014)")
    fu = rows("fecr_phonon_dos_fultz1995.csv")
    fd = np.array([[float(r["nu_THz"]), float(r["Fe"]), float(r["Cr"])] for r in fu])
    ax.plot(fd[:, 0], fd[:, 1], ":", color="#D55E00", lw=1.0, label="Fultz 1995")
    ax.set(xlim=(0, 11.5), ylim=(0, None), xlabel=r"$\nu$ (THz)", ylabel=r"$g(\nu)$ (THz$^{-1}$)")
    ax.set_title("(c) Fe phonon DOS", fontsize=8, loc="left"); ax.legend(fontsize=5.5, loc="upper left")
    ax = axs[1, 1]
    for case, col, lab in (("Cr_nm_v5_k24", "#2A9DF4", "QE DFPT, NM"), ("Cr_afm", "#0B5D3B", "QE DFPT, AFM (2-atom)")):
        if (ENG / case / "dos.dat").exists():
            f, g = dfpt_dos(case); ax.plot(f, g, color=col, lw=1.1, label=lab)
    ax.plot(fd[:, 0], fd[:, 2], ":", color="#D55E00", lw=1.0, label="Fultz 1995")
    ax.set(xlim=(0, 11.5), ylim=(0, None), xlabel=r"$\nu$ (THz)", ylabel=r"$g(\nu)$ (THz$^{-1}$)")
    ax.set_title("(d) Cr phonon DOS", fontsize=8, loc="left"); ax.legend(fontsize=5.5, loc="upper left")
    fig.tight_layout()
    fig.savefig(outfile, dpi=600); fig.savefig(Path(outfile).with_suffix(".pdf")); fig.savefig(Path(outfile).with_suffix(".svg"))
    print("wrote", outfile)


if __name__ == "__main__":
    main(sys.argv[1], fe_fd="--fe-fd" in sys.argv)
