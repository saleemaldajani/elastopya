#!/usr/bin/env python3
"""Spin-density wave (SDW) and accompanying charge-density wave (CDW) of bcc Cr in a 1x1x20 supercell of the
conventional cubic cell (40 atoms, layers n = 0..39 at z = n a/2): incommensurability delta = 1/20, close to
the experimental 0.048 (Fawcett, RMP 60, 209 (1988)). The same cell is computed non-magnetic (NM) and
commensurate antiferromagnetic (AFM) for reference. Settings as the other QE runs (PBE-PAW pslibrary,
60/480 Ry, MV 0.02 Ry), room-temperature experimental lattice parameter a = 2.8846 A, k = 16x16x1.

  sdw_cr.py inputs <dir>                    nm/, afm/, sdw/ with scf.in (+ pp.in, avg.in for the planar average)
  sdw_cr.py plot <outfile> <dir>           moments m(z), site charges, planar-averaged density and forces
"""
import re
import sys
from pathlib import Path

import numpy as np

A = 2.8846
NL = 40                       # layers (2 per conventional cell)
PSEUDO = "/orcd/pool/004/sdajani/fecr_qe/pseudo"
UPF = "Cr.pbe-spn-kjpaw_psl.1.0.0.UPF"


def envelope(n):
    """Starting SDW: AF alternation times a cosine envelope of period 20 a (nodes between layers)."""
    return (-1) ** n * np.cos(2 * np.pi * (n + 0.5) / NL)


def write(d, kind):
    d = Path(d); d.mkdir(parents=True, exist_ok=True)
    n = np.arange(NL)
    if kind == "nm":
        species, start = ["Cr"] * NL, {}
    elif kind == "afm":
        species = ["Cr1" if k % 2 == 0 else "Cr2" for k in n]
        start = {"Cr1": 0.3, "Cr2": -0.3}
    else:   # sdw: 5 magnitude levels x 2 signs = 10 species (QE ntyp <= 10)
        s = envelope(n)
        lev = np.clip(np.ceil(np.abs(s) * 5), 1, 5).astype(int)
        species = [f"Cr{'p' if v > 0 else 'm'}{l}" for v, l in zip(s, lev)]
        start = {sp: (0.3 if sp[2] == "p" else -0.3) * int(sp[3]) / 5 for sp in sorted(set(species))}
    sp_list = sorted(set(species))
    mag = "" if kind == "nm" else "".join(f"  starting_magnetization({i + 1}) = {start[s]:.3f}\n" for i, s in enumerate(sp_list))
    pos = "".join(f"{sp} {0.5 * (k % 2):.4f} {0.5 * (k % 2):.4f} {k / NL:.6f}\n" for k, sp in zip(n, species))
    (d / "scf.in").write_text(
        f"&control\n  calculation = 'scf', prefix = 'sdw', outdir = './tmp', pseudo_dir = '{PSEUDO}',\n"
        f"  tprnfor = .true., tstress = .true.\n/\n"
        f"&system\n  ibrav = 0, nat = {NL}, ntyp = {len(sp_list)}, ecutwfc = 60, ecutrho = 480,\n"
        f"  occupations = 'smearing', smearing = 'mv', degauss = 0.02, nspin = {1 if kind == 'nm' else 2}\n{mag}/\n"
        f"&electrons\n  conv_thr = 1e-9, mixing_beta = 0.1, mixing_mode = 'local-TF', electron_maxstep = 500\n/\n"
        "ATOMIC_SPECIES\n" + "".join(f"{s} 51.996 {UPF}\n" for s in sp_list) +
        f"CELL_PARAMETERS angstrom\n{A} 0 0\n0 {A} 0\n0 0 {A * NL / 2}\nATOMIC_POSITIONS crystal\n{pos}"
        "K_POINTS automatic\n16 16 1 0 0 0\n")
    (d / "pp.in").write_text("&inputpp\n  prefix = 'sdw', outdir = './tmp', filplot = 'rho.dat', plot_num = 0\n/\n")
    (d / "avg.in").write_text("1\nrho.dat\n1.0\n4000\n3\n3.0\n")   # planar average along z (average.x)


def inputs(d):
    for kind in ("nm", "afm", "sdw"):
        write(Path(d) / kind, kind)
    print("wrote", d)


def sites(out):
    t = Path(out).read_text()
    blk = t.split("Magnetic moment per site")[-1].split("\n\n")[0] if "Magnetic moment per site" in t else ""
    ch = np.array([float(v) for v in re.findall(r"charge=\s*([-\d.]+)", blk)])
    mg = np.array([float(v) for v in re.findall(r"magn=\s*([-\d.]+)", blk)])
    fz = np.array([float(v) for v in re.findall(r"atom\s+\d+\s+type\s+\d+\s+force =\s+\S+\s+\S+\s+(\S+)",
                                                 t.split("Forces acting on atoms")[-1])][:NL]) if "Forces acting" in t else None
    e = float(re.findall(r"^!\s+total energy\s+=\s+([-\d.]+)", t, re.M)[-1])
    return ch, mg, fz, e


def plot(outfile, d):
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    d = Path(d)
    z = np.arange(NL) * A / 2
    fig, axs = plt.subplots(3, 1, figsize=(5.2, 5.0), sharex=True)
    res = {k: sites(d / k / "scf.out") for k in ("nm", "afm", "sdw") if (d / k / "scf.out").exists()}
    if "sdw" in res:
        ch, mg, fz, e = res["sdw"]
        axs[0].plot(z, mg, "o-", color="#0072B2", ms=3, lw=0.8, label="SDW")
        axs[0].plot(z, np.abs(mg), "-", color="0.6", lw=0.6, label="|m|")
        if "afm" in res:
            axs[0].plot(z, res["afm"][1], ".", color="0.4", ms=2, label="commensurate AFM")
        axs[0].set_ylabel(r"$m$ ($\mu_{\rm B}$)")
        axs[0].set_ylim(-1.35, 2.0)                          # head-room so the legend does not cover the moments
        axs[0].legend(fontsize=6, ncol=3, loc="upper center", frameon=False)
        # the NM run prints no site charges (no magnetic-moment block), so the CDW is measured from the mean
        # site charge of the SDW cell itself (2026-10-15)
        has_nm = "nm" in res and len(res["nm"][0]) == len(ch)
        ref = res["nm"][0] if has_nm else ch.mean()
        axs[1].plot(z, 1e3 * (ch - ref), "s-", color="#D55E00", ms=3, lw=0.8)
        axs[1].set_ylabel(r"$\Delta q$ (10$^{-3}$ e)")
        axs[1].set_ylim(-7.5, 6.5); axs[1].text(0.01, 0.88, "site charge, SDW $-$ " + ("NM" if has_nm else "mean") + " (period $\\lambda_{\\rm SDW}/2$: CDW)",
                    transform=axs[1].transAxes, fontsize=6)
        if fz is not None:
            axs[2].plot(z, 1e3 * fz * 13.605693 / 0.529177, "^-", color="#009E73", ms=3, lw=0.8)
            axs[2].set_ylabel(r"$F_z$ (meV/$\AA$)")
            axs[2].set_ylim(-75, 95); axs[2].text(0.01, 0.88, "force along $Q$ (drives the strain wave)", transform=axs[2].transAxes, fontsize=6)
        if "afm" in res and "nm" in res:
            de = lambda k: (res[k][3] - res["nm"][3]) / NL * 13.605693e3
            axs[0].set_title(f"E(AFM) $-$ E(NM) = {de('afm'):.1f}, E(SDW) $-$ E(NM) = {de('sdw'):.1f} meV/atom", fontsize=7)
    axs[2].set_xlabel(r"$z$ ($\AA$) along [001]")
    for ax in axs:
        ax.axhline(0, color="0.7", lw=0.5)
    fig.tight_layout(h_pad=0.4)
    fig.savefig(outfile, dpi=600); fig.savefig(Path(outfile).with_suffix(".pdf")); fig.savefig(Path(outfile).with_suffix(".svg"))
    print("wrote", outfile)


if __name__ == "__main__":
    if sys.argv[1] == "inputs":
        inputs(sys.argv[2])
    elif sys.argv[1] == "plot":
        plot(sys.argv[2], sys.argv[3])
