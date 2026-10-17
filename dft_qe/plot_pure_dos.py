#!/usr/bin/env python3
"""SM figure S_dos_pure (2026-10-16): spin-resolved total DOS of the pure metals from engaging/pure_dos/ (pure_dos.sbatch:
nscf with the optimised tetrahedron method, 36^3 mesh for the bcc cells, 24^3 for the 2-atom AFM Cr cell), at the
PBE-relaxed (solid) and experimental (dashed) lattice parameter. N(E_F) per atom is compared with the measured
electronic specific-heat coefficient, gamma = (pi^2/3) k_B^2 N(E_F) (1 + lambda): gamma/2.357 mJ mol^-1 K^-2 gives
N(E_F)(1 + lambda) in states/eV/atom (Fe 4.98, Cr 1.40 mJ mol^-1 K^-2, Kittel, Table 6.3).

  plot_pure_dos.py  -> manuscript/figures/kdense_scripta/supplementary/S_dos_pure.{png,pdf,svg}; prints N(E_F)
"""
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ENG = HERE / "engaging" / "pure_dos"
OUT = HERE.parent / "manuscript" / "figures" / "kdense_scripta" / "supplementary" / "S_dos_pure"
GAMMA = {"Fe": 4.98, "Cr": 1.40}          # mJ mol^-1 K^-2 (Kittel)
CASES = (("Fe_fm", "Fe, ferromagnetic", 1), ("Cr_nm", "Cr, non-magnetic", 1), ("Cr_afm", "Cr, commensurate AFM", 2))


def load(d, nat):
    t = (d / "dos.dat").read_text()
    ef = float(re.search(r"EFermi =\s*([-\d.]+)", t).group(1))
    a = np.loadtxt(d / "dos.dat")
    if a.shape[1] == 3:                       # nspin = 1 (non-magnetic Cr): E, dos (both spins), integrated dos
        return a[:, 0] - ef, a[:, 1] / nat / 2, a[:, 1] / nat / 2
    return a[:, 0] - ef, a[:, 1] / nat, a[:, 2] / nat


def nef(e, up, dw, w=0.05):
    m = np.abs(e) <= w
    return float(up[m].mean()), float(dw[m].mean())


def main():
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(HERE.parent / "scripts" / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.7), sharey=True)
    rows = []
    for ax, (name, title, nat) in zip(axs, CASES):
        for suf, ls, lab in (("", "-", "PBE-relaxed $a$"), ("_aexp", "--", "experimental $a$")):
            d = ENG / f"{name}{suf}"
            if not (d / "dos.dat").exists():
                continue
            e, up, dw = load(d, nat)
            ax.plot(e, up, ls, color="#C02F00", lw=0.9, label=f"spin up, {lab}" if name == "Fe_fm" else None)
            ax.plot(e, -dw, ls, color="#0072B2", lw=0.9, label=f"spin down, {lab}" if name == "Fe_fm" else None)
            nu, nd = nef(e, up, dw)
            el = name[:2]
            rows.append((name + suf, nu, nd, nu + nd, GAMMA[el] / 2.357 / (nu + nd)))
        if name == "Fe_fm":   # Korzhavyi et al. 2009, Fig. 4: Fe site in Fe98Cr02 (EMTO-CPA, LDA), digitized, per Ry -> per eV
            import csv
            rows_k = [r for r in csv.reader(l for l in open(HERE / "refdata" / "korzhavyi2009_fig4_pdos_fm.csv") if not l.startswith("#"))]
            hdr = rows_k[0]; ix = {h: i for i, h in enumerate(hdr)}
            for spin, sg in (("up", 1), ("down", -1)):
                pts = np.array([(float(r[ix["E_minus_EF_Ry"]]), float(r[ix["DOS_states_per_Ry"]])) for r in rows_k[1:]
                                if r[ix["series"]].startswith("Fe(98%)") and r[ix["spin"]] == spin])
                if len(pts):
                    ax.plot(pts[:, 0] * 13.605693, sg * pts[:, 1] / 13.605693, ".", color="0.35", ms=1.2,
                            label="Fe site in Fe$_{98}$Cr$_2$, EMTO-CPA [Korzhavyi]" if spin == "up" else None, zorder=0)
        ax.axvline(0, color="0.5", lw=0.6, ls=":"); ax.axhline(0, color="0.5", lw=0.5)
        ax.set(xlim=(-8, 4), xlabel=r"$E-E_{\rm F}$ (eV)")
        ax.set_title(title, loc="left", fontsize=7)
    axs[0].set(ylim=(-3.2, 3.2), ylabel="DOS (states/eV/atom/spin)")
    axs[0].legend(fontsize=5.5, loc="lower left", frameon=False)
    for ax, t in zip(axs, "abc"):
        ax.set_title(f"({t}) " + ax.get_title(loc="left"), loc="left", fontsize=7)
    fig.tight_layout(w_pad=1.6)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT.with_suffix(f".{ext}"), dpi=600)
    print("wrote", OUT)
    print("case, N_up(E_F), N_dw(E_F), N(E_F) states/eV/atom, implied 1+lambda = gamma_exp/(2.357 N)")
    for r in rows:
        print(f"{r[0]:14s} {r[1]:.3f} {r[2]:.3f} {r[3]:.3f} {r[4]:.2f}")


if __name__ == "__main__":
    main()
