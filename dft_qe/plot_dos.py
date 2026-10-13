#!/usr/bin/env python3
"""Spin- and atom-projected DOS of the 40-atom SQS (projwfc.x output of dos_sqs.py), in the layout of
Korzhavyi et al., PRB 79, 054202 (2009), Fig. 4: Fe- and Cr-projected DOS per atom, spin up above and
spin down below zero, energies relative to E_F.

  plot_dos.py <outfile> <dos_dir>...      e.g. dft_qe/engaging/dos/x0.20 dft_qe/engaging/dos/x0.95
"""
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "scripts" / "kdense_scripta" / "scripts"))
COL = {"Fe": "0.1", "Cr": "#009E73"}


def element_dos(d):
    """Per-atom projected DOS by element: {el: (E - E_F, up, dw)} in states/eV/atom."""
    d = Path(d)
    ef = None
    for out in ("nscf.out", "scf.out"):          # tetrahedra nscf does not print E_F; fall back to the SCF value
        m = re.findall(r"the Fermi energy is\s+([-\d.]+) ev", (d / out).read_text())
        if m:
            ef = float(m[-1]); break
    acc, atoms = {}, defaultdict(set)
    for f in d.glob("pdos.pdos_atm#*"):
        m = re.match(r"pdos\.pdos_atm#(\d+)\((\w+)\)_wfc", f.name)
        n, el = int(m.group(1)), m.group(2)
        a = np.loadtxt(f)
        atoms[el].add(n)
        if el not in acc:
            acc[el] = [a[:, 0], np.zeros(len(a)), np.zeros(len(a))]
        acc[el][1] += a[:, 1]; acc[el][2] += a[:, 2]
    return ({el: (e - ef, up / len(atoms[el]), dw / len(atoms[el])) for el, (e, up, dw) in acc.items()},
            {el: len(v) for el, v in atoms.items()})


def plot(outfile, dirs):
    import matplotlib.pyplot as plt
    from style_presets import apply_publication_style
    apply_publication_style("scripta")
    fig, axs = plt.subplots(1, len(dirs), figsize=(2.3 * len(dirs), 3.0), squeeze=False, sharey=True)
    for ax, d in zip(axs[0], dirs):
        dos, nat = element_dos(d)
        x = nat.get("Cr", 0) / sum(nat.values())
        for el in ("Fe", "Cr"):
            if el in dos:
                e, up, dw = dos[el]
                ax.plot(e, up, color=COL[el], lw=0.8, label=el)
                ax.plot(e, -dw, color=COL[el], lw=0.8)
        ax.axvline(0, color="0.5", lw=0.5, ls=":"); ax.axhline(0, color="0.5", lw=0.5)
        ax.set(xlim=(-7, 3), xlabel=r"$E - E_{\rm F}$ (eV)",
               title=f"Fe$_{{{100 - round(100 * x)}}}$Cr$_{{{round(100 * x)}}}$")
        ax.text(0.03, 0.97, "spin up", transform=ax.transAxes, va="top", fontsize=7)
        ax.text(0.03, 0.03, "spin down", transform=ax.transAxes, va="bottom", fontsize=7)
    axs[0, 0].set_ylabel("PDOS (states/eV/atom)")
    axs[0, -1].legend(loc="upper right", fontsize=7)
    fig.tight_layout()
    fig.savefig(outfile, dpi=600); fig.savefig(Path(outfile).with_suffix(".pdf"))
    print("wrote", outfile)


if __name__ == "__main__":
    plot(sys.argv[1], sys.argv[2:])
