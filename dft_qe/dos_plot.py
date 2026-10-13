#!/usr/bin/env python3
"""Spin- and site-projected densities of states of the 40-atom SQS (projwfc.x output of dos_sqs.py), drawn as in
Korzhavyi et al., PRB 79, 054202 (2009), Fig. 4: Fe- and Cr-projected DOS per atom, spin up above and spin down
below the axis, versus E - E_F.

  dos_plot.py <outfile> <dos_dir>...      one panel per composition (dirs named x0.20, x0.95, ...)
"""
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
COL = {"Fe": "0.1", "Cr": "#009E73"}


def site_dos(d):
    """E - E_F (eV) and {element: (up, down)} per atom of that element (states/eV/atom)."""
    d = Path(d)
    ef = float(re.findall(r"the Fermi energy is\s+([-\d.]+) ev", (d / "nscf.out").read_text())[-1])
    acc, atoms, e = {}, {}, None
    prefix = "pdosT" if any(d.glob("pdosT.pdos_atm#*")) else "pdos"     # pdosT = tetrahedron projwfc (fixed)
    for f in d.glob(f"{prefix}.pdos_atm#*"):
        m = re.match(rf"{prefix}\.pdos_atm#(\d+)\((\w+)\)_wfc#\d+\(\w\)", f.name)
        if not m:
            continue
        idx, el = int(m.group(1)), m.group(2)
        a = np.loadtxt(f)
        e = a[:, 0]
        up, dw = acc.get(el, (0.0, 0.0))
        acc[el] = (up + a[:, 1], dw + a[:, 2])
        atoms.setdefault(el, set()).add(idx)
    return e - ef, {el: (u / len(atoms[el]), w / len(atoms[el])) for el, (u, w) in acc.items()}


def main(outfile, dirs):
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(HERE.parent / "scripts" / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta")
    fig, axs = plt.subplots(1, len(dirs), figsize=(2.0 * len(dirs) + 0.6, 3.0), squeeze=False, sharey=True)
    EWIN = (-7, 1.5)    # upper limit below the top of the computed empty bands (~+1.8 eV)
    ymax = 0.0
    for d in dirs:                                   # common y-range from the plotted energy window only
        e, dos = site_dos(d)
        w = (e >= EWIN[0]) & (e <= EWIN[1])
        ymax = max([ymax] + [float(max(u[w].max(), dw[w].max())) for u, dw in dos.values()])
    for ax, d in zip(axs[0], dirs):
        e, dos = site_dos(d)
        for el in ("Fe", "Cr"):
            if el in dos:
                up, dw = dos[el]
                ax.plot(e, up, color=COL[el], lw=0.9, label=el)
                ax.plot(e, -dw, color=COL[el], lw=0.9)
        ax.axhline(0, color="0.5", lw=0.5); ax.axvline(0, color="0.5", lw=0.6, ls=":")
        ax.set(xlim=EWIN, ylim=(-1.08 * ymax, 1.08 * ymax), xlabel=r"$E - E_{\rm F}$ (eV)")
        x = float(Path(d).name[1:])
        ax.set_title(f"{round(100 * x)}% Cr", fontsize=8)
        ax.text(0.03, 0.97, "spin up", transform=ax.transAxes, va="top", fontsize=7, color="0.35")
        ax.text(0.03, 0.03, "spin down", transform=ax.transAxes, va="bottom", fontsize=7, color="0.35")
    axs[0, 0].set_ylabel("DOS (states/eV/atom)")
    axs[0, -1].legend(loc="upper right", fontsize=7)
    fig.tight_layout(w_pad=0.6)
    fig.savefig(outfile, dpi=600); fig.savefig(Path(outfile).with_suffix(".pdf"))
    print("wrote", outfile)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
