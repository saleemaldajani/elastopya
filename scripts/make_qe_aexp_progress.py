"""QE elastic series at the PBE-relaxed vs the room-temperature experimental lattice parameter, against EMTO-CPA:
f_SAW (6.4 um), K and G vs x_Cr, with the pure-element end points.

  cd scripts && python3 make_qe_aexp_progress.py   ->  manuscript/figures/kdense_scripta/supplementary/S_qe_aexp_progress.png
"""
import json, sys, numpy as np
from pathlib import Path
sys.path.insert(0, ".")
import make_figures as mf
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R = Path("../dft_qe/engaging/runs"); lit = mf.literature()
def series(suf):
    out = []
    for d in sorted(R.glob(f"sqs40_x*_k442_e60{suf}")) + (sorted(R.glob("sqs40_x0.80_k442_e60_aexp_crafm")) if suf == "_aexp" else []):
        if suf == "" and d.name.split("e60")[1]: continue
        j = d / "result.json"
        if j.exists():
            r = json.loads(j.read_text()); x = r["x_cr"]; c11, c12, c44 = r["C11"], r["C12"], r["C44"]
            b, cp = (c11 + 2 * c12) / 3, (c11 - c12) / 2
            g = ((2 * cp + 3 * c44) / 5 + 5 * (2 * cp) * c44 / (4 * c44 + 3 * (2 * cp))) / 2
            e = 9 * b * g / (3 * b + g)
            out.append((x, float(mf.rayleigh_f(b, e, x * mf.RHO_CR + (1 - x) * mf.RHO_FE)) * 1e3, b, g))
    pure = {"": ("Fe_fm_k16x16x16_e60", "Cr_afm_k16x16x16_e60"), "_aexp": ("Fe_fm_k16x16x16_e60_aexp", "Cr_afm_k16x16x16_e60_aexp")}[suf]
    for x, n in zip((0.0, 1.0), pure):                     # pure-element end points (2-atom cells, 16^3)
        j = R / n / "result.json"
        if j.exists():
            r = json.loads(j.read_text()); c11, c12, c44 = r["C11"], r["C12"], r["C44"]
            b, cp = (c11 + 2 * c12) / 3, (c11 - c12) / 2
            g = ((2 * cp + 3 * c44) / 5 + 5 * (2 * cp) * c44 / (4 * c44 + 3 * (2 * cp))) / 2
            e = 9 * b * g / (3 * b + g)
            out.append((x, float(mf.rayleigh_f(b, e, x * mf.RHO_CR + (1 - x) * mf.RHO_FE)) * 1e3, b, g))
    return np.array(sorted(out))
rel, aex = series(""), series("_aexp")
fig, ax = plt.subplots(1, 3, figsize=(7.4, 2.4))
ax[0].plot(lit["x"], np.array(lit["f"]) * 1e3, "o-", color="0.15", ms=3, lw=1, label="EMTO-CPA (Razumovskiy)")
ax[0].plot(rel[:, 0], rel[:, 1], "D--", color="#D55E00", mfc="white", ms=3.5, lw=0.8, label="QE, relaxed PBE volume")
ax[0].plot(aex[:, 0], aex[:, 1], "D-", color="#0072B2", ms=3.5, lw=1.1, label="QE, experimental lattice")
ax[0].set(xlabel="$x_{Cr}$", ylabel="$f_{SAW}$ (MHz, 6.4 µm)", xlim=(0, 1)); ax[0].legend(fontsize=5.5, frameon=False)
for k, (i, lab) in enumerate(((2, "$K$ (GPa)"), (3, "$G$ (GPa)")), start=1):
    ax[k].plot(rel[:, 0], rel[:, i], "D--", color="#D55E00", mfc="white", ms=3.5, lw=0.8)
    ax[k].plot(aex[:, 0], aex[:, i], "D-", color="#0072B2", ms=3.5, lw=1.1)
    ax[k].set(xlabel="$x_{Cr}$", ylabel=lab, xlim=(0, 1))
for a in ax: a.axvspan(0.30, 0.40, color="#FFB000", alpha=0.2, lw=0)
fig.tight_layout()
out = "../manuscript/figures/kdense_scripta/supplementary/S_qe_aexp_progress.png"
fig.savefig(out, dpi=300); print("wrote", out)
