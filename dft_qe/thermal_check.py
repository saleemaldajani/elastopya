#!/usr/bin/env python3
"""Validation of the thermal-transport estimates of SM Section S1 with every elastic tensor of this work.

Section S1 estimates, from the elastic moduli alone,
  Debye temperature   Theta_D = (h/k_B) (3n/4 pi)^(1/3) v_m,  v_m^-3 = (2 v_t^-3 + v_l^-3)/3        (Anderson 1963)
  minimum conductivity  kappa_min = 0.87 k_B (N_A m/M)^(2/3) rho^(1/6) E^(1/2),  m = 1 atom per formula unit (Clarke 2003)
  thermal diffusivity  D = kappa / (rho c_p),  c_p = 460.5 J/(kg K) (as in the working sheet)
These are applied here to EMTO-CPA, VASP, QE (both lattice parameters) and MD (random and decomposed), and checked against
  * the calorimetric Debye temperatures of Fe (470 K) and Cr (630 K) (Kittel), and Debye temperatures from measured or
    computed phonon DOS (second moment, k_B Theta_D = hbar sqrt(5<w^2>/3)): Fultz 1995 (0-100 % Cr), NRIXS Fe, ours;
  * the TGS thermal diffusivity of CF8 (Fig. S7d): 3.50e-6 m^2/s unaged, 3.64 / 4.19 / 4.30 / 4.53e-6 after aging at
    290 / 330 / 360 / 400 C.
Note: the working sheet (and hence Fig. S7e,f of the original draft) evaluates kappa with Planck's h x 1e11 in place of
k_B and m = 2. That multiplies kappa and D by 7.62 relative to the equation printed in S1. Relative trends are unaffected.

  thermal_check.py   -> prints the comparison, writes S_thermal.{png,pdf,svg} and thermal_results.json
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "scripts"))
import make_figures as mf  # noqa: E402

H, KB, HBAR, NA = 6.62607015e-34, 1.380649e-23, 1.054571817e-34, 6.02214076e23
CP = 460.548
CM1_HZ = 2.99792458e10
CALORIMETRIC = {0.0: 470.0, 1.0: 630.0}
CF8 = [(25, 3.50), (290, 3.64), (330, 4.19), (360, 4.30), (400, 4.53)]      # Fig. S7d (aging temperature in C, 1e-6 m^2/s)
SHEET_FACTOR = (H * 1e11 / KB) * 2 ** (2 / 3)


def vrh(c11, c12, c44):
    b = (c11 + 2 * c12) / 3; cp = (c11 - c12) / 2
    g = ((2 * cp + 3 * c44) / 5 + 5 * (2 * cp) * c44 / (4 * c44 + 3 * (2 * cp))) / 2
    return b, g, 9 * b * g / (3 * b + g)


def props(c11, c12, c44, x, rho=None):
    b, g, e = vrh(c11, c12, c44)
    rho = rho or (x * mf.RHO_CR + (1 - x) * mf.RHO_FE)
    M = ((1 - x) * 55.845 + x * 51.996) * 1e-3
    n = rho * NA / M
    vt, vl = np.sqrt(g * 1e9 / rho), np.sqrt((b + 4 * g / 3) * 1e9 / rho)
    vm = (((2 / vt ** 3) + (1 / vl ** 3)) / 3) ** (-1 / 3)
    theta = H / KB * (3 * n / (4 * np.pi)) ** (1 / 3) * vm
    kappa = 0.87 * KB * (NA / M) ** (2 / 3) * rho ** (1 / 6) * np.sqrt(e * 1e9)
    return dict(x=x, E=e, G=g, B=b, theta=theta, kappa=kappa, D=kappa / (rho * CP) * 1e6)


def theta_dos(nu_thz, g):
    nu, g = np.asarray(nu_thz, float), np.clip(np.asarray(g, float), 0, None)
    m = nu > 0
    w2 = np.trapezoid(g[m] * (2 * np.pi * nu[m] * 1e12) ** 2, nu[m]) / np.trapezoid(g[m], nu[m])
    return float(HBAR * np.sqrt(5 * w2 / 3) / KB)


def rows(p):
    return list(csv.DictReader(l for l in open(p) if not l.startswith("#")))


def tensor_sets():
    sets = {}
    lit = rows(ROOT / "data" / "sheet_raw_razumovskiy.csv")
    sets["EMTO-CPA"] = [props(float(r["bulk_gpa"]) + 4 / 3 * float(r["cprime_gpa"]), float(r["bulk_gpa"]) - 2 / 3 * float(r["cprime_gpa"]),
                              float(r["c44_gpa"]), min(float(r["x_cr"]), 1.0)) for r in lit]
    sets["VASP"] = [props(c11, c12, c44, x) for t, (x, c11, c12, c44) in mf.saw_tensors().items() if t.startswith("dft_")]
    for suf, lab in (("", "QE relaxed"), ("_aexp", "QE exp. a")):
        pts = []
        for el, mag, x in (("Fe", "fm", 0.0), ("Cr", "afm", 1.0)):
            j = HERE / "engaging" / "runs" / f"{el}_{mag}_k16x16x16_e60{suf}" / "result.json"
            if j.exists():
                r = json.loads(j.read_text()); pts.append(props(r["C11"], r["C12"], r["C44"], x))
        for d in sorted((HERE / "engaging" / "runs").glob(f"sqs40_x*_k442_e60{suf}")):
            if suf == "" and d.name.split("e60")[1]:
                continue
            if (d / "result.json").exists():
                r = json.loads((d / "result.json").read_text()); pts.append(props(r["C11"], r["C12"], r["C44"], r["x_cr"]))
        sets[lab] = sorted(pts, key=lambda p: p["x"])
    md = rows(ROOT / "data" / "md_random_alloys.csv")
    sets["MD random"] = [props(float(r["c11"]), float(r["c12"]), float(r["c44"]), float(r["x_cr"])) for r in md]
    return {k: sorted(v, key=lambda p: p["x"]) for k, v in sets.items()}


def decomposition_change():
    """Relative change of E and D on full decomposition (phi 0.5 -> 1) of each MD composition (S1 scaling, fixed rho)."""
    out = {}
    dd = rows(ROOT / "data" / "md_decomposition.csv")
    for x0 in sorted({float(r["x_cr"]) for r in dd}):
        sel = [(float(r["phi"]), float(r["c11"]), float(r["c12"]), float(r["c44"])) for r in dd
               if abs(float(r["x_cr"]) - x0) < 1e-9 and r["L"] == "25"]
        mean = lambda ph: np.mean([s_[1:] for s_ in sel if abs(s_[0] - ph) < 1e-6], axis=0)   # over MD alloys (seeds)
        if {round(s_[0], 3) for s_ in sel} >= {0.5, 1.0}:
            a, b = props(*mean(0.5), x0), props(*mean(1.0), x0)
            out[x0] = dict(dE_pct=100 * (b["E"] / a["E"] - 1), dD_pct=100 * (b["D"] / a["D"] - 1), dG_pct=100 * (b["G"] / a["G"] - 1))
    return out


def main():
    sets = tensor_sets()
    res = {"sets": sets, "sheet_factor": SHEET_FACTOR, "dos_theta": {}, "decomposition": decomposition_change()}
    fu = rows(HERE / "refdata" / "fecr_phonon_dos_fultz1995.csv"); nu = np.array([float(r["nu_THz"]) for r in fu])
    for col, x in (("Fe", 0.0), ("Fe70Cr30", 0.30), ("Fe53Cr47", 0.47), ("Fe30Cr70", 0.70), ("Cr", 1.0)):
        res["dos_theta"].setdefault("Fultz 1995", []).append((x, theta_dos(nu, [float(r[col]) for r in fu])))
    mau = rows(HERE / "refdata" / "fe_phonon_dos_mauger2014.csv")
    for T in ("30", "300"):
        d = np.array([[float(r["nu_THz"]), float(r["g_perTHz"])] for r in mau if r["T_K"] == T])
        res["dos_theta"].setdefault("NRIXS", []).append((0.0, theta_dos(d[:, 0], d[:, 1])))
    for lab, d, x in (("QE Fe (finite displ., exp. a)", HERE / "phonopy_derived" / "Fe_cub_aexp", 0.0),
                      ("QE Cr (DFPT AFM, exp. a)", HERE / "engaging" / "phonon" / "Cr_afm_aexp", 1.0)):
        if (d / "dos.dat").exists():
            a = np.loadtxt(d / "dos.dat"); res["dos_theta"][lab] = [(x, theta_dos(a[:, 0] * CM1_HZ * 1e-12, a[:, 1]))]
    print(f"sheet / S1-as-written factor for kappa and D: {SHEET_FACTOR:.2f}")
    for k, pts in sets.items():
        for p in pts:
            if p["x"] in (0.0, 1.0) or abs(p["x"] - 0.36) < 0.02:
                print(f"{k:12s} x={p['x']:.2f}  Theta_D={p['theta']:.0f} K  kappa_min={p['kappa']:.2f} W/mK  D={p['D']:.3f}e-6 m2/s")
    print("Theta_D from phonon DOS:", {k: [(x, round(t)) for x, t in v] for k, v in res["dos_theta"].items()})
    print("decomposition (MD, S1 scaling):", {k: {q: round(w, 2) for q, w in v.items()} for k, v in res["decomposition"].items()})
    m = [(T, D / CF8[0][1] - 1) for T, D in CF8[1:]]
    print("CF8 measured dD/D on aging:", [(T, round(100 * v, 1)) for T, v in m])
    res["cf8_rel"] = m
    (HERE / "thermal_results.json").write_text(json.dumps(res, indent=1, default=float))
    return res


def plot(res):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ROOT / "scripts" / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    st = {"EMTO-CPA": ("o-", "0.15", "0.15"), "VASP": ("o", "#D55E00", "#D55E00"), "QE relaxed": ("D--", "#8F3A00", "white"),
          "QE exp. a": ("D-", "#0072B2", "#0072B2"), "MD random": ("s-", "#009E73", "#009E73")}
    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.7))
    for lab, pts in res["sets"].items():
        x = [p["x"] for p in pts]; mk, col, fc = st[lab]
        axs[0].plot(x, [p["theta"] for p in pts], mk, color=col, mfc=fc, ms=3, lw=0.8, label=f"elastic, {lab}")
        axs[1].plot(x, [p["D"] for p in pts], mk, color=col, mfc=fc, ms=3, lw=0.8, label=lab)
    for lab, mk in (("Fultz 1995", "^"), ("NRIXS", "v")):
        p = np.array(res["dos_theta"][lab]); axs[0].plot(p[:, 0], p[:, 1], mk, color="#CC79A7", ms=5, label=f"phonon DOS, {lab}")
    for lab in [k for k in res["dos_theta"] if k.startswith("QE")]:
        p = np.array(res["dos_theta"][lab]); axs[0].plot(p[:, 0], p[:, 1], "*", color="#0072B2", ms=7, label="phonon DOS, QE" if "Fe" in lab else None)
    c = np.array(sorted(CALORIMETRIC.items())); axs[0].plot(c[:, 0], c[:, 1], "X", color="#C00000", ms=7, label="calorimetric")
    axs[0].set(xlabel="$x_{Cr}$", ylabel=r"$\Theta_D$ (K)", xlim=(-0.03, 1.03)); axs[0].set_title("(a) Debye temperature", loc="left", fontsize=7)
    axs[0].legend(fontsize=4.4, loc="upper left", frameon=False, ncol=1)
    axs[1].set(xlabel="$x_{Cr}$", ylabel=r"$D_{\min}$ ($10^{-6}$ m$^2$/s)", xlim=(-0.03, 1.03))
    axs[1].set_title(r"(b) S1 estimate, $\kappa_{\min}/\rho c_p$", loc="left", fontsize=7); axs[1].legend(fontsize=4.6, loc="upper left", frameon=False)
    # (c) relative change: measured on aging (CF8) vs S1 scaling on decomposition (MD)
    T = [t for t, _ in res["cf8_rel"]]; v = [100 * w for _, w in res["cf8_rel"]]
    axs[2].plot(T, v, "o", color="0.1", mfc="white", ms=5, label="measured, CF8 (TGS)")
    dec = res["decomposition"]
    for x0, col in ((0.36, "#009E73"), (0.5, "#0B5D3B")):
        if x0 in dec:
            axs[2].axhline(dec[x0]["dD_pct"], color=col, lw=1.0, ls="--", label=f"S1 scaling, MD Fe$-${round(100 * x0)}Cr decomposed")
    axs[2].set(xlabel="aging temperature (°C)", ylabel=r"$\Delta D/D$ (%)", xlim=(250, 420))
    axs[2].set_title("(c) change on aging / decomposition", loc="left", fontsize=7); axs[2].legend(fontsize=4.8, loc="upper left", frameon=False)
    fig.tight_layout()
    out = ROOT / "manuscript" / "figures" / "kdense_scripta" / "supplementary" / "S_thermal"
    fig.savefig(f"{out}.png", dpi=600); fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.svg"); print("wrote", out)


if __name__ == "__main__":
    plot(main())
