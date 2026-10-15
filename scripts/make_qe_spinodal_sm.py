#!/usr/bin/env python3
"""SM figure S_qe_spinodal: spinodal window implied by the complete QE mixing-enthalpy sets.

The QE enthalpies of mixing of the 40-atom SQS (PBE-relaxed and experimental lattice parameter, references FM Fe
and AFM Cr at the same convention) are fitted with Redlich-Kister polynomials
    dH_mix(x) = x (1 - x) sum_k L_k (1 - 2x)^k ,   k = 0..n-1,  n = 3 or 4 terms,
and combined with the ideal configurational entropy. The demixing curvature d2(dG_demix)/dx2 = -d2(dG_mix)/dx2
(Eq. demix convention of the paper) is positive where the alloy is spinodally unstable.

Why in the SM and not in Fig. 1b (decided 2026-10-15): the lower edge of the window depends on the fit order
(3 terms put Fe-36Cr on the edge at 748 K, 4 terms well inside), and the upper edge (0.96-0.99) lies far above the
CALPHAD window because these are 0 K ferromagnetic energies without paramagnetism. The CALPHAD window of
Xiong et al. (2011) at 748 K (0.30-0.71) remains the reference in the main text.

  make_qe_spinodal_sm.py  -> manuscript/figures/kdense_scripta/supplementary/S_qe_spinodal.{png,pdf,svg}, printed windows
"""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_figures as mf  # noqa: E402

KB = 8.617333e-5
OUT = HERE.parent / "manuscript" / "figures" / "kdense_scripta" / "supplementary" / "S_qe_spinodal"
CALPHAD_748 = (0.30, 0.71)          # Xiong et al. 2011, chemical spinodal at 748 K
X_FE36CR = 0.36


def rk_fit(x, dh, n):
    A = np.column_stack([x * (1 - x) * (1 - 2 * x) ** k for k in range(n)])
    L, *_ = np.linalg.lstsq(A, dh, rcond=None)
    rms = float(np.sqrt(np.mean((A @ L - dh) ** 2)))
    return L, rms


def rk_h(L, t):
    return t * (1 - t) * sum(L[k] * (1 - 2 * t) ** k for k in range(len(L)))


def rk_h2(L, t):
    """Exact second derivative of the Redlich-Kister enthalpy."""
    p = np.polynomial.Polynomial([0.0])
    for k, l in enumerate(L):
        p = p + l * np.polynomial.Polynomial([0, 1, -1]) * np.polynomial.Polynomial([1, -2]) ** k
    return p.deriv(2)(t)


def demix_curv(L, t, T):
    return -(rk_h2(L, t) + KB * T * (1 / t + 1 / (1 - t)))     # > 0: unstable


def windows(t, c):
    u = t[c > 0]
    if not len(u):
        return []
    br = np.where(np.diff(u) > 1.5 * (t[1] - t[0]))[0]
    return list(zip(np.r_[u[0], u[br + 1]], np.r_[u[br], u[-1]]))


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(HERE / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    sets = {"PBE-relaxed $a$": (mf.qe_mixing(), "o", "#D55E00"), "experimental $a$": (mf.qe_mixing(aexp=True), "D", "#0072B2")}
    t = np.linspace(0.002, 0.998, 4000)
    fig, axs = plt.subplots(1, 2, figsize=(7.0, 3.3))
    for lab, (q, mk, col) in sets.items():
        x, dh = np.asarray(q["x"]), np.asarray(q["dh"])
        m = (x > 0) & (x < 1)
        axs[0].plot(x[m], 1e3 * dh[m], mk, mfc="none", mec=col, ms=4.5, mew=1.0, label=f"QE, {lab}")
        for n, ls in ((3, "--"), (4, "-")):
            L, rms = rk_fit(x[m], dh[m], n)
            axs[0].plot(t, 1e3 * rk_h(L, t), ls, color=col, lw=0.9)
            c = demix_curv(L, t, 748.0)
            axs[1].plot(t, c, ls, color=col, lw=0.9, label=f"{lab}, {n} terms")
            w = [(round(float(a), 2), round(float(b), 2)) for a, b in windows(t, c)]
            print(f"{lab:18s} {n} terms  rms {1e3 * rms:.1f} meV/atom  unstable at 748 K: {w}")
    axs[0].set(xlabel=r"$x_{\rm Cr}$", ylabel=r"$\Delta H_{\rm mix}$ (meV/atom)", xlim=(0, 1))
    axs[0].axhline(0, color="0.7", lw=0.5)
    axs[0].legend(fontsize=6, loc="upper left", frameon=False)
    axs[0].text(0.98, 0.04, "dashed: 3 terms\nsolid: 4 terms", transform=axs[0].transAxes, ha="right", fontsize=6)
    axs[0].set_title("(a) QE enthalpy of mixing, Redlich$-$Kister fits", loc="left", fontsize=7)
    axs[1].axvspan(*CALPHAD_748, facecolor="0.85", edgecolor="0.6", lw=0.5, zorder=0, label="CALPHAD spinodal, 748 K")
    axs[1].axvline(X_FE36CR, color="0.3", lw=0.7, ls=":")
    axs[1].text(X_FE36CR + 0.01, 0.92, "Fe$-$36Cr", transform=axs[1].get_xaxis_transform(), fontsize=6)
    axs[1].axhline(0, color="0.5", lw=0.6)
    axs[1].set(xlabel=r"$x_{\rm Cr}$", ylabel=r"$\partial^2\Delta G_{\rm demix}/\partial x^2$ (eV/atom), 748 K",
               xlim=(0, 1), ylim=(-1.5, 1.2))
    axs[1].legend(fontsize=5.5, loc="upper center", bbox_to_anchor=(0.5, -0.2), frameon=False, ncol=2)   # below the axes: the curves fill the panel
    axs[1].set_title("(b) demixing curvature; $>0$: spinodally unstable", loc="left", fontsize=7)
    fig.tight_layout(w_pad=1.2)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT.with_suffix(f".{ext}"), dpi=600)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
