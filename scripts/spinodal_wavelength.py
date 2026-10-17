#!/usr/bin/env python3
"""Spinodal wavelength of the decomposing MD structures from the Fourier transform of the Cr concentration field
(2026-10-17; the 3D analogue of the 2D FFT image analysis of Acta Sec. 3.8, Table 4).

Each 25^3 bcc cell (31,250 atoms, L = 71.65 A) is coarse-grained onto its 25^3 conventional cubic cells (two atoms
each), giving c(r) in {0, 1/2, 1}. S(k) = |FFT[c - <c>]|^2 / N is averaged over spherical shells of |k|. Two measures:
  k_peak  position of the maximum of S(k) for k > 0 (quantised in steps of 2 pi / L, so lambda_peak <= L)
  <k>     first moment of the excess S(k) - S_inf over 0 < k < K_CUT, where S_inf is the mean S(k) at large k (the
          uncorrelated, random-alloy level); lambda = 2 pi / <k> is continuous and is the quantity plotted
Real space (user suggestion, same date): on the same black/white map,
  chord   mean length of the Cr-rich and of the Fe-rich runs along every lattice line (x, y, z, periodic);
          lambda_chord = <l_Cr> + <l_Fe> is one Fe/Cr period and is not quantised like k (plotted quantity)
  C(r)    two-point correlation of the map, radially averaged (shown for Fe-36Cr; its first minimum sits near L/2,
          so it is box-limited)
The run compositions and phi values are those of lammps/decomposition/runs (12 x_Cr, 5 alloys, phi 0.5-1).

  spinodal_wavelength.py  -> data/spinodal_wavelength.csv, manuscript/figures/kdense_scripta/supplementary/S_spinodal_wavelength.*
"""
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RUNS = ROOT / "lammps" / "decomposition" / "runs"
OUT_CSV = ROOT / "data" / "spinodal_wavelength.csv"
OUT_FIG = ROOT / "manuscript" / "figures" / "kdense_scripta" / "supplementary" / "S_spinodal_wavelength"
K_CUT, K_INF = 4.0, 6.0          # nm^-1: first-moment window, and start of the large-k (uncorrelated) level
ACTA = dict(phi=[0.5, 0.6, 0.7, 0.8, 0.9, 1.0], lam=[4.19, 4.28, 4.45, 4.61, 4.70, 4.81], ss304=4.76)   # Acta Table 4


def field(path):
    t = path.read_text().split("\n")
    L = float(next(l for l in t if l.endswith("xlo xhi")).split()[1])
    i0 = next(i for i, l in enumerate(t) if l.startswith("Atoms")) + 2
    nat = int(re.search(r"(\d+) atoms", "\n".join(t[:5])).group(1))
    a = np.array([l.split()[:5] for l in t[i0:i0 + nat]], float)      # the Velocities section follows
    n = int(round(L / 2.866))
    idx = np.floor(a[:, 2:5] / (L / n) + 1e-6).astype(int) % n
    cr = (a[:, 1] == 2).astype(float)
    c = np.zeros((n, n, n)); cnt = np.zeros((n, n, n))
    np.add.at(c, tuple(idx.T), cr); np.add.at(cnt, tuple(idx.T), 1)
    return c / np.maximum(cnt, 1), L / 10      # Cr fraction per cubic cell, box length in nm


def binarize(c, x):
    """Black/white map as in the Acta image analysis (thresholded before the FFT): the cell field is smoothed over one
    cell (periodic Gaussian, sigma = 1 cell) and thresholded at the alloy composition: 1 = Cr-rich, 0 = Fe-rich."""
    n = c.shape[0]
    k = np.fft.fftfreq(n)
    kx, ky, kz = np.meshgrid(k, k, k, indexing="ij")
    g = np.exp(-2 * (np.pi ** 2) * (kx ** 2 + ky ** 2 + kz ** 2))      # Fourier transform of a sigma = 1 Gaussian
    sm = np.real(np.fft.ifftn(np.fft.fftn(c) * g))
    return (sm > x).astype(float)


def chords(b, cell):
    """Mean run lengths (nm) of the Cr-rich (1) and Fe-rich (0) phases along all lattice lines, periodic."""
    out = {0: [], 1: []}
    n = b.shape[0]
    for ax in range(3):
        for ln in np.moveaxis(b, ax, -1).reshape(-1, n):
            if ln.min() == ln.max():
                out[int(ln[0])].append(n); continue
            r = np.roll(ln, -int(np.argmax(ln != ln[0])))           # start at a boundary
            for seg in np.split(r, np.flatnonzero(np.diff(r) != 0) + 1):
                out[int(seg[0])].append(len(seg))
    return np.mean(out[1]) * cell, np.mean(out[0]) * cell


def correlation(b, Lnm):
    """Radially averaged two-point correlation of the map, normalised to 1 at r = 0, for r <= L/2."""
    n = b.shape[0]; d = b - b.mean()
    C = np.real(np.fft.ifftn(np.abs(np.fft.fftn(d)) ** 2)) / d.size
    C = (C / C.flat[0]).ravel()
    i = np.arange(n); i = np.minimum(i, n - i)
    r = np.sqrt(i[:, None, None] ** 2 + i[None, :, None] ** 2 + i[None, None, :] ** 2).ravel() * Lnm / n
    sh = np.rint(r / (Lnm / n)).astype(int); cnt = np.bincount(sh)
    rs = np.bincount(sh, r)[cnt > 0] / cnt[cnt > 0]; Cs = np.bincount(sh, C)[cnt > 0] / cnt[cnt > 0]
    m = rs <= Lnm / 2
    return rs[m], Cs[m]


def spectrum(c, Lnm):
    n = c.shape[0]
    F = np.fft.fftn(c - c.mean())
    S = (np.abs(F) ** 2 / c.size).ravel()
    k1 = 2 * np.pi * np.fft.fftfreq(n, d=Lnm / n)
    kx, ky, kz = np.meshgrid(k1, k1, k1, indexing="ij")
    k = np.sqrt(kx ** 2 + ky ** 2 + kz ** 2).ravel()
    dk = 2 * np.pi / Lnm
    shell = np.rint(k / dk).astype(int)
    m = shell > 0
    ks = np.bincount(shell[m], weights=k[m]) / np.maximum(np.bincount(shell[m]), 1)
    Ss = np.bincount(shell[m], weights=S[m]) / np.maximum(np.bincount(shell[m]), 1)
    good = np.bincount(shell[m]) > 0
    return ks[good], Ss[good]


def measures(ks, Ss):
    Sinf = Ss[ks > K_INF].mean()
    w = np.clip(Ss - Sinf, 0, None) * (ks < K_CUT)
    kmean = (ks * w).sum() / w.sum() if w.sum() > 0 else np.nan
    win = ks < K_CUT
    kpeak = ks[win][np.argmax(Ss[win])]
    return kmean, kpeak, Ss[ks < K_CUT].max() / Sinf


def main():
    rows = []
    for d in sorted(RUNS.glob("x*_L25_a*")):
        x = float(re.match(r"x([\d.]+)_", d.name).group(1)); alloy = int(d.name.split("_a")[1])
        for f in sorted(d.glob("phi*.data")):
            phi = float(f.stem[3:])
            c, L = field(f)
            km, kp, ratio = measures(*spectrum(c, L))
            bw = binarize(c, x)
            kb, kbp, _ = measures(*spectrum(bw, L))
            lc, lf = chords(bw, L / c.shape[0])
            rows.append((x, alloy, phi, L, km, 2 * np.pi / km, kp, 2 * np.pi / kp, ratio, kb, 2 * np.pi / kb, 2 * np.pi / kbp,
                         lc, lf, lc + lf))
    OUT_CSV.write_text("x_cr,alloy,phi,box_nm,k_mean_nm1,lambda_mean_nm,k_peak_nm1,lambda_peak_nm,peak_over_random,"
                       "k_mean_bw_nm1,lambda_mean_bw_nm,lambda_peak_bw_nm,chord_cr_nm,chord_fe_nm,lambda_chord_nm\n"
                       + "".join(",".join(f"{v:.5g}" for v in r) + "\n" for r in rows))
    print("wrote", OUT_CSV, len(rows), "structures")


def plot():
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(HERE))
    import make_figures_kdense as mk  # noqa: F401  (Scripta style)
    from make_figures_kdense import mf
    import make_response_figures as rf
    rf.apply_publication_style("scripta", palette_name="okabe_ito_on_white"); plt.rcParams["svg.fonttype"] = "none"
    d = np.genfromtxt(OUT_CSV, delimiter=",", names=True)
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(7.4, 2.6))
    cmap = rf._trunc("plasma", 0.0, 0.78)
    for phi in (0.5, 0.6, 0.7, 0.8, 0.9, 1.0):
        cc, L = field(RUNS / "x0.36_L25_a1" / f"phi{phi:.3f}.data")
        r, C = correlation(binarize(cc, 0.36), L)
        a.plot(r, C, "-", lw=0.9, color=cmap((phi - 0.5) / 0.5), label=rf"$\phi$ = {phi:.1f}")
    a.axhline(0, color="0.5", lw=0.5)
    a.set(xlabel=r"$r$ (nm)", ylabel=r"Two-point correlation $C(r)$", xlim=(0, L / 2), ylim=(-0.12, 1.02))
    a.text(0.97, 0.55, "Fe$-$36Cr, one alloy", transform=a.transAxes, ha="right", fontsize=6.5)
    a.legend(loc="upper right", ncol=2)
    win = d["phi"] > 0.55                                   # phi = 0.5 is the random alloy: no spinodal structure
    cols = {0.2: mf.BLUE, 0.36: mf.RED, 0.5: "0.2", 0.7: mf.GREEN}
    for x, col in cols.items():
        m = np.isclose(d["x_cr"], x) & win
        ph = np.unique(d["phi"][m])
        mu = [d["lambda_chord_nm"][m & np.isclose(d["phi"], p)].mean() for p in ph]
        sd = [d["lambda_chord_nm"][m & np.isclose(d["phi"], p)].std(ddof=1) for p in ph]
        b.errorbar(ph, mu, yerr=sd, fmt="-s", ms=3, lw=0.8, capsize=1.5, color=col, label=f"{round(100 * x)}% Cr, 25$^3$ (5 alloys)")
    b.plot(ACTA["phi"], ACTA["lam"], "o", mfc="none", mec="0.15", ms=4.5, mew=0.9, label="Fe$_{50}$Cr$_{50}$, Acta (2D FFT)")
    b.axhline(ACTA["ss304"], color="#8A5A00", lw=0.8, ls="--")
    b.text(0.49, ACTA["ss304"] + 0.06, "STEM-EDX, aged CF8 (Acta)", fontsize=6, color="#8A5A00", va="bottom")
    b.set(xlabel=r"Spinodal conversion fraction, $\phi$", ylabel=r"Wavelength $\lambda$ (nm), chord length",
          xlim=(0.47, 1.03), ylim=(3.0, 7.0))
    b.legend(loc="upper left")
    xs = np.unique(d["x_cr"])
    for key, mk_, col, lab in (("lambda_chord_nm", "s", mf.GREEN, r"chord length, $\phi$ = 1"),
                               ("lambda_mean_bw_nm", "^", "0.5", r"FFT first moment, $\phi$ = 1")):
        mu = [d[key][np.isclose(d["x_cr"], x) & np.isclose(d["phi"], 1.0)].mean() for x in xs]
        sd = [d[key][np.isclose(d["x_cr"], x) & np.isclose(d["phi"], 1.0)].std(ddof=1) for x in xs]
        c.errorbar(xs, mu, yerr=sd, fmt="-" + mk_, ms=3, lw=0.8, capsize=1.5, color=col, label=lab)
    c.axvspan(0.30, 0.71, color="0.9", lw=0, zorder=0, label="spinodal, 748 K")
    c.axhline(ACTA["ss304"], color="#8A5A00", lw=0.8, ls="--")
    c.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"Wavelength $\lambda$ (nm)", xlim=(0, 1), ylim=(3.0, 6.4))
    c.legend(loc="upper center", bbox_to_anchor=(0.5, 1.0))
    for ax, t in zip((a, b, c), "abc"):
        mk.style(ax, f"({t})")
    fig.set_size_inches(mk.PROFILE["widths_mm"]["full"] / 25.4, 2.6 * mk.PROFILE["widths_mm"]["full"] / 25.4 / 7.4)
    fig.tight_layout(w_pad=1.0)
    mk.normalize_text(fig)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT_FIG.with_suffix(f".{ext}"), dpi=600 if ext == "png" else None)
    print("wrote", OUT_FIG)


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "plot":
        plot()
    else:
        main()
