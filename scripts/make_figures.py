#!/usr/bin/env python3
"""Regenerate Figs. 1-3 of the Scripta manuscript from the working Google Sheet.

Usage
  python scripts/make_figures.py --xlsx <sheet export .xlsx>   # refresh data/ snapshots, then plot
  python scripts/make_figures.py                               # plot from data/ snapshots

The figures reproduce the panels and series of the draft reviewed by M. P. Short. The
calculations use the sheet's own formulas ('Analyzed data' row 2):
  C11 = (3B + 4C')/3,  C12 = C11 - 2C',  G_V, G_R, G (Hill),  E = 9BG/(3B+G),
  nu = (3B - E)/(6B),  rho = [x rho_Cr + (1-x) rho_Fe],
  v_R = (0.874 + 0.196 nu - 0.042 nu^2 - 0.055 nu^3) sqrt(E / (2(1+nu) rho)),
  f_SAW = v_R / 6.4 um.
Three spreadsheet referencing errors are corrected (CORRECTION notes below); all other
numbers are used as they appear in the sheet.
"""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
FIG = ROOT / "manuscript" / "figures"
LAMBDA = 6.4e-6                       # m, TGS grating period (Acta Sec. 3.5)
RHO_FE, RHO_CR = 7870.0, 7200.0       # kg/m^3, 'Analyzed data'!B38:B39
X_FE36CR = 0.36
DFT_ROWS = (2, 12, 17, 22, 23)        # rows of 'Analyzed data' holding the author VASP values
MD_ROWS = range(28, 34)               # rows of 'Analyzed data' holding the Acta EAM moduli
KB, T = 8.61733031e-5, 298.15

# Google-Sheets palette of the original figures
BLUE, RED, YELLOW, GREEN, ORANGE = "#4285F4", "#EA4335", "#F4B400", "#34A853", "#FF6D01"
PURPLE = "#8E44AD"


# ----------------------------------------------------------------------------- snapshot
def write_csv(name, header, rows):
    with open(DATA / name, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def snapshot(xlsx):
    import openpyxl
    wb = openpyxl.load_workbook(xlsx, data_only=True)
    raw, ad, ed, sp = wb["Raw data"], wb["Analyzed data"], wb["Enthalpy of demixing"], wb["Spinodal"]
    write_csv("sheet_raw_razumovskiy.csv", ["x_cr", "cprime_gpa", "c44_gpa", "bulk_gpa"],
              [[raw.cell(r, 1).value, raw.cell(r, 2).value, raw.cell(r, 5).value, raw.cell(r, 8).value]
               for r in range(2, 23)])
    write_csv("sheet_dft.csv", ["x_cr", "c11_gpa", "c12_gpa", "c44_gpa", "fsaw_6p4um_ghz"],
              [[ad.cell(r, 1).value, ad.cell(r, 4).value, ad.cell(r, 7).value, ad.cell(r, 10).value,
                ad.cell(r, 54).value] for r in DFT_ROWS])
    write_csv("sheet_md.csv", ["phi", "bulk_gpa", "young_gpa", "c11_phi1_gpa"],
              [[ad.cell(r, 1).value, ad.cell(r, 11).value, ad.cell(r, 23).value, ad["C2"].value]
               for r in MD_ROWS])
    write_csv("sheet_band.csv", ["unallowed_min_ghz", "unallowed_max_ghz"], [[ad["DY2"].value, ad["DZ2"].value]])
    write_csv("sheet_demixing.csv",
              ["x_cr_a", "muffin_tin", "gga", "experimental", "x_cr_b", "pwscf_paw", "siesta_nc",
               "sheet_second_derivative_selected"],
              [[ed.cell(r, c).value for c in (7, 8, 9, 10, 12, 13, 14, 16)] for r in range(2, 35)])
    write_csv("sheet_levesque_parameters.csv", ["series", "alpha", "beta", "gamma", "delta"],
              [["PWSCF-PAW"] + [sp.cell(r, 2).value for r in range(2, 6)],
               ["SIESTA-NC"] + [sp.cell(r, 3).value for r in range(2, 6)]])


def read_csv(name):
    with open(DATA / name) as fh:
        rows = list(csv.DictReader(fh))
    out = {}
    for k in rows[0]:
        try:
            out[k] = np.array([np.nan if r[k] in ("", "None") else float(r[k]) for r in rows])
        except ValueError:
            out[k] = [r[k] for r in rows]
    return out


# ----------------------------------------------------------------------------- physics (sheet formulas)
def rayleigh_f(bulk, young, rho):
    nu = (3 * bulk - young) / (6 * bulk)
    vr = (0.874 + 0.196 * nu - 0.042 * nu**2 - 0.055 * nu**3) * np.sqrt(young * 1e9 / (2 * (1 + nu) * rho))
    return vr / LAMBDA / 1e9


def literature():
    """CORRECTION 1: the sheet's C11 and B columns read 'Raw data'!H one row down from x~0.11
    (and wrap at x=1). Here every quantity is taken from the same 'Raw data' row."""
    d = read_csv("sheet_raw_razumovskiy.csv")
    x, cp, c44, b = d["x_cr"], d["cprime_gpa"], d["c44_gpa"], d["bulk_gpa"]
    c11 = (3 * b + 4 * cp) / 3
    gv = (2 * cp + 3 * c44) / 5
    gr = 5 * (2 * cp) * c44 / (4 * c44 + 3 * (2 * cp))
    g = (gv + gr) / 2
    e = 9 * b * g / (3 * b + g)
    rho = x * RHO_CR + (1 - x) * RHO_FE
    return dict(x=x, c11=c11, f=rayleigh_f(b, e, rho))


def md():
    """CORRECTION 2: density held at the equiatomic value (the sheet used phi as x_Cr).
    CORRECTION 3: Poisson ratio uses each row's own bulk modulus (the sheet used $K$23)."""
    d = read_csv("sheet_md.csv")
    rho = 0.5 * RHO_CR + 0.5 * RHO_FE
    return dict(phi=d["phi"], f=rayleigh_f(d["bulk_gpa"], d["young_gpa"], rho), c11_phi1=d["c11_phi1_gpa"][0])


def md_random():
    """Random (homogeneous) alloys, LAMMPS + Bonny 2011 EAM (scripts/collect_md_random.py).
    Same cubic VRH and Rayleigh formulas, and same density interpolation, as the literature series."""
    path = DATA / "md_random_alloys.csv"
    if not path.exists():
        return None
    d = read_csv(path.name)
    x, c11, c12, c44 = d["x_cr"], d["c11"], d["c12"], d["c44"]
    b = (c11 + 2 * c12) / 3
    gv = (c11 - c12 + 3 * c44) / 5
    gr = 5 * (c11 - c12) * c44 / (4 * c44 + 3 * (c11 - c12))
    g = (gv + gr) / 2
    e = 9 * b * g / (3 * b + g)
    return dict(x=x, c11=c11, c11_sd=d["c11_sd"], f=rayleigh_f(b, e, x * RHO_CR + (1 - x) * RHO_FE))


def md_decomposition():
    """Small-cell (25^3) decomposition series from scripts/collect_md_decomposition.py.
    Returns {x_cr: dict(phi, f_6p4)} using the same cubic VRH/Rayleigh formulas and density
    interpolation as the random-alloy series."""
    path = DATA / "md_decomposition.csv"
    if not path.exists():
        return {}
    d = read_csv(path.name)
    out = {}
    for x in sorted(set(d["x_cr"])):
        m = d["x_cr"] == x
        phi, c11, c12, c44 = d["phi"][m], d["c11"][m], d["c12"][m], d["c44"][m]
        b = (c11 + 2 * c12) / 3
        gv = (c11 - c12 + 3 * c44) / 5
        gr = 5 * (c11 - c12) * c44 / (4 * c44 + 3 * (c11 - c12))
        g = (gv + gr) / 2
        e = 9 * b * g / (3 * b + g)
        o = np.argsort(phi)
        out[float(x)] = dict(phi=phi[o], f=rayleigh_f(b, e, x * RHO_CR + (1 - x) * RHO_FE)[o])
    return out


def levesque_model(x, a, b, c, d):
    """dG_demix = (x-a)(bx^2+cx+d) x(1-x) + kT[-x ln x - (1-x) ln(1-x)] and its exact curvature."""
    h = np.poly1d([b, c, d]) * np.poly1d([1.0, -a]) * np.poly1d([-1.0, 1.0, 0.0])
    xc = np.clip(x, 1e-12, 1 - 1e-12)
    g = h(x) - KB * T * (xc * np.log(xc) + (1 - xc) * np.log(1 - xc))
    return g, np.polyder(h, 2)(xc) - KB * T * (1 / xc + 1 / (1 - xc))


def numerical_curvature(x, g, targets):
    """Three-point second derivative on the sheet's non-uniform grid (unique x only)."""
    keep = ~np.isnan(x)
    xu, idx = np.unique(x[keep], return_index=True)
    gu = g[keep][idx]
    out = []
    for t in targets:
        i = int(np.argmin(abs(xu - t)))
        if i == 0 or i == len(xu) - 1:
            out.append(np.nan)
            continue
        h1, h2 = xu[i] - xu[i - 1], xu[i + 1] - xu[i]
        out.append(2 * ((gu[i + 1] - gu[i]) / h2 - (gu[i] - gu[i - 1]) / h1) / (h1 + h2))
    return np.array(out)


# ----------------------------------------------------------------------------- figures
def style(ax, label):
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True, labelsize=7)
    ax.tick_params(which="major", length=3.5)
    ax.tick_params(which="minor", length=2)
    ax.set_title(label, loc="left", fontsize=9, fontweight="bold", pad=4)


def save(fig, stem):
    fig.tight_layout(w_pad=1.6)
    fig.savefig(FIG / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIG / f"{stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fig1():
    d = read_csv("sheet_demixing.csv")
    p = read_csv("sheet_levesque_parameters.csv")
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.0))
    m = dict(ms=3.6, ls="none")
    a.plot(d["x_cr_a"], d["muffin_tin"], "o", color=BLUE, label="Muffin-tin approximation", **m)
    a.plot(d["x_cr_a"], d["gga"], "o", color=RED, label="GGA (this work)", **m)
    a.plot(d["x_cr_a"], d["experimental"], "o", color=YELLOW, label="Experimental", **m)
    a.plot(d["x_cr_b"], d["pwscf_paw"], "o", color=GREEN, label="PWSCF-PAW", **m)
    a.plot(d["x_cr_b"], d["siesta_nc"], "o", color=ORANGE, label="SIESTA-NC", **m)
    a.axhline(0, color="0.5", lw=0.7)
    a.set(xlabel=r"$x_{\rm Cr}$", ylabel="Demixing free energy (eV/atom)", xlim=(0, 1), ylim=(-0.25, 0.05))
    a.legend(frameon=False, fontsize=6.5, loc="lower left", handletextpad=0.2)
    style(a, "(a)")

    # PWSCF-PAW: the series whose curvature the sheet differentiated
    a0, b0, c0, d0 = p["alpha"][0], p["beta"][0], p["gamma"][0], p["delta"][0]
    xb = d["x_cr_b"]
    g_nodes, _ = levesque_model(xb, a0, b0, c0, d0)
    resid = np.nanmax(abs(g_nodes - d["pwscf_paw"]))
    xs = np.linspace(0.005, 0.995, 800)
    _, curv = levesque_model(xs, a0, b0, c0, d0)
    sel = ~np.isnan(d["sheet_second_derivative_selected"])
    xsel = np.unique(xb[sel])
    pts = numerical_curvature(xb, d["pwscf_paw"], xsel)
    sheet = d["sheet_second_derivative_selected"]
    xo, po = xb[sel], sheet[sel]          # the sheet's own column P (two-point scheme)
    quad = np.polyval(np.polyfit(xo, po, 2), xs)
    b.fill_between(xs, 0, curv, where=curv > 0, color=GREEN, alpha=0.15, lw=0)
    b.axhline(0, color="black", lw=0.8)
    b.plot(xs, quad, "--", color="0.6", lw=1.2, label="quadratic fit to two-point values")
    b.plot(xo, po, "o", mfc="white", mec="0.5", mew=1.0, ms=4.5, ls="none", label="two-point differences")
    b.plot(xs, curv, color="0.2", lw=1.4, label="exact (PWSCF-PAW model)")
    b.plot(xsel, pts, "o", color="black", ms=4, ls="none", label="three-point differences")
    b.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"$\partial^2\Delta G/\partial x_{\rm Cr}^2$ (eV/atom)",
          xlim=(0, 1), ylim=(-10, 5))
    b.legend(frameon=False, fontsize=6.0, loc="lower right")
    b.text(0.03, 4.3, "unstable: curvature > 0 (shaded)", ha="left", va="top", fontsize=6.3)
    style(b, "(b)")
    save(fig, "fig1_free_energy")
    return resid, dict(zip(np.round(xsel, 3), np.round(pts, 3)))


def fig2(lit, dft, mdd, band, mdr, lam=5.5, stem="fig2_elastic_mapping", band_lam=5.5, ylim_b=None):
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.0))
    a.plot(lit["x"], lit["c11"], "-", color=BLUE, alpha=0.45, lw=0.8, zorder=1)
    a.plot(lit["x"], lit["c11"], "o", color=BLUE, ms=3.8, ls="none", label="Razumovskiy et al. [8]", zorder=3)
    if mdr is not None:
        a.plot(mdr["x"], mdr["c11"], "-", color=RED, alpha=0.45, lw=0.8, zorder=1)
        a.errorbar(mdr["x"], mdr["c11"], yerr=mdr["c11_sd"], fmt="s", color=RED, ms=4.5, capsize=2,
                   ls="none", label="MD, random alloy", zorder=3)
        # fully decomposed (phi = 1): Acta 50^3 at 50/50, present 25^3 cells elsewhere
        dpath = DATA / "md_decomposition.csv"
        xd, cd = [0.5], [mdd["c11_phi1"]]
        if dpath.exists():
            dd = read_csv(dpath.name)
            for x, ph, c in zip(dd["x_cr"], dd["phi"], dd["c11"]):
                if abs(ph - 1.0) < 1e-9 and abs(x - 0.5) > 1e-9:
                    xd.append(x); cd.append(c)
        a.plot(xd, cd, "s", mfc="white", mec=RED, mew=1.2, ms=5.5, ls="none",
               label=r"MD, fully decomposed ($\phi$ = 1)", zorder=4)
    else:
        a.axhline(mdd["c11_phi1"], color=RED, ls="--", lw=1.3,
                  label=r"MD spinodal, Fe$_{50}$Cr$_{50}$, $\phi$ = 1")
    a.plot(dft["x_cr"], dft["c11_gpa"], "D", color=YELLOW, mec="0.3", mew=0.5, ms=5.5, ls="none",
           label="DFT (this work)")
    a.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"$C_{11}$ (GPa)", xlim=(-0.03, 1.03), ylim=(200, 600))
    h, l = a.get_legend_handles_labels()
    order = [l.index(k) for k in ("Razumovskiy et al. [8]", "DFT (this work)", "MD, random alloy",
                                  r"MD, fully decomposed ($\phi$ = 1)") if k in l]
    order += [i for i in range(len(l)) if i not in order]
    a.legend([h[i] for i in order], [l[i] for i in order], frameon=False, fontsize=6.5, loc="upper left")
    style(a, "(a)")

    # Fig. 2b: the original draft's 5.5-um basis (sheet columns EB/EC = 6.4-um values x 6.4/5.5).
    # The band edges are Acta Fig. 9 values, defined at 6.4 um (pass band_lam=6.4 for a consistent plot).
    k = 6.4 / lam
    kb = band_lam / lam
    lo, hi = band["unallowed_min_ghz"][0] * kb, band["unallowed_max_ghz"][0] * kb
    b.axhspan(lo, hi, color="0.9", lw=0)
    b.axhline(hi, color=BLUE, lw=1.2, label="unallowed range maximum")
    b.axhline(lo, color=RED, lw=1.2, label="unallowed range minimum")
    b.plot(lit["c11"], lit["f"] * k, "o", color=YELLOW, mec="0.3", mew=0.4, ms=4, ls="none",
           label="Razumovskiy et al. [8]")
    b.plot(dft["c11_gpa"], dft["fsaw_6p4um_ghz"] * k, "D", color=GREEN, ms=5, ls="none", label="DFT (this work)")
    b.text(590, 0.822 * kb, "unallowed for a\nhomogeneous\nsolid solution", ha="right", va="top",
           fontsize=6.5, color="0.3")
    b.set(xlabel=r"$C_{11}$ (GPa)", ylabel=rf"$f_{{\rm SAW}}$ (GHz), $\lambda$ = {lam:g} $\mu$m",
          xlim=(200, 600), ylim=ylim_b or (0.50 * kb, 0.85 * kb))
    b.legend(frameon=False, fontsize=6.3, loc="upper left", bbox_to_anchor=(0.0, 0.9))
    style(b, "(b)")
    save(fig, stem)


def saw_ranges(prefix):
    """Orientation-resolved (P)SAW frequencies at 6.4 um from data/saw/<prefix>_x*.csv
    (Du & Zhao method, 10k random grains, PSAW kept if intensity >= 2e-3 of SAW, as in Acta)."""
    import glob
    out = []
    for path in sorted(glob.glob(str(DATA / "saw" / f"{prefix}_x*.csv"))):
        d = np.loadtxt(path, delimiter=",", skiprows=1)
        keep = d[:, 3] / d[:, 2] >= 2e-3
        modes = np.r_[d[:, 0], d[keep, 1]] / 6.4e-6 / 1e9
        x = float(Path(path).stem.split("_x")[1])
        out.append((x, modes.min(), np.percentile(modes, 5), np.percentile(modes, 95), modes.max()))
    return np.array(out)


def fig2_orientation(lit, dft, mdd, mdr):
    """Fig. 2 with panel (b) orientation-resolved, consistent with Acta Fig. 6 (6.4 um)."""
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.0))
    a.plot(lit["x"], lit["c11"], "-", color=BLUE, alpha=0.45, lw=0.8, zorder=1)
    a.plot(lit["x"], lit["c11"], "o", color=BLUE, ms=3.8, ls="none", label="Razumovskiy et al. [8]", zorder=3)
    if mdr is not None:
        a.plot(mdr["x"], mdr["c11"], "-", color=RED, alpha=0.45, lw=0.8, zorder=1)
        a.errorbar(mdr["x"], mdr["c11"], yerr=mdr["c11_sd"], fmt="s", color=RED, ms=4.5, capsize=2,
                   ls="none", label="MD, random alloy", zorder=3)
    a.plot(dft["x_cr"], dft["c11_gpa"], "D", color=YELLOW, mec="0.3", mew=0.5, ms=5.5, ls="none",
           label="DFT (this work)")
    a.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"$C_{11}$ (GPa)", xlim=(-0.03, 1.03), ylim=(200, 600))
    a.legend(frameon=False, fontsize=6.5, loc="upper left")
    style(a, "(a)")

    b.axhspan(0.609, 0.833, color="0.9", lw=0, zorder=0)
    b.axhline(0.609, color="0.35", lw=0.9, ls="--", zorder=1)
    series = [("lit", "black", "o", -0.012, "Razumovskiy et al. [8]"),
              ("dft", RED, "D", 0.0, "DFT (this work)"),
              ("md", GREEN, "s", 0.012, "MD, random alloy")]
    for prefix, col, mk, dx, lab in series:
        r = saw_ranges(prefix)
        if not len(r):
            continue
        xs = r[:, 0] + dx
        b.vlines(xs, r[:, 1], r[:, 4], color=col, lw=0.8, zorder=2)
        b.vlines(xs, r[:, 2], r[:, 3], color=col, lw=3.2, alpha=0.75, zorder=3)
        b.plot(xs, r[:, 4], mk, color=col, ms=3.5, ls="none", zorder=4, label=lab)
    b.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"(P)SAW frequency (GHz), $\lambda$ = 6.4 $\mu$m",
          xlim=(-0.05, 1.05), ylim=(0.33, 0.86))
    from matplotlib.patches import Patch
    h, l = b.get_legend_handles_labels()
    h.append(Patch(facecolor="0.9", edgecolor="0.35", ls="--", lw=0.8)); l.append("unallowed for pristine CASS [1]")
    b.legend(h, l, frameon=False, fontsize=6.3, loc="upper left", title="marker = fastest grain",
             title_fontsize=6.0)
    style(b, "(b)")
    save(fig, "fig2_elastic_mapping_orientation")


def saw_c11_lookup():
    """C11 (GPa) of every tensor run through the SAW calculation, keyed by file tag."""
    tab = {"dft_x0.00": 324.434353, "dft_x0.20": 400.025537, "dft_x0.50": 347.46247,
           "dft_x0.95": 499.051533, "dft_x1.00": 568.728963}
    for name in ("lit_tensors.csv", "lit_tensors_extra.csv", "md_tensors.csv"):
        path = DATA / "saw" / name
        if path.exists():
            for line in path.read_text().split("\n"):
                if line.strip():
                    tag, _x, c11, *_ = line.split(",")
                    tab[tag] = float(c11)
    return tab


def fig2_orientation_c11(lit, dft, mdd, mdr):
    """Fig. 2 variant: panel (b) as orientation-resolved (P)SAW frequency vs C11 (6.4 um)."""
    import glob
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.0))
    a.plot(lit["x"], lit["c11"], "-", color=BLUE, alpha=0.45, lw=0.8, zorder=1)
    a.plot(lit["x"], lit["c11"], "o", color=BLUE, ms=3.8, ls="none", label="Razumovskiy et al. [8]", zorder=3)
    if mdr is not None:
        a.plot(mdr["x"], mdr["c11"], "-", color=RED, alpha=0.45, lw=0.8, zorder=1)
        a.errorbar(mdr["x"], mdr["c11"], yerr=mdr["c11_sd"], fmt="s", color=RED, ms=4.5, capsize=2,
                   ls="none", label="MD, random alloy", zorder=3)
    a.plot(dft["x_cr"], dft["c11_gpa"], "D", color=YELLOW, mec="0.3", mew=0.5, ms=5.5, ls="none",
           label="DFT (this work)")
    a.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"$C_{11}$ (GPa)", xlim=(-0.03, 1.03), ylim=(200, 600))
    a.legend(frameon=False, fontsize=6.5, loc="upper left")
    style(a, "(a)")

    c11 = saw_c11_lookup()
    b.axhspan(0.609, 0.833, color="0.9", lw=0, zorder=0)
    b.axhline(0.609, color="0.35", lw=0.9, ls="--", zorder=1)
    for prefix, col, mk, lab in (("lit", "black", "o", "Razumovskiy et al. [8]"),
                                 ("dft", RED, "D", "DFT (this work)"),
                                 ("md", GREEN, "s", "MD, random alloy")):
        pts = []
        for path in sorted(glob.glob(str(DATA / "saw" / f"{prefix}_x*.csv"))):
            tag = Path(path).stem
            if tag not in c11:
                continue
            d = np.loadtxt(path, delimiter=",", skiprows=1)
            keep = d[:, 3] / d[:, 2] >= 2e-3
            modes = np.r_[d[:, 0], d[keep, 1]] / 6.4e-6 / 1e9
            pts.append((c11[tag], modes.min(), *np.percentile(modes, [5, 95]), modes.max()))
        if not pts:
            continue
        r = np.array(sorted(pts))
        b.vlines(r[:, 0], r[:, 1], r[:, 4], color=col, lw=0.8, zorder=2)
        b.vlines(r[:, 0], r[:, 2], r[:, 3], color=col, lw=3.2, alpha=0.75, zorder=3)
        b.plot(r[:, 0], r[:, 4], mk, color=col, ms=3.5, ls="none", zorder=4, label=lab)
    b.set(xlabel=r"$C_{11}$ (GPa)", ylabel=r"(P)SAW frequency (GHz), $\lambda$ = 6.4 $\mu$m",
          xlim=(200, 600), ylim=(0.33, 0.86))
    from matplotlib.patches import Patch
    h, l = b.get_legend_handles_labels()
    h.append(Patch(facecolor="0.9", edgecolor="0.35", ls="--", lw=0.8)); l.append("unallowed for pristine CASS [1]")
    b.legend(h, l, frameon=False, fontsize=6.3, loc="upper left", title="marker = fastest grain",
             title_fontsize=6.0)
    style(b, "(b)")
    save(fig, "fig2_elastic_mapping_orientation_c11")


def fig2_orientation_violin(lit, dft, mdd, mdr):
    """Fig. 2 variant: (b) violins of orientation-resolved (P)SAW frequency at each tensor's C11."""
    import glob
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1, 1.25]})
    a.plot(lit["x"], lit["c11"], "-", color=BLUE, alpha=0.45, lw=0.8, zorder=1)
    a.plot(lit["x"], lit["c11"], "o", color=BLUE, ms=3.8, ls="none", label="Razumovskiy et al. [8]", zorder=3)
    if mdr is not None:
        a.plot(mdr["x"], mdr["c11"], "-", color=RED, alpha=0.45, lw=0.8, zorder=1)
        a.errorbar(mdr["x"], mdr["c11"], yerr=mdr["c11_sd"], fmt="s", color=RED, ms=4.5, capsize=2,
                   ls="none", label="MD, random alloy", zorder=3)
    a.plot(dft["x_cr"], dft["c11_gpa"], "D", color=YELLOW, mec="0.3", mew=0.5, ms=5.5, ls="none",
           label="DFT (this work)")
    a.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"$C_{11}$ (GPa)", xlim=(-0.03, 1.03), ylim=(200, 600))
    a.legend(frameon=False, fontsize=6.5, loc="upper left")
    style(a, "(a)")

    c11 = saw_c11_lookup()
    b.axhspan(0.609, 0.86, color="0.92", lw=0, zorder=0)
    b.axhline(0.609, color="0.35", lw=0.9, ls="--", zorder=1)
    from matplotlib.patches import Patch
    handles = []
    for prefix, col, lab in (("lit", "0.15", "Razumovskiy et al. [8]"), ("dft", RED, "DFT (this work)"),
                             ("md", GREEN, "MD, random alloy")):
        for path in sorted(glob.glob(str(DATA / "saw" / f"{prefix}_x*.csv"))):
            tag = Path(path).stem
            if tag not in c11:
                continue
            d = np.loadtxt(path, delimiter=",", skiprows=1)
            keep = d[:, 3] / d[:, 2] >= 2e-3
            modes = np.r_[d[:, 0], d[keep, 1]] / 6.4e-6 / 1e9
            pos = c11[tag]
            v = b.violinplot(modes, positions=[pos], widths=8, showextrema=False, bw_method=0.15)
            for body in v["bodies"]:
                body.set_facecolor(col); body.set_edgecolor(col); body.set_alpha(0.35); body.set_linewidth(0.5)
            q1, med, q3 = np.percentile(modes, [25, 50, 75])
            b.vlines(pos, np.percentile(modes, 5), np.percentile(modes, 95), color=col, lw=0.8, zorder=3)
            b.add_patch(plt.Rectangle((pos - 2.2, q1), 4.4, q3 - q1, facecolor="white", edgecolor=col,
                                      lw=0.8, zorder=4))
            b.hlines(med, pos - 2.2, pos + 2.2, color=col, lw=1.2, zorder=5)
            frac = 100 * np.mean(modes > 0.609)
            if frac >= 1:
                b.text(pos, modes.max() + 0.008, f"{frac:.0f}%", ha="center", va="bottom", fontsize=5.8, color=col)
        handles.append(Patch(facecolor=col, alpha=0.35, edgecolor=col, label=lab))
    handles.append(Patch(facecolor="0.92", edgecolor="0.35", ls="--", label="unallowed for pristine CASS [1]"))
    b.set(xlabel=r"$C_{11}$ (GPa)", ylabel=r"(P)SAW frequency (GHz), $\lambda$ = 6.4 $\mu$m",
          xlim=(225, 595), ylim=(0.34, 0.86))
    b.legend(handles=handles, frameon=False, fontsize=6.0, loc="upper left")
    b.text(0.98, 0.03, "box: quartiles; whiskers: 5-95%\n%: share of modes in band",
           transform=b.transAxes, ha="right", va="bottom", fontsize=5.6, color="0.35")
    style(b, "(b)")
    save(fig, "fig2_elastic_mapping_orientation_violin")


def saw_tensors():
    """tag -> (x, C11, C12, C44) for every tensor in the orientation-resolved SAW runs."""
    tab = {"dft_x0.00": (0.0, 324.434353, 169.838953, 127.7242), "dft_x0.20": (0.2, 400.025537, 164.57738, 120.574293),
           "dft_x0.50": (0.5, 347.46247, 186.14903, 123.354977), "dft_x0.95": (0.95, 499.051533, 145.051597, 99.27015),
           "dft_x1.00": (1.0, 568.728963, 158.14884, 100.913353)}
    for name in ("lit_tensors.csv", "lit_tensors_extra.csv", "md_tensors.csv", "mddec_tensors.csv"):
        path = DATA / "saw" / name
        if path.exists():
            for line in path.read_text().split("\n"):
                if line.strip():
                    tag, x, c11, c12, c44 = line.split(",")
                    tab[tag] = (float(x), float(c11), float(c12), float(c44))
    return tab


def bulk_longitudinal_ghz(x, c11, c12, c44, n=4000):
    """Quasi-longitudinal bulk-wave frequency (6.4 um) along the in-plane propagation direction of the
    first n grains of the 10k set (same direction as saw_du.get_saw with deg = 0)."""
    import saw_du as sd
    cache = DATA / "saw" / f"bulkL_{x:.3f}_{c11:.2f}_{c12:.2f}_{c44:.2f}.npy"
    if cache.exists():
        return np.load(cache)
    ang = np.loadtxt(DATA / "saw" / "angles_10k.csv", delimiter=",")[:n]
    C = sd.cubic_cijkl(c11 * 1e9, c12 * 1e9, c44 * 1e9)
    rho = 7870 * (1 - x) + 7200 * x
    out = []
    for e in ang:
        nv = sd.euler2matrix(*e) @ np.array([0, 1.0, 0])
        out.append(np.sqrt(np.linalg.eigvalsh(np.einsum("ijkl,j,l->ik", C, nv, nv) / rho)[-1]))
    f = np.array(out) / 6.4e-6 / 1e9
    np.save(cache, f)
    return f


def fig2_orientation_bulk(lit, dft, mdd, mdr, kind="violin", with_decomp=False):
    """Orientation-resolved (P)SAW (solid) and quasi-longitudinal bulk waves (hatched) vs C11."""
    import glob
    sys_path = str(ROOT / "scripts")
    import sys
    if sys_path not in sys.path:
        sys.path.insert(0, sys_path)
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.2), gridspec_kw={"width_ratios": [1, 1.25]})
    a.plot(lit["x"], lit["c11"], "-", color=BLUE, alpha=0.45, lw=0.8, zorder=1)
    a.plot(lit["x"], lit["c11"], "o", color=BLUE, ms=3.8, ls="none", label="Razumovskiy et al. [8]", zorder=3)
    if mdr is not None:
        a.plot(mdr["x"], mdr["c11"], "-", color=RED, alpha=0.45, lw=0.8, zorder=1)
        a.errorbar(mdr["x"], mdr["c11"], yerr=mdr["c11_sd"], fmt="s", color=RED, ms=4.5, capsize=2,
                   ls="none", label="MD, random alloy", zorder=3)
    a.plot(dft["x_cr"], dft["c11_gpa"], "D", color=YELLOW, mec="0.3", mew=0.5, ms=5.5, ls="none",
           label="DFT (this work)")
    a.set(xlabel=r"$x_{\rm Cr}$", ylabel=r"$C_{11}$ (GPa)", xlim=(-0.03, 1.03), ylim=(200, 600))
    a.legend(frameon=False, fontsize=6.5, loc="upper left")
    style(a, "(a)")

    tens = saw_tensors()
    b.axhspan(0.609, 0.833, color="0.92", lw=0, zorder=0)
    b.axhline(0.609, color="0.35", lw=0.9, ls="--", zorder=1)
    b.axhline(0.833, color="0.35", lw=0.9, ls=":", zorder=1)
    from matplotlib.patches import Patch
    handles = []
    series = [("lit", "0.15", "Razumovskiy et al. [8]"), ("dft", RED, "DFT (this work)"),
              ("md", GREEN, "MD, random alloy")]
    if with_decomp:
        series.append(("mddec", "#0B5D3B", r"MD, fully decomposed ($\phi$ = 1)"))
    for prefix, col, lab in series:
        for path in sorted(glob.glob(str(DATA / "saw" / f"{prefix}_x*.csv"))):
            tag = Path(path).stem
            if tag not in tens:
                continue
            x, c11, c12, c44 = tens[tag]
            d = np.loadtxt(path, delimiter=",", skiprows=1)
            keep = d[:, 3] / d[:, 2] >= 2e-3
            modes = np.r_[d[:, 0], d[keep, 1]] / 6.4e-6 / 1e9
            bulk = bulk_longitudinal_ghz(x, c11, c12, c44)
            for data, hatch, alpha in ((modes, None, 0.35), (bulk, "////", 0.15)):
                if kind == "violin":
                    v = b.violinplot(data, positions=[c11], widths=8, showextrema=False, bw_method=0.15)
                    for body in v["bodies"]:
                        body.set_facecolor(col); body.set_edgecolor(col); body.set_alpha(alpha)
                        body.set_linewidth(0.5); body.set_hatch(hatch)
                    q1, med, q3 = np.percentile(data, [25, 50, 75])
                    b.add_patch(plt.Rectangle((c11 - 2.2, q1), 4.4, q3 - q1, facecolor="white", edgecolor=col,
                                              lw=0.7, zorder=4))
                    b.hlines(med, c11 - 2.2, c11 + 2.2, color=col, lw=1.1, zorder=5)
                else:
                    b.vlines(c11, data.min(), data.max(), color=col, lw=0.8, zorder=2,
                             linestyles="-" if hatch is None else ":")
                    b.vlines(c11, *np.percentile(data, [5, 95]), color=col, lw=3.2,
                             alpha=0.75 if hatch is None else 0.35, zorder=3)
        handles.append(Patch(facecolor=col, alpha=0.35, edgecolor=col, label=lab))

    b.set(xlabel=r"$C_{11}$ (GPa)", ylabel=r"Frequency (GHz), $\lambda$ = 6.4 $\mu$m",
          xlim=(225, 595), ylim=(0.34, 1.32))
    b.legend(handles=handles, frameon=False, fontsize=6.0, loc="upper left", handlelength=1.4)
    b.text(590, 0.95, "bulk (quasi-longitudinal)\nwaves", ha="right", va="center", fontsize=6.3, color="0.35",
           style="italic")
    b.text(272, 0.72, "unallowed for\npristine CASS [1]", ha="center", va="center", fontsize=6.0, color="0.35")
    b.text(590, 0.36, "(P)SAW", ha="right", va="bottom", fontsize=6.3, color="0.35", style="italic")
    style(b, "(b)")
    save(fig, f"fig2_elastic_mapping_orientation_{kind}_bulk" + ("_mddecomp" if with_decomp else ""))


def fig3(lit, dft, mdd, mdr, lam=6.4, stem="fig3_fsaw", rel_b=False, crcolor=False, md_decomp_a=False,
         alpha_marks=False):
    k = 6.4 / lam  # f scales as 1/lambda; layout scales with it
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.0))
    a.plot(lit["x"], lit["f"] * k, "-", color="0.6", lw=0.8, zorder=1)
    a.plot(lit["x"], lit["f"] * k, "o", color="black", ms=3.6, ls="none", label="Razumovskiy et al. [8]", zorder=3)
    a.plot(dft["x_cr"], dft["fsaw_6p4um_ghz"] * k, "o", color=RED, mec="white", mew=0.6, ms=6, ls="none",
           label="DFT (this work)", zorder=4)
    if mdr is not None:
        a.plot(mdr["x"], mdr["f"] * k, "-", color=GREEN, alpha=0.5, lw=0.8, zorder=1)
        a.plot(mdr["x"], mdr["f"] * k, "s", color=GREEN, mec="white", mew=0.6, ms=5.5, ls="none",
               label="MD, random alloy", zorder=4)
    if md_decomp_a and mdr is not None:
        dec_a = md_decomposition()
        first = True
        for x, ser in dec_a.items():
            o = np.argsort(ser["phi"])
            ph, fs = ser["phi"][o], ser["f"][o] * k
            f_rand = float(np.interp(x, mdr["x"], mdr["f"])) * k     # 50^3 random alloy at this x
            shift = f_rand - fs[0]                                    # anchor the ladder to that point
            ys = fs + shift
            a.plot([x, x], [ys[0], ys[-1]], ls=":", color=GREEN, lw=1.0, zorder=3)
            a.scatter(np.full(len(ph) - 2, x), ys[1:-1], s=7, c=GREEN, alpha=0.75, zorder=4, lw=0,
                      label=r"MD, $\phi$ = 0.6$-$0.9" if first else None)
            a.plot([x], [ys[-1]], "s", mfc="white", mec=GREEN, mew=1.2, ms=5.5, ls="none", zorder=5,
                   label=r"MD, fully decomposed ($\phi$ = 1)" if first else None)
            first = False
    # Partitioning schematic in its own strip below the data (no overlap with any series)
    xa, xb = 0.15, 0.85
    if alpha_marks:   # MD-measured alpha / alpha' compositions of decomposed Fe-36Cr (phi = 1)
        pc = read_csv("md_phase_compositions.csv")
        m = (np.abs(pc["x_cr"] - X_FE36CR) < 1e-6) & (pc["phi"] == 1.0)
        xa, xb = float(pc["x_alpha"][m][0]), float(pc["x_alpha_prime"][m][0])
        for xv, lab in ((xa, r"$\alpha$"), (xb, r"$\alpha'$")):
            fv = k * float(np.interp(xv, lit["x"], lit["f"]))
            a.plot([xv], [fv], "v", color=BLUE, mec="white", mew=0.6, ms=7, zorder=6)
    y0 = 0.428 * k
    for xv in (xa, X_FE36CR, xb):
        fv = k * float(np.interp(xv, lit["x"], lit["f"]))
        if mdr is not None:  # stop below the lowest series at this composition
            fv = min(fv, k * float(np.interp(xv, mdr["x"], mdr["f"])))
        a.plot([xv, xv], [y0 + 0.013 * k, fv - 0.009 * k], ls=(0, (1, 2)), color=BLUE, lw=0.9, alpha=0.8, zorder=0)
    arrow = dict(arrowstyle="-|>", lw=1.4, color=BLUE, shrinkA=7, shrinkB=0, mutation_scale=10)
    a.annotate("", xy=(xa, y0), xytext=(X_FE36CR, y0), arrowprops=arrow)
    a.annotate("", xy=(xb, y0), xytext=(X_FE36CR, y0), arrowprops=arrow)
    a.plot([X_FE36CR], [y0], "*", color=BLUE, ms=11, zorder=5, clip_on=False)
    a.text(X_FE36CR, y0 - 0.008 * k, "Fe–36Cr", ha="center", va="top", fontsize=6.6, color=BLUE)
    if alpha_marks:
        a.text(xa, y0 - 0.008 * k, rf"$\alpha$: {round(100 * xa)}% Cr", ha="center", va="top", fontsize=6.6, color=BLUE)
        a.text(xb, y0 - 0.008 * k, rf"$\alpha'$: {round(100 * xb)}% Cr", ha="center", va="top", fontsize=6.6,
               color=BLUE)
    else:
        a.text(xa, y0 - 0.008 * k, r"Fe-rich $\alpha$", ha="center", va="top", fontsize=6.6, color=BLUE)
        a.text(xb, y0 - 0.008 * k, r"Cr-rich $\alpha'$", ha="center", va="top", fontsize=6.6, color=BLUE)
    a.set(xlabel=r"$x_{\rm Cr}$", ylabel=rf"$f_{{\rm SAW}}$ (GHz), $\lambda$ = {lam:g} $\mu$m",
          xlim=(-0.03, 1.03), ylim=(0.405 * k, 0.64 * k))
    a.legend(frameon=False, fontsize=6.5, loc="upper left")
    style(a, "(a)")

    dec = md_decomposition()

    def yb(phi, f):  # absolute frequency, or % change relative to the random alloy (phi = 0.5)
        return 100 * (f / f[np.argmin(phi)] - 1) if rel_b else f * k
    if crcolor:
        cmap = plt.get_cmap("viridis")
        cx = lambda x: cmap(0.0 + 0.78 * x)          # colour by Cr fraction (dark = Fe-rich, yellow = Cr-rich)
        ends = []                                     # (y at phi=1, x) for direct labels
        ya = yb(mdd["phi"], mdd["f"])
        b.plot(mdd["phi"], ya, "-", color=cx(0.5), lw=1.1, zorder=2)
        b.plot(mdd["phi"], ya, "o", color=cx(0.5), mec="black", mew=0.4, ms=5, ls="none", zorder=4,
               label=r"50$^3$ cell (Acta)")
        ends.append((ya[np.argmax(mdd["phi"])], 0.5))
        for x, ser in dec.items():
            ys = yb(ser["phi"], ser["f"])
            if abs(x - 0.5) < 1e-9:   # small-cell validation of the Acta series
                b.plot(ser["phi"], ys, "o", mfc="none", mec=cx(0.5), mew=1.0, ms=6.5, ls="none", zorder=3,
                       label=r"25$^3$ cell (validation)")
                continue
            b.plot(ser["phi"], ys, "-", color=cx(x), lw=1.1, zorder=2)
            b.plot(ser["phi"], ys, "s", color=cx(x), mec="black", mew=0.4, ms=4.5, ls="none", zorder=4)
            ends.append((ys[np.argmax(ser["phi"])], x))
        # direct labels at phi = 1, nudged apart vertically
        ends.sort()
        yr = np.ptp([e[0] for e in ends]) or 1.0
        placed = []
        for yv, x in ends:
            yl = yv if not placed or yv - placed[-1] > 0.06 * yr else placed[-1] + 0.06 * yr
            placed.append(yl)
            b.text(1.012, yl, f"{round(100 * x)}% Cr", color=cx(x), fontsize=6.3, va="center", ha="left",
                   fontweight="bold")
    else:
        b.plot(mdd["phi"], yb(mdd["phi"], mdd["f"]), "-", color="0.6", lw=0.8, zorder=1)
        b.plot(mdd["phi"], yb(mdd["phi"], mdd["f"]), "o", color="black", ms=5, ls="none", zorder=3,
               label=r"Fe$_{50}$Cr$_{50}$ (Acta, 50$^3$ cell)")
        colors = {0.2: BLUE, 0.36: PURPLE, 0.95: ORANGE}   # compositions named in the caption
        for x, ser in dec.items():
            if round(x, 2) not in colors and abs(x - 0.5) > 1e-9:
                continue   # the other compositions are shown in the _crcolor variants
            if abs(x - 0.5) < 1e-9:   # small-cell validation of the Acta series
                b.plot(ser["phi"], yb(ser["phi"], ser["f"]), "o", mfc="none", mec="0.55", ms=6, ls="none", zorder=2,
                       label=r"Fe$_{50}$Cr$_{50}$, 25$^3$ cell (validation)")
                continue
            col = colors.get(round(x, 2), RED)
            b.plot(ser["phi"], yb(ser["phi"], ser["f"]), "-", color=col, alpha=0.5, lw=0.8, zorder=1)
            b.plot(ser["phi"], yb(ser["phi"], ser["f"]), "s", color=col, ms=4.5, ls="none", zorder=3,
                   label=rf"Fe$_{{{100 - round(100 * x)}}}$Cr$_{{{round(100 * x)}}}$, 25$^3$ cell")
    b.set(xlabel=r"Spinodal conversion fraction, $\phi$", ylabel=rf"$f_{{\rm SAW}}$ (GHz), $\lambda$ = {lam:g} $\mu$m",
          xlim=(0.48, 1.13 if crcolor else 1.02))
    b.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    if rel_b:
        b.set_ylabel(r"$\Delta f_{\rm SAW}/f_{\rm SAW}(\phi=0.5)$ (%)")
    if dec:
        allf = np.concatenate([yb(mdd["phi"], mdd["f"])] + [yb(v["phi"], v["f"]) for v in dec.values()])
        rng = allf.max() - allf.min()
        b.set_ylim(allf.min() - 0.08 * rng, allf.max() + ((0.30 if rel_b else 0.45) if crcolor else 0.55) * rng)
        if crcolor:
            b.plot([], [], "s", color="0.5", mec="black", mew=0.4, ms=4.5, ls="none", label=r"25$^3$ cell")
            b.legend(frameon=False, fontsize=6.2, loc="upper left")
        else:
            import re
            h, l = b.get_legend_handles_labels()
            key = lambda t: (int(re.search(r"Cr\$_\{(\d+)\}\$", t).group(1)) if re.search(r"Cr\$_\{(\d+)\}\$", t) else 0,
                             "validation" in t)
            o = sorted(range(len(l)), key=lambda i: key(l[i]))
            b.legend([h[i] for i in o], [l[i] for i in o], frameon=False, fontsize=6.0, loc="upper left")
    else:
        b.set_ylim(0.46 * k, 0.50 * k)
        b.text(0.03, 0.04, r"MD, Fe$_{50}$Cr$_{50}$", transform=b.transAxes, fontsize=6.8)
    style(b, "(b)")
    save(fig, stem)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--xlsx", help="Google Sheet .xlsx export; refreshes data/sheet_*.csv")
    args = ap.parse_args()
    if args.xlsx:
        snapshot(args.xlsx)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "pdf.fonttype": 42,
                         "axes.linewidth": 0.7})
    lit, dft, mdd = literature(), read_csv("sheet_dft.csv"), md()
    band = read_csv("sheet_band.csv")
    resid, pts = fig1()
    mdr = md_random()
    fig2_orientation(lit, dft, mdd, mdr)      # candidate: (b) vs x_Cr
    fig2_orientation_c11(lit, dft, mdd, mdr)  # candidate: (b) vs C11
    fig2_orientation_violin(lit, dft, mdd, mdr)  # candidate: violins vs C11
    fig2_orientation_bulk(lit, dft, mdd, mdr, kind="violin")  # + bulk waves
    fig2_orientation_bulk(lit, dft, mdd, mdr, kind="bars")    # + bulk waves
    if (DATA / "saw" / "mddec_tensors.csv").exists():
        fig2_orientation_bulk(lit, dft, mdd, mdr, kind="violin", with_decomp=True)
        fig2_orientation_bulk(lit, dft, mdd, mdr, kind="bars", with_decomp=True)
    fig2(lit, dft, mdd, band, mdr)
    fig2(lit, dft, mdd, band, mdr, lam=6.4, stem="fig2_elastic_mapping_6p4um_bandunchanged",
         band_lam=6.4, ylim_b=(0.44, 0.85))  # check version: band not converted
    fig2(lit, dft, mdd, band, mdr, lam=5.5, stem="fig2_elastic_mapping_5p5um_bandconverted",
         band_lam=6.4, ylim_b=(0.50, 1.0))  # band edges (Acta Fig. 9, 6.4 um) converted to 5.5 um
    fig3(lit, dft, mdd, mdr, lam=5.5)                              # manuscript (5.5 um, as Fig. 2b)
    fig3(lit, dft, mdd, mdr, lam=6.4, stem="fig3_fsaw_6p4um")      # alternative at the TGS 6.4 um
    fig3(lit, dft, mdd, mdr, lam=5.5, stem="fig3_fsaw_relative", rel_b=True)  # option: 3b as % change
    fig3(lit, dft, mdd, mdr, lam=5.5, stem="fig3_fsaw_crcolor", crcolor=True)            # (b) coloured by x_Cr
    fig3(lit, dft, mdd, mdr, lam=5.5, stem="fig3_fsaw_relative_crcolor", rel_b=True, crcolor=True)
    fig3(lit, dft, mdd, mdr, lam=5.5, stem="fig3_fsaw_mddecomp", md_decomp_a=True)       # (a) + decomposed MD
    fig3(lit, dft, mdd, mdr, lam=5.5, stem="fig3_fsaw_alpha", alpha_marks=True)          # (a) MD alpha/alpha' ends
    if mdr is not None:
        print('MD random x, C11, f:', np.round(mdr['x'], 2), np.round(mdr['c11'], 1), np.round(mdr['f'], 4))
    print(f"max |sheet PWSCF-PAW - analytic model| = {resid:.2e} eV/atom")
    print("Fig 1b numerical curvature:", pts)
    print("lit f (x, GHz):", [(round(a, 3), round(b, 4)) for a, b in zip(lit["x"], lit["f"])])
    print("lit C11 x=1:", round(lit["c11"][-1], 1))
    print("DFT f:", dft["fsaw_6p4um_ghz"])
    print("MD f:", np.round(mdd["f"], 4), "phi", mdd["phi"])
    print("Fe-36Cr lit f:", round(float(np.interp(X_FE36CR, lit["x"], lit["f"])), 4))


if __name__ == "__main__":
    main()
