#!/usr/bin/env python3
"""Backup validation plots: every variant of our calculations replotted on top of the literature it should match,
with the deviation printed on the plot. These are for checking the supplementary figures and choosing which variant
to keep; they are not themselves manuscript figures.

  validation_overlays.py phonons      V_phonon_disp_Fe, V_phonon_disp_Cr, V_phonon_dos_Fe, V_phonon_dos_Cr
  validation_overlays.py electrons    V_moments, V_bulk_modulus, V_nef, V_dos_Fe, V_bands_Fe (needs the digitized
                                      Korzhavyi 2009 / Klaver 2006 tables in refdata/ and pure_dos/ from Engaging)
  validation_overlays.py all

Output: manuscript/figures/kdense_scripta/supplementary/validation/*.png and validation_summary.csv (one row per
comparison: what, variant, reference, metric, value).
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ENG, REF = HERE / "engaging", HERE / "refdata"
OUT = HERE.parent / "manuscript" / "figures" / "kdense_scripta" / "supplementary" / "validation"
CM1_THZ = 0.0299792458
RY_EV = 13.605693
S2, S3 = np.sqrt(2), np.sqrt(3)
XB = dict(G=0.0, H=1.0, N=1.0 + S2 / 2, G2=1.0 + S2, P=1.0 + S2 + S3 / 2, H2=1.0 + S2 + S3)   # bcc G-H-N-G-P-H
SUMMARY = []


def rows(name):
    with open(REF / name) as fh:
        return list(csv.DictReader(l for l in fh if not l.startswith("#")))


def note(what, variant, ref, metric, value):
    SUMMARY.append(dict(what=what, variant=variant, reference=ref, metric=metric, value=f"{value:.4g}"))


def style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(HERE.parent / "scripts" / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta")
    return plt


# ----------------------------------------------------------------------------------------------------- phonons
def path_x(direction, zeta):
    """Position on the bcc G-H-N-G-P-H path (units 2 pi/a) of a point given in the experimental convention."""
    d = direction.replace(" ", "")
    if d in ("[z00]", "[00z]"):
        return XB["G"] + zeta
    if d == "[zz1]":
        return XB["H"] + zeta * S2
    if d in ("[zz0]", "[0zz]"):
        return XB["G2"] - zeta * S2
    if d == "[zzz]":
        return XB["G2"] + zeta * S3
    return None


def exp_points(name):
    pts = [(path_x(r["direction"], float(r["zeta"])), float(r["nu_THz"])) for r in rows(name)]
    return np.array([p for p in pts if p[0] is not None])


def bcc_branches(case):
    """(x on the bcc path, branches in THz) for a 1-atom bcc DFPT run."""
    gp = np.loadtxt(ENG / "phonon" / case / "freq.gp")
    return [(gp[:, 0], gp[:, 1:] * CM1_THZ)]


def afm_branches(case):
    """Commensurate AFM Cr (2-atom simple cubic, path G-X-M-G-R-X) mapped onto the bcc path where the sc segments
    coincide with bcc lines: sc G-X = bcc [z00] for z <= 1/2 and, folded, z >= 1/2 (drawn mirrored); sc M-G = bcc N-G;
    sc G-R = bcc G-P. The AFM cell also carries branches folded in from q + (1,0,0); these are all drawn."""
    gp = np.loadtxt(ENG / "phonon" / case / "freq.gp")
    x, w = gp[:, 0], gp[:, 1:] * CM1_THZ
    out = []
    m = x <= 0.5 + 1e-9
    out += [(x[m], w[m]), (1.0 - x[m], w[m])]
    m = (x >= 1.0 - 1e-9) & (x <= 1.0 + S2 / 2 + 1e-9)
    out.append((XB["N"] + (x[m] - 1.0), w[m]))
    m = (x >= 1.0 + S2 / 2 - 1e-9) & (x <= 1.0 + S2 / 2 + S3 / 2 + 1e-9)
    out.append((XB["G2"] + (x[m] - 1.0 - S2 / 2), w[m]))
    return out


def deviation(segs, pts):
    """Nearest-branch deviation at every experimental point covered by the calculated segments: RMS and mean (THz)."""
    d = []
    for xe, ne in pts:
        best = None
        for xs, ws in segs:
            if xs.min() - 1e-6 <= xe <= xs.max() + 1e-6:
                o = np.argsort(xs)
                cand = [np.interp(xe, xs[o], ws[o, b]) - ne for b in range(ws.shape[1])]
                c = min(cand, key=abs)
                best = c if best is None or abs(c) < abs(best) else best
        if best is not None:
            d.append(best)
    d = np.array(d)
    return (float(np.sqrt(np.mean(d ** 2))), float(np.mean(d)), len(d)) if len(d) else (np.nan, np.nan, 0)


def plot_dispersion(plt, out, title, refname, reflabel, variants):
    pts = exp_points(refname)
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    for case, getter, col, ls, lab in variants:
        if not (ENG / "phonon" / case / "freq.gp").exists():
            continue
        segs = getter(case)
        rms, mean, n = deviation(segs, pts)
        note(title, case, refname, "RMS nearest-branch deviation (THz)", rms)
        note(title, case, refname, "mean deviation calc-exp (THz)", mean)
        for i, (xs, ws) in enumerate(segs):
            for b in range(ws.shape[1]):
                ax.plot(xs, ws[:, b], ls, color=col, lw=0.9, alpha=0.9,
                        label=f"{lab}: RMS {rms:.2f}, mean {mean:+.2f} THz (n={n})" if (i == 0 and b == 0) else None)
    ax.plot(pts[:, 0], pts[:, 1], "o", mfc="white", mec="0.05", ms=3, mew=0.7, label=reflabel, zorder=5)
    for v in XB.values():
        ax.axvline(v, color="0.85", lw=0.5)
    ax.set_xticks(list(XB.values())); ax.set_xticklabels(["Γ", "H", "N", "Γ", "P", "H"]); ax.minorticks_off()
    ax.set(xlim=(0, XB["H2"]), ylim=(0, 12), ylabel="frequency (THz)", title=title)
    ax.legend(fontsize=5.6, loc="upper center", ncol=1, framealpha=0.92)
    fig.tight_layout(); fig.savefig(OUT / f"{out}.png", dpi=300); plt.close(fig)
    print("wrote", out)


def dfpt_dos(case, sigma=0.15):
    a = np.loadtxt(ENG / "phonon" / case / "dos.dat")
    f, g = a[:, 0] * CM1_THZ, a[:, 1]
    grid = np.linspace(0, 12, 1201)
    gi = np.interp(grid, f, g, left=0, right=0)
    k = np.exp(-0.5 * (np.arange(-60, 61) * 0.01 / sigma) ** 2)
    gs = np.convolve(gi, k / k.sum(), mode="same")
    return grid, gs / np.trapezoid(gs, grid)


def moments_of(f, g):
    """Mean frequency and the frequency of the highest peak (THz) of a normalised DOS."""
    return float(np.trapezoid(f * g, f) / np.trapezoid(g, f)), float(f[np.argmax(g)])


def plot_phdos(plt, out, title, refs, variants):
    fig, ax = plt.subplots(figsize=(4.6, 3.2))
    for lab, (f, g), col, ls in refs:
        mf, pk = moments_of(f, g)
        ax.plot(f, g, ls, color=col, lw=1.0, label=f"{lab}: <ν> {mf:.2f}, peak {pk:.2f} THz")
    rmean = moments_of(*refs[0][1])[0]
    for case, col, ls, lab in variants:
        if not (ENG / "phonon" / case / "dos.dat").exists():
            continue
        f, g = dfpt_dos(case)
        mf, pk = moments_of(f, g)
        note(title, case, refs[0][0], "mean frequency (THz)", mf)
        note(title, case, refs[0][0], "mean frequency vs first reference (%)", 100 * (mf / rmean - 1))
        ax.plot(f, g, ls, color=col, lw=1.1, label=f"{lab}: <ν> {mf:.2f} ({100 * (mf / rmean - 1):+.1f}%), peak {pk:.2f}")
    ax.set(xlim=(0, 11.5), ylim=(0, None), xlabel="ν (THz)", ylabel="g(ν) (THz$^{-1}$), unit area", title=title)
    ax.legend(fontsize=5.0, loc="upper left", framealpha=0.92)
    fig.tight_layout(); fig.savefig(OUT / f"{out}.png", dpi=300); plt.close(fig)
    print("wrote", out)


def phonons():
    plt = style()
    plot_dispersion(plt, "V_phonon_disp_Fe", "bcc Fe phonons vs neutron scattering (295 K)", "fe_phonon_minkiewicz1967.csv",
                    "neutron, 295 K (Minkiewicz 1967, tabulated)",
                    [("Fe_fm", bcc_branches, "#2A9DF4", "-", "QE DFPT FM, relaxed a"),
                     ("Fe_fm_aexp", bcc_branches, "#D55E00", "--", "QE DFPT FM, experimental a")])
    plot_dispersion(plt, "V_phonon_disp_Cr", "bcc Cr phonons vs neutron scattering (300 K, AFM)", "cr_phonon_shaw1971.csv",
                    "neutron, 300 K (Shaw & Muhlestein 1971, digitized)",
                    [("Cr_nm_v5_k24", bcc_branches, "#2A9DF4", "-", "QE DFPT NM ($k=24^3$), relaxed a"),
                     ("Cr_nm_aexp", bcc_branches, "#D55E00", "--", "QE DFPT NM, experimental a"),
                     ("Cr_afm", afm_branches, "#009E73", "-", "QE DFPT AFM (2-atom, folded), relaxed a"),
                     ("Cr_afm_aexp", afm_branches, "#CC79A7", "--", "QE DFPT AFM (2-atom, folded), experimental a")])
    mau = rows("fe_phonon_dos_mauger2014.csv")
    fu = rows("fecr_phonon_dos_fultz1995.csv")
    fd = np.array([[float(r["nu_THz"]), float(r["Fe"]), float(r["Cr"])] for r in fu])
    norm = lambda f, g: (f, g / np.trapezoid(g, f))
    refs_fe = []
    for T, ls in (("300", "-"), ("30", "--")):
        d = np.array([[float(r["nu_THz"]), float(r["g_perTHz"])] for r in mau if r["T_K"] == T])
        refs_fe.append((f"NRIXS {T} K (Mauger 2014)", norm(d[:, 0], d[:, 1]), "0.05", ls))
    refs_fe.append(("BvK from neutron (Fultz 1995)", norm(fd[:, 0], fd[:, 1]), "0.5", ":"))
    plot_phdos(plt, "V_phonon_dos_Fe", "Fe phonon DOS", refs_fe,
               [("Fe_fm", "#2A9DF4", "-", "QE FM, relaxed a"), ("Fe_fm_aexp", "#D55E00", "--", "QE FM, experimental a")])
    plot_phdos(plt, "V_phonon_dos_Cr", "Cr phonon DOS", [("BvK from neutron (Fultz 1995)", norm(fd[:, 0], fd[:, 2]), "0.05", "-")],
               [("Cr_nm_v5_k24", "#2A9DF4", "-", "QE NM, relaxed a"), ("Cr_nm_aexp", "#D55E00", "--", "QE NM, experimental a"),
                ("Cr_afm", "#009E73", "-", "QE AFM, relaxed a"), ("Cr_afm_aexp", "#CC79A7", "--", "QE AFM, experimental a")])


# --------------------------------------------------------------------------------------------------- electrons
def sqs_runs(tag):
    """{x: run dir} of the 40-atom SQS elastic runs: tag '' = relaxed PBE volume, 'aexp' = experimental lattice."""
    out = {}
    for d in sorted((ENG / "runs").glob("sqs40_x*_k442_e60*")):
        suf = d.name.split("_e60")[1].lstrip("_")
        if suf == tag and (d / "result.json").exists():
            out[json.loads((d / "result.json").read_text())["x_cr"]] = d
    return out


def site_moments(d):
    """Mean and std of the Fe and Cr site moments (muB) from the last 'Magnetic moment per site' block of relax.out."""
    txt = (d / "relax.out").read_text().splitlines()
    i = max(n for n, l in enumerate(txt) if "Magnetic moment per site" in l)
    mag = []
    for l in txt[i + 1:]:
        if "magn=" not in l:
            break
        mag.append(float(l.split("magn=")[1]))
    inp = (d / "relax.in").read_text().split("ATOMIC_POSITIONS")[1].splitlines()[1:]
    sym = [l.split()[0] for l in inp if l.strip() and l.split()[0] in ("Fe", "Cr")][:len(mag)]
    mag = np.array(mag)
    return {el: (float(mag[[s == el for s in sym]].mean()), float(mag[[s == el for s in sym]].std()))
            for el in ("Fe", "Cr") if el in sym}


def vrh_k(c11, c12):
    return (c11 + 2 * c12) / 3


def plot_moments(plt):
    kz = rows("korzhavyi2009_fig8_moments.csv"); kl = rows("klaver2006_fig4_moments.csv")
    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    cols = {"Fe": "#D55E00", "Cr": "#0072B2"}
    for el in ("Fe", "Cr"):
        d = np.array(sorted((float(r["c_Cr_at_pct"]) / 100, float(r["moment_muB"])) for r in kz if r["series"] == f"M_{el} FM"))
        ax.plot(d[:, 0], d[:, 1], "-" if el == "Cr" else "--", color="0.15", lw=1.0,
                label=f"EMTO-CPA FM, M_{el} (Korzhavyi 2009)")
        k = np.array([(float(r["c_Cr"]), float(r["moment_muB"])) for r in kl if r["series"] == el])
        ax.plot(k[:, 0], k[:, 1], "x", color="0.55", ms=3, mew=0.7, label=f"VASP supercells, {el} (Klaver 2006)")
        for tag, mfc, lab in (("", cols[el], "relaxed a"), ("aexp", "white", "experimental a")):
            runs = sqs_runs(tag)
            pts = [(x, *site_moments(dd)[el]) for x, dd in sorted(runs.items()) if el in site_moments(dd)]
            if not pts:
                continue
            p = np.array(pts)
            ax.errorbar(p[:, 0], p[:, 1], yerr=p[:, 2], fmt="D", color=cols[el], mfc=mfc, ms=4, lw=0.6, capsize=1.5,
                        label=f"QE SQS <m_{el}> ± site s.d., {lab}")
            dev = p[:, 1] - np.interp(p[:, 0], d[:, 0], d[:, 1])
            note("site moments", f"QE {lab} {el}", "korzhavyi2009_fig8 FM", "mean |QE - EMTO| (muB)", float(np.mean(np.abs(dev))))
            for (x, m, s), dv in zip(p, dev):
                note("site moments", f"QE {lab} {el} x={x:.3f}", "korzhavyi2009_fig8 FM", "QE - EMTO (muB)", float(dv))
    ax.axhline(0, color="0.7", lw=0.5)
    ax.set(xlim=(0, 1), ylim=(-1.8, 2.6), xlabel="$x_{Cr}$", ylabel="site moment (μB)", title="Fe and Cr moments")
    ax.legend(fontsize=4.6, loc="center right", bbox_to_anchor=(1.0, 0.42), ncol=1, framealpha=0.92)
    fig.tight_layout(); fig.savefig(OUT / "V_moments.png", dpi=300); plt.close(fig); print("wrote V_moments")


def plot_bulk(plt):
    sys.path.insert(0, str(HERE.parent / "scripts"))
    import make_figures as mf
    kz = rows("korzhavyi2009_fig3b_bulk_modulus.csv")
    fig, ax = plt.subplots(figsize=(4.8, 3.4))
    mk = {"Expt. RT FM": ("^", "0.05", "0.05"), "Expt. LT FM": ("^", "white", "0.05"), "Expt. RT AFM": ("D", "0.05", "0.05"),
          "Expt. LT AFM": ("D", "white", "0.05"), "Expt. RT PM": ("s", "#C00000", "#C00000")}
    for s, (m, fc, ec) in mk.items():
        d = np.array([(float(r["c_Cr_at_pct"]) / 100, float(r["B_GPa"])) for r in kz if r["series"] == s])
        if len(d):
            ax.plot(d[:, 0], d[:, 1], m, mfc=fc, mec=ec, ms=4.5, ls="none", label=f"{s} (via Korzhavyi 2009)")
    th = np.array(sorted((float(r["c_Cr_at_pct"]) / 100, float(r["B_GPa"])) for r in kz if r["series"] == "Theory FM"))
    ax.plot(th[:, 0], th[:, 1], "-", color="0.5", lw=0.9, label="EMTO-CPA FM (Korzhavyi 2009)")
    exp = {(0.0, 168.0), (1.0, 200.7)}           # RT FM Fe, RT PM Cr: the two end points used for the bias below
    for tag, col, mfc, lab in (("", "#D55E00", "white", "QE SQS, relaxed a"), ("aexp", "#0072B2", "#0072B2", "QE SQS, experimental a")):
        runs = sqs_runs(tag)
        p = np.array([(x, vrh_k(json.loads((d / "result.json").read_text())["C11"], json.loads((d / "result.json").read_text())["C12"]))
                      for x, d in sorted(runs.items())])
        if len(p):
            ax.plot(p[:, 0], p[:, 1], "D-", color=col, mfc=mfc, ms=4, lw=0.8, label=lab)
            pm = np.array([(float(r["c_Cr_at_pct"]) / 100, float(r["B_GPa"])) for r in kz if r["series"] == "Expt. RT PM"])
            ref95 = float(np.interp(0.95, pm[:, 0], pm[:, 1]))
            if 0.95 in runs or any(abs(x - 0.95) < 1e-6 for x in runs):
                v = float(p[np.argmin(abs(p[:, 0] - 0.95)), 1])
                note("bulk modulus", lab, "Expt. RT PM interpolated to x=0.95", "K(QE)/K(exp) - 1 (%)", 100 * (v / ref95 - 1))
    vt = mf.saw_tensors()
    v = np.array([(x, vrh_k(c11, c12)) for t, (x, c11, c12, c44) in vt.items() if t.startswith("dft_")])
    ax.plot(v[:, 0], v[:, 1], "o", color="#009E73", ms=4, label="VASP (this work)")
    note("bulk modulus", "VASP x=0", "Expt. RT FM Fe 168 GPa", "K/K_exp - 1 (%)", 100 * (float(v[v[:, 0] == 0, 1][0]) / 168.0 - 1))
    note("bulk modulus", "VASP x=1", "Expt. RT PM Cr 200.7 GPa", "K/K_exp - 1 (%)", 100 * (float(v[v[:, 0] == 1, 1][0]) / 200.7 - 1))
    ax.set(xlim=(-0.02, 1.02), ylim=(130, 300), xlabel="$x_{Cr}$", ylabel="bulk modulus K (GPa)", title="Bulk modulus vs experiment")
    ax.legend(fontsize=4.6, loc="upper left", ncol=2, framealpha=0.92)
    fig.tight_layout(); fig.savefig(OUT / "V_bulk_modulus.png", dpi=300); plt.close(fig); print("wrote V_bulk_modulus")


def sqs_nef(dos_root):
    import re
    out = []
    for d in sorted(dos_root.glob("x*")):
        if "_" in d.name or not (d / "nscf.out").exists():
            continue
        pre = next((p for p in ("pdosW", "pdosT") if (d / f"{p}.pdos_tot").exists()), None)
        if pre is None:
            continue
        ef = float(re.findall(r"the Fermi energy is\s+([-\d.]+) ev", (d / "nscf.out").read_text())[-1])
        a = np.loadtxt(d / f"{pre}.pdos_tot"); m = np.abs(a[:, 0] - ef) < 0.1
        out.append((float(d.name[1:]), (a[m, 1].mean() + a[m, 2].mean()) / 40, a[m, 2].mean() / 40))
    return np.array(out)


def plot_nef(plt):
    kz = rows("korzhavyi2009_fig3d_nef.csv")
    fig, ax = plt.subplots(figsize=(4.4, 3.0))
    for s, ls in (("Total FM", "-"), ("Spin down FM", "-."), ("Total DLM", "--")):
        d = np.array(sorted((float(r["c_Cr_at_pct"]) / 100, float(r["NEF_states_per_Ry"]) / RY_EV) for r in kz if r["series"] == s))
        ax.plot(d[:, 0], d[:, 1], ls, color="0.15" if "FM" in s else "#C00000", lw=1.0, label=f"EMTO-CPA {s} (Korzhavyi 2009)")
        if s == "Total FM":
            tot = d
    q = sqs_nef(ENG / "dos")
    if len(q):
        ax.plot(q[:, 0], q[:, 1], "D", color="#0072B2", ms=4.5, label="QE SQS total, relaxed a (±0.1 eV)")
        ax.plot(q[:, 0], q[:, 2], "v", color="#0072B2", mfc="white", ms=4, label="QE SQS spin down")
        for x, n, _ in q:
            note("N(E_F)", f"QE relaxed x={x:.2f}", "korzhavyi2009_fig3d Total FM", "QE/EMTO - 1 (%)",
                 100 * (n / float(np.interp(x, tot[:, 0], tot[:, 1])) - 1))
    ax.set(xlim=(0, 1), ylim=(0, 1.8), xlabel="$x_{Cr}$", ylabel="$N(E_F)$ (states/eV/atom)", title="Fermi-level DOS")
    ax.legend(fontsize=4.8, loc="upper right", framealpha=0.92)
    fig.tight_layout(); fig.savefig(OUT / "V_nef.png", dpi=300); plt.close(fig); print("wrote V_nef")


def plot_bands_fe(plt):
    sys.path.insert(0, str(HERE))
    import fermi_slices as fs
    import bands_path as bp
    kz = rows("korzhavyi2009_fig7_bands_fe.csv")
    fig, axs = plt.subplots(1, 2, figsize=(5.6, 3.0), sharey=True)
    _, dist, _ = bp.path_points(bp.PATH_BCC)
    m = dist <= 1.0 + 1e-9                                               # Gamma-H = first segment, |k| = dist
    for ax, (panel, spin, s) in zip(axs, (("a", "up", 0), ("d", "down", 1))):
        e = np.array([(float(r["x_2pi_over_a"]), float(r["E_minus_EF_Ry"]) * RY_EV) for r in kz if r["panel"] == panel])
        for case, col, ls, lab in (("Fe_fm", "#D55E00", "-", "QE, relaxed a"), ("Fe_fm_aexp", "#0072B2", "--", "QE, experimental a")):
            d = ENG / "bandpath" / case
            if not (d / "bands_data.xml").exists():
                continue
            k, ev, ef = fs.read_xml(d / "bands_data.xml")
            bands = ev[s][m] - ef
            dev = [min(abs(np.interp(x, dist[m], bands[:, b]) - y) for b in range(bands.shape[1])) for x, y in e]
            rms = float(np.sqrt(np.mean(np.square(dev))))
            note(f"Fe bands Gamma-H spin {spin}", case, "korzhavyi2009_fig7", "RMS nearest-band |dE| (eV)", rms)
            for b in range(bands.shape[1]):
                ax.plot(dist[m], bands[:, b], ls, color=col, lw=0.9, label=f"{lab}: RMS {rms * 1000:.0f} meV" if b == 0 else None)
        ax.plot(e[:, 0], e[:, 1], ".", color="0.1", ms=1.5, label="EMTO-CPA, 0% Cr (Korzhavyi 2009 Fig. 7)")
        ax.axhline(0, color="0.6", lw=0.5, ls=":")
        ax.set(xlim=(0, 1), ylim=(-0.2 * RY_EV, 0.2 * RY_EV), title=f"Fe Γ–H, spin {spin}")
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Γ", "H"])
        ax.legend(fontsize=4.6, loc="lower left", framealpha=0.92)
    axs[0].set_ylabel("$E-E_F$ (eV)")
    fig.tight_layout(); fig.savefig(OUT / "V_bands_Fe.png", dpi=300); plt.close(fig); print("wrote V_bands_Fe")


def electrons():
    plt = style()
    plot_moments(plt); plot_bulk(plt); plot_nef(plt); plot_bands_fe(plt)


def write_summary():
    if not SUMMARY:
        return
    with open(OUT / "validation_summary.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(SUMMARY[0]))
        w.writeheader(); w.writerows(SUMMARY)
    print("wrote validation_summary.csv")


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    what = sys.argv[1] if len(sys.argv) > 1 else "all"
    if what in ("phonons", "all"):
        phonons()
    if what in ("electrons", "all"):
        electrons()
    write_summary()
