#!/usr/bin/env python3
"""Electronic topological transition (ETT, Lifshitz 1960) on the Fe-rich side of bcc Fe-Cr, and its imprint on the
elastic constants and the SAW frequency.

Theory. When a band extremum (here the majority-spin level at H, triply degenerate, 0.15 eV above E_F in FM Fe) crosses
the Fermi level, the electronic grand potential acquires a non-analytic term
    dOmega(z) = -c |z|^{5/2} theta(z),      z = E_F - E_c,
so the density of states of the new Fermi-surface sheet is dN(z) ~ z^{1/2} theta(z). An elastic constant is a second
strain derivative of the energy; its singular part is dominated by
    dC_ij = d^2 dOmega / d eps_i d eps_j  ~  -(15/4) c (dz/d eps_i)(dz/d eps_j) z^{1/2} theta(z)
(the band term C_band = -sum (d e_k/d eps)^2 delta(e_k - E_F) of the main text, evaluated for the new sheet). Alloying
shifts z roughly linearly, z = kappa (x - x_c), so near the transition every modulus behaves as
    Q(x) = a + b x + A g(x),   g(x) = Re sqrt(x - x_c + i Gamma) - Re sqrt(-x_c + i Gamma),
a smooth background plus a square-root anomaly whose onset x_c is the same for every quantity; the substitutional
disorder of the alloy broadens it (Gamma). The sign and size of A differ between moduli.

What this script does
  1. x_c from the electronic structure: z(x) of the majority H level, read from our band structures (pure Fe, both
     lattice parameters) and the unfolded SQS bands (Fe-Cr), linearly interpolated to z = 0.
  2. x_c and A from fits of Q(x) to the EMTO-CPA moduli (Razumovskiy 2011), the measured bulk moduli (compiled by
     Korzhavyi 2009) and our QE SQS series; consistency of x_c across them is the test of the description.
  3. Link to stiffening on decomposition: the Fe-rich alpha phase (x ~ 0.12-0.14) lies beyond x_c, so its f_SAW carries
     the ETT term A g(x_alpha). Its share of the lever-rule gain of Fe-36Cr is (1 - w) A g(x_alpha) / f(x0).

  ett_analysis.py            -> prints the numbers, writes manuscript/.../supplementary/S_ett.{png,pdf,svg} and
                                dft_qe/ett_results.json
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(HERE))
import make_figures as mf  # noqa: E402

ENG = HERE / "engaging"
OUT = ROOT / "manuscript" / "figures" / "kdense_scripta" / "supplementary"
XMAX = 0.20                        # Fe-rich fitting window


def g(x, xc, gam):
    x = np.asarray(x, float)
    s = lambda u: np.real(np.sqrt(u + 1j * gam))
    return s(x - xc) - s(-xc)


def model(x, a, b, A, xc, gam):
    return a + b * np.asarray(x, float) + A * g(x, xc, gam)


def fit(x, y, label, gam_bounds=(1e-4, 0.03)):
    x, y = np.asarray(x, float), np.asarray(y, float)
    best = None
    for xc0 in np.arange(0.02, 0.16, 0.01):                 # multistart on the onset
        for A0 in (-50, 50):
            try:
                p, cov = curve_fit(model, x, y, p0=[y[0], 0.0, A0, xc0, 0.005],
                                   bounds=([-1e4, -1e4, -1e4, 0.0, gam_bounds[0]], [1e4, 1e4, 1e4, XMAX, gam_bounds[1]]),
                                   maxfev=20000)
            except RuntimeError:
                continue
            r = y - model(x, *p)
            sse = float(r @ r)
            if best is None or sse < best[2]:
                best = (p, cov, sse)
    p, cov, sse = best
    dof = max(len(x) - 5, 1)
    err = np.sqrt(np.clip(np.diag(cov), 0, None))
    return dict(label=label, a=p[0], b=p[1], A=p[2], xc=p[3], gamma=p[4], xc_err=err[3], A_err=err[2],
                rms=float(np.sqrt(sse / len(x))), n=len(x))


def vrh(c11, c12, c44):
    b = (c11 + 2 * c12) / 3; cp = (c11 - c12) / 2
    gg = ((2 * cp + 3 * c44) / 5 + 5 * (2 * cp) * c44 / (4 * c44 + 3 * (2 * cp))) / 2
    return b, gg


def f_saw(c11, c12, c44, x):
    b, gg = vrh(c11, c12, c44)
    e = 9 * b * gg / (3 * b + gg)
    return float(mf.rayleigh_f(b, e, x * mf.RHO_CR + (1 - x) * mf.RHO_FE)) * 1e3     # MHz


def qe_series(suffix):
    """(x, C44, B, f) of the QE pure Fe (2-atom) + SQS series in the Fe-rich window."""
    rows = []
    fe = ENG / "runs" / f"Fe_fm_k16x16x16_e60{suffix}" / "result.json"
    if fe.exists():
        r = json.loads(fe.read_text()); rows.append((0.0, r["C11"], r["C12"], r["C44"]))
    for d in sorted((ENG / "runs").glob(f"sqs40_x*_k442_e60{suffix}")):
        if suffix == "" and d.name.split("e60")[1]:
            continue
        j = d / "result.json"
        if j.exists():
            r = json.loads(j.read_text())
            if r["x_cr"] <= XMAX + 1e-9:
                rows.append((r["x_cr"], r["C11"], r["C12"], r["C44"]))
    rows.sort()
    return np.array([(x, c44, vrh(c11, c12, c44)[0], f_saw(c11, c12, c44, x)) for x, c11, c12, c44 in rows])


def critical_level():
    """z = E_F - E_c (eV) of the majority H level: pure Fe from the primitive-cell bands, alloys from unfolded bands."""
    import fermi_slices as fs
    import bands_path as bp
    out = {}
    _, dist, _ = bp.path_points(bp.PATH_BCC)
    iH = int(np.argmin(abs(dist - 1.0)))
    for case, key in (("Fe_fm", "relaxed"), ("Fe_fm_aexp", "aexp")):
        k, ev, ef = fs.read_xml(ENG / "bandpath" / case / "bands_data.xml")
        lev = np.sort(ev[0, iH] - ef)
        ec = float(lev[(lev > -0.05) & (lev < 0.6)][0])
        out.setdefault(key, []).append((0.0, -ec, 0.0))
    for d in sorted((ENG / "bandpath").glob("sqs_x*")):
        f = d / "ebs_up.npz"
        if not f.exists() or d.name.startswith("sqs_x0.00"):   # pure Fe comes from the primitive cell above;
            continue                                             # its 40-atom copy would count x = 0 twice
        u = np.load(f)["bands"]
        z = np.load(d / "pathmap.npz", allow_pickle=True)
        iH2 = int(np.argmin(abs(np.asarray(z["dist"]) - 1.0)))
        m = u[:, 0].astype(int) == iH2
        e, w = u[m, 2], u[m, 3]
        sel = (e > -0.3) & (e < 0.3) & (w > 0.05)               # the critical level, smeared by the disorder
        if sel.any():
            ec = float(np.average(e[sel], weights=w[sel]))
            spread = float(np.sqrt(np.average((e[sel] - ec) ** 2, weights=w[sel])))
            key = "aexp" if d.name.endswith("_aexp") else "relaxed"
            out.setdefault(key, []).append((float(z["x_cr"]), -ec, spread))
    return {k: sorted(v) for k, v in out.items()}


def main():
    res = {"fits": [], "critical_level": {}, "xc_electronic": {}, "stiffening": {}}
    lit = list(csv.DictReader(open(ROOT / "data" / "sheet_raw_razumovskiy.csv")))
    L = np.array([[float(r["x_cr"]), float(r["c44_gpa"]), float(r["bulk_gpa"]), float(r["cprime_gpa"])] for r in lit])
    L = L[L[:, 0] <= XMAX + 1e-9]
    lf = mf.literature(); xl, fl = np.asarray(lf["x"]), np.asarray(lf["f"]) * 1e3
    sel = xl <= XMAX + 1e-9
    for y, lab in ((L[:, 1], "EMTO-CPA C44"), (L[:, 2], "EMTO-CPA B"), (L[:, 3], "EMTO-CPA C'")):
        res["fits"].append(fit(L[:, 0], y, lab))
    res["fits"].append(fit(xl[sel], fl[sel], "EMTO-CPA f_SAW"))
    # measured bulk moduli (room temperature, FM), as compiled by Korzhavyi 2009 (digitized)
    ex = [r for r in csv.DictReader(l for l in open(HERE / "refdata" / "korzhavyi2009_fig3b_bulk_modulus.csv") if not l.startswith("#"))
          if r["series"] == "Expt. RT FM"]
    E = np.array(sorted((float(r["c_Cr_at_pct"]) / 100, float(r["B_GPa"])) for r in ex))
    # five points: the four-parameter fit is underdetermined, so fix Gamma and the slope to the EMTO-CPA B fit
    fb = next(f for f in res["fits"] if f["label"] == "EMTO-CPA B")
    pe, ce = curve_fit(lambda x, a, A, xc: model(x, a, fb["b"], A, xc, fb["gamma"]), E[:, 0], E[:, 1],
                       p0=[E[0, 1], fb["A"], fb["xc"]], bounds=([0, -1e4, 0], [1e4, 1e4, XMAX]))
    res["fits"].append(dict(label="measured B (RT FM)", a=pe[0], b=fb["b"], A=pe[1], xc=pe[2], gamma=fb["gamma"],
                            xc_err=float(np.sqrt(ce[2, 2])), A_err=float(np.sqrt(ce[1, 1])),
                            rms=float(np.sqrt(np.mean((E[:, 1] - model(E[:, 0], pe[0], fb["b"], pe[1], pe[2], fb["gamma"])) ** 2))), n=len(E)))
    for suf, tag in (("", "relaxed"), ("_aexp", "exp. a")):
        q = qe_series(suf)
        if len(q) >= 4:
            # 4-5 points: fit the onset with the slope b fixed to zero and Gamma to the EMTO-CPA C44 value
            fc = next(f for f in res["fits"] if f["label"] == "EMTO-CPA C44")
            for col, nm in ((1, "C44"), (3, "f_SAW")):
                p, c = curve_fit(lambda x, a, A, xc: model(x, a, 0.0, A, xc, fc["gamma"]), q[:, 0], q[:, col],
                                 p0=[q[0, col], 30, 0.05], bounds=([0, -1e4, 0], [1e4, 1e4, XMAX]))
                res["fits"].append(dict(label=f"QE {tag} {nm}", a=p[0], b=0.0, A=p[1], xc=p[2], gamma=fc["gamma"],
                                        xc_err=float(np.sqrt(c[2, 2])), A_err=float(np.sqrt(c[1, 1])), n=len(q),
                                        rms=float(np.sqrt(np.mean((q[:, col] - model(q[:, 0], p[0], 0, p[1], p[2], fc["gamma"])) ** 2)))))
    # electronic x_c
    zc = critical_level()
    res["critical_level"] = zc
    for k, pts in zc.items():
        if len(pts) >= 2:
            x, z, sp = np.array(pts).T
            s, i0 = np.polyfit(x, z, 1)
            xc = float(-i0 / s)
            res["xc_electronic"][k] = dict(xc=xc, dz_dx_eV=float(s), level_spread_eV=float(sp.max()),
                                           xc_window=float(sp.max() / s))      # disorder width converted to composition
    # joint fit of the EMTO-CPA moduli with one shared onset x_c and broadening Gamma (separate a, b, A per modulus)
    sets = [(L[:, 0], L[:, 1], "C44"), (L[:, 0], L[:, 2], "B"), (L[:, 0], L[:, 3], "C'"), (xl[sel], fl[sel], "f_SAW")]
    def joint_sse(xc, gam):
        tot = 0.0; per = {}
        for x, y, nm in sets:
            G = np.column_stack([np.ones_like(x), x, g(x, xc, gam)])     # linear in a, b, A for fixed (xc, Gamma)
            c, *_ = np.linalg.lstsq(G, y, rcond=None)
            r = y - G @ c; tot += float(r @ r); per[nm] = (c, float(np.sqrt(np.mean(r ** 2))))
        return tot, per
    xcs, gams = np.linspace(0.0, 0.18, 181), np.geomspace(1e-3, 0.1, 41)
    S = np.array([[joint_sse(a, b)[0] for b in gams] for a in xcs])
    i, j = np.unravel_index(np.argmin(S), S.shape)
    prof = S.min(axis=1)
    n_dat = sum(len(x) for x, _, _ in sets); n_par = 3 * len(sets) + 2
    s2 = S[i, j] / (n_dat - n_par)
    ok = xcs[(prof - S[i, j]) / s2 <= 1.0]                                 # Delta chi^2 <= 1 (profile, scaled)
    best_sse, best_per = joint_sse(xcs[i], gams[j])
    res["joint_fit"] = dict(xc=float(xcs[i]), gamma=float(gams[j]), xc_lo=float(ok.min()), xc_hi=float(ok.max()),
                            rms={k: v[1] for k, v in best_per.items()}, A={k: float(v[0][2]) for k, v in best_per.items()})
    res["joint_profile"] = dict(xc=xcs.tolist(), dchi2=((prof - S[i, j]) / s2).tolist())
    # prediction test: onset and width fixed to the electronic-structure values (experimental lattice parameter)
    ee = res["xc_electronic"].get("aexp")
    if ee:
        pred_sse, pred_per = joint_sse(ee["xc"], max(ee["xc_window"], 1e-3))
        res["prediction_test"] = dict(xc=ee["xc"], gamma=ee["xc_window"], rms={k: v[1] for k, v in pred_per.items()},
                                      dchi2_vs_best=float((pred_sse - S[i, j]) / s2))
    # model-independent anomaly locations
    res["anomaly_locations"] = {
        "EMTO-CPA B minimum": float(L[np.argmin(L[:, 2]), 0]),
        "measured B minimum (RT FM)": float(E[np.argmin(E[:, 1]), 0]),
        "EMTO-CPA C44 knee (max curvature)": float(L[1 + np.argmin(np.diff(L[:, 1], 2) / np.diff(L[:, 0])[1:]), 0]),
    }
    # stiffening link: split of the lever-rule gain of Fe-36Cr into the alpha and alpha' terms, electronic vs classical
    md = mf.md_random(); xm, fm = np.asarray(md["x"]), np.asarray(md["f"]) * 1e3
    for name, (xa, xb) in {"MD phases": (0.12, 0.74), "CALPHAD 748 K": (0.14, 0.88)}.items():
        w = (0.36 - xa) / (xb - xa)
        row = {}
        for lab, (X, F) in (("electronic (EMTO-CPA)", (xl, fl)), ("classical (MD/EAM)", (xm, fm))):
            f0, fa, fb2 = (float(np.interp(v, X, F)) for v in (0.36, xa, xb))
            row[lab] = dict(alpha_term_pct=100 * (1 - w) * (fa - f0) / f0, alpha_prime_term_pct=100 * w * (fb2 - f0) / f0,
                            total_pct=100 * ((1 - w) * fa + w * fb2 - f0) / f0)
        res["lever_split"] = res.get("lever_split", {}); res["lever_split"][name] = row
    # (superseded first-pass estimate kept below for the record; it extrapolates the sqrt term beyond its validity)
    # stiffening link: share of the electronic lever-rule gain carried by the ETT term at alpha
    ff = next(f for f in res["fits"] if f["label"] == "EMTO-CPA f_SAW")
    x0 = 0.36
    f0 = float(np.interp(x0, xl, fl))
    for name, (xa, xb) in {"MD phases": (0.12, 0.74), "CALPHAD 748 K": (0.14, 0.88)}.items():
        w = (x0 - xa) / (xb - xa)
        gain = ((1 - w) * np.interp(xa, xl, fl) + w * np.interp(xb, xl, fl)) / f0 - 1
        ett = (1 - w) * ff["A"] * float(g(xa, ff["xc"], ff["gamma"])) / f0
        res["stiffening"][name] = dict(alpha=xa, alpha_prime=xb, w_alpha_prime=w, lever_gain_pct=100 * gain,
                                       ett_share_pct=100 * ett, ett_fraction_of_gain=ett / gain)
    for f in res["fits"]:
        print(f"{f['label']:24s} x_c = {f['xc']:.3f} ± {f['xc_err']:.3f}   A = {f['A']:+8.1f} ± {f['A_err']:.1f}   Γ = {f['gamma']:.4f}   rms = {f['rms']:.2f} (n={f['n']})")
    print("critical level z = E_F - E_c (eV):", zc)
    print("x_c from electronic structure:", res["xc_electronic"])
    print("joint fit:", json.dumps(res["joint_fit"], indent=1))
    print("prediction test:", json.dumps(res.get("prediction_test"), indent=1))
    print("anomaly locations:", res["anomaly_locations"])
    print("lever split:", json.dumps(res["lever_split"], indent=1))
    (HERE / "ett_results.json").write_text(json.dumps(res, indent=1, default=float))
    return res, L, (xl, fl), E


if __name__ == "__main__" and "--plot" not in sys.argv:
    main()


def plot(res, L, lit_f, E):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(ROOT / "scripts" / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    xl, fl = lit_f
    jf, ee = res["joint_fit"], res["xc_electronic"].get("aexp")
    fig, axs = plt.subplots(2, 3, figsize=(7.4, 4.6))
    win = lambda ax: (ax.axvspan(max(ee["xc"] - ee["xc_window"], 0), ee["xc"] + ee["xc_window"], color="#FFB000", alpha=0.10, lw=0),
                      ax.axvspan(jf["xc_lo"], jf["xc_hi"], color="0.55", alpha=0.25, lw=0),
                      ax.axvline(ee["xc"], color="#B07800", lw=0.8, ls="--"))
    # (a) critical level
    ax = axs[0, 0]; win(ax)
    xx = np.linspace(0, 0.2, 50)
    for key, mk, col, lab, dx, ls in (("relaxed", "o", "#D55E00", "PBE-relaxed $a$", -0.003, "--"),
                                      ("aexp", "D", "#0072B2", "experimental $a$", 0.003, "-")):
        pts = np.array(res["critical_level"].get(key, []))
        if len(pts):   # the two series are offset by +-0.003 in x so that the error bars do not hide each other
            ax.errorbar(pts[:, 0] + dx, pts[:, 1], yerr=pts[:, 2], fmt=mk, color=col, ms=4, capsize=2, lw=0.8, label=lab)
        e_k = res["xc_electronic"].get(key)
        if e_k:
            ax.plot(xx, e_k["dz_dx_eV"] * (xx - e_k["xc"]), ls, color=col, lw=0.8)
    ax.axhline(0, color="0.5", lw=0.5)
    ax.set(xlim=(0, 0.2), xlabel="$x_{Cr}$", ylabel="$z=E_F-E_c$ (eV)", title="(a) majority level at H")
    ax.legend(fontsize=5.5, loc="upper left", frameon=False)
    # (b-d) moduli with the joint Lifshitz fit
    def fitline(ax, x, y, nm, col, lab):
        G = np.column_stack([np.ones_like(x), x, g(x, jf["xc"], jf["gamma"])]); c, *_ = np.linalg.lstsq(G, y, rcond=None)
        xs = np.linspace(0, 0.2, 200)
        ax.plot(x, y, "o", color=col, ms=3, label=lab)
        ax.plot(xs, c[0] + c[1] * xs + c[2] * g(xs, jf["xc"], jf["gamma"]), "-", color=col, lw=1.0)
    for ax, col_i, nm, ylab in ((axs[0, 1], 1, "C44", "$C_{44}$ (GPa)"), (axs[0, 2], 2, "B", "$B$ (GPa)")):
        win(ax); fitline(ax, L[:, 0], L[:, col_i], nm, "0.15", "EMTO-CPA")
        ax.set(xlim=(0, 0.2), xlabel="$x_{Cr}$", ylabel=ylab, title=f"({'b' if nm == 'C44' else 'c'}) {ylab.split(' (')[0]}")
    axs[0, 2].plot(E[:, 0], E[:, 1], "^", color="#C00000", ms=4, label="measured, RT FM")
    for suf, mk, col, lab in (("", "o", "#D55E00", "QE, relaxed $a$"), ("_aexp", "D", "#0072B2", "QE, exp. $a$")):
        q = qe_series(suf)
        if len(q):
            axs[0, 1].plot(q[:, 0], q[:, 1], mk, mfc="white", mec=col, ms=4, label=lab)
    axs[0, 1].legend(fontsize=5.2, loc="lower right", frameon=False); axs[0, 2].legend(fontsize=5.2, loc="upper right", frameon=False)
    ax = axs[1, 0]; win(ax); sel = xl <= XMAX + 1e-9
    fitline(ax, xl[sel], fl[sel], "f", "0.15", "EMTO-CPA")
    for suf, mk, col, lab in (("", "o", "#D55E00", "QE, relaxed $a$"), ("_aexp", "D", "#0072B2", "QE, exp. $a$")):
        q = qe_series(suf)
        if len(q):
            ax.plot(q[:, 0], q[:, 3], mk, mfc="white", mec=col, ms=4, label=lab)
    ax.set(xlim=(0, 0.2), xlabel="$x_{Cr}$", ylabel="$f_{SAW}$ (MHz, 6.4 µm)", title="(d) SAW frequency")
    ax.legend(fontsize=5.2, loc="lower right", frameon=False)
    # (e) profile of the joint fit
    ax = axs[1, 1]; pr = res["joint_profile"]; win(ax)
    ax.plot(pr["xc"], pr["dchi2"], color="0.15", lw=1.0); ax.axhline(1, color="0.5", lw=0.5, ls="--")
    for lab, v in res["anomaly_locations"].items():
        ax.plot([v], [0.3], "v", ms=5, color="#C00000" if "measured" in lab else "0.4")
    ax.set(xlim=(0, 0.18), ylim=(0, 12), xlabel="onset $x_c$", ylabel=r"$\Delta\chi^2$ (joint fit)", title="(e) where the moduli put $x_c$")
    # (f) lever-rule split
    ax = axs[1, 2]; ls = res["lever_split"]["CALPHAD 748 K"]
    labs = list(ls); xpos = np.arange(len(labs))
    a_ = [ls[k]["alpha_term_pct"] for k in labs]; b_ = [ls[k]["alpha_prime_term_pct"] for k in labs]
    ax.bar(xpos, a_, 0.55, color="#FFB000", label=r"$\alpha$ (Fe-rich, beyond $x_c$)")
    ax.bar(xpos, b_, 0.55, bottom=a_, color="#0072B2", label=r"$\alpha'$ (Cr-rich)")
    for i, k in enumerate(labs):
        ax.text(i, a_[i] + b_[i] + 0.15, f"{a_[i] + b_[i]:.1f}%", ha="center", fontsize=6)
    ax.set_xticks(xpos); ax.set_xticklabels(["electronic\n(EMTO-CPA)", "classical\n(MD/EAM)"], fontsize=6)
    ax.set(ylabel="lever-rule gain, Fe$-$36Cr (%)", ylim=(0, 9.5), title="(f) stiffening on decomposition")
    ax.legend(fontsize=5.2, loc="upper right", frameon=False)
    for a in axs.flat:
        t = a.get_title(); a.set_title(""); a.set_title(t, loc="left", fontsize=7)
    fig.tight_layout()
    out = OUT / "S_ett"
    fig.savefig(f"{out}.png", dpi=600); fig.savefig(f"{out}.pdf"); fig.savefig(f"{out}.svg")
    print("wrote", out)


if __name__ == "__main__" and "--plot" in sys.argv:
    plot(*main())
