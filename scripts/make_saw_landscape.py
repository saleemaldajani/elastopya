#!/usr/bin/env python3
"""SM figure S_saw_landscape: where in composition the (P)SAW of Fe-Cr reaches the band that pristine CASS cannot
produce, and what spinodal decomposition of Fe-36Cr does to that (added 2026-10-15).

(a) Highest (P)SAW of each of 10^4 grains (6.4 um) against x_Cr: 5-95 % range (bar), median (dot) and maximum (tick)
    for QE at the experimental lattice parameter, QE at the PBE-relaxed lattice parameter, and MD (25^3 cell, five
    alloys) random (phi = 0.5) and fully decomposed (phi = 1). The band above 609 MHz is shaded. Decomposition of
    Fe-36Cr moves material to the alpha (0.14) and alpha' (0.88) compositions (CALPHAD, 748 K).
(b) Fraction of grains in the band against x_Cr, and the lever-rule estimate for decomposed Fe-36Cr:
    f_band(dec) = (1 - w) f_band(alpha) + w f_band(alpha'),  w = (0.36 - 0.14) / (0.88 - 0.14).

  make_saw_landscape.py  -> manuscript/figures/kdense_scripta/supplementary/S_saw_landscape.{png,pdf,svg}
"""
import glob
import json
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_figures as mf  # noqa: E402

SAW = mf.DATA / "saw"
OUT = HERE.parent / "manuscript" / "figures" / "kdense_scripta" / "supplementary" / "S_saw_landscape"
EDGE = 609.0
X0, XA, XB = 0.36, 0.14, 0.88          # Fe-36Cr and the CALPHAD alpha / alpha' compositions at 748 K (Xiong 2011)
SERIES = (("qeaexp", "QE, exp. $a$ (random SQS)", "#0072B2", "D", -0.009),
          ("qe", "QE, PBE-relaxed (random SQS)", "#8F3A00", "o", -0.003),
          ("md5r", r"MD random ($\phi$ = 0.5)", "#009E73", "s", 0.003),
          ("md5d", r"MD decomposed ($\phi$ = 1)", "#0B5D3B", "s", 0.009))


def grain_top(path):
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    keep = d[:, 3] / d[:, 2] >= 2e-3                     # PSAW counted only where it carries surface intensity
    return np.where(keep, np.maximum(d[:, 0], d[:, 1]), d[:, 0]) / 6.4e-6 / 1e6


def load(prefix):
    out = []
    for p in sorted(glob.glob(str(SAW / f"{prefix}_x*.csv"))):
        x = float(re.search(r"_x([\d.]+)\.csv", p).group(1))
        g = grain_top(p)
        out.append((x, np.percentile(g, 5), np.median(g), np.percentile(g, 95), g.max(), 100 * np.mean(g > EDGE)))
    return np.array(sorted(out))


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, str(HERE / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    fig, (a, b) = plt.subplots(2, 1, figsize=(6.8, 5.6), sharex=True, gridspec_kw=dict(height_ratios=[1.6, 1], hspace=0.12))
    a.axhspan(EDGE, 860, color="0.90", lw=0); a.axhline(EDGE, color="0.4", lw=0.8, ls="--")
    a.text(0.62, 800, "unallowed in pristine CASS\n(above the highest (P)SAW of its grains)", ha="center", va="center", fontsize=6.5, color="0.35")
    data, summary = {}, {}
    for pre, lab, col, mk, dx in SERIES:
        s = load(pre)
        if not len(s):
            continue
        data[pre] = s
        x = s[:, 0] + dx
        a.vlines(x, s[:, 1], s[:, 3], color=col, lw=2.2, alpha=0.55)
        a.plot(x, s[:, 2], mk, color=col, ms=3.0, mec="white", mew=0.3, label=lab)
        a.plot(x, s[:, 4], "_", color=col, ms=5, mew=1.2)
        b.plot(s[:, 0], s[:, 5], mk + "-", color=col, ms=3.0, lw=0.9, label=lab)
        fa, fb, f0 = (float(np.interp(v, s[:, 0], s[:, 5])) for v in (XA, XB, X0))
        w = (X0 - XA) / (XB - XA)
        summary[pre] = dict(random_36=f0, alpha=fa, alpha_prime=fb, lever_decomposed=(1 - w) * fa + w * fb)
    for ax in (a, b):
        for v, t in ((XA, r"$\alpha$"), (X0, "Fe$-$36Cr"), (XB, r"$\alpha'$")):
            ax.axvline(v, color="#FF3B00" if v == X0 else "0.55", lw=0.6, ls=":")
    a.annotate("", xy=(XA, 372), xytext=(X0, 372), arrowprops=dict(arrowstyle="-|>", color="#FF3B00", lw=1.0))
    a.annotate("", xy=(XB, 372), xytext=(X0, 372), arrowprops=dict(arrowstyle="-|>", color="#FF3B00", lw=1.0))
    a.text(X0, 380, "decomposition (CALPHAD, 748 K)", ha="center", va="bottom", fontsize=6, color="#FF3B00")
    for v, t in ((XA, r"$\alpha$ (0.14)"), (X0, "Fe$-$36Cr"), (XB, r"$\alpha'$ (0.88)")):
        a.text(v + 0.006, 745, t, ha="left", va="center", fontsize=6.3, color="#FF3B00" if v == X0 else "0.3")
    a.set(ylabel="highest (P)SAW per grain (MHz)", ylim=(350, 860), xlim=(-0.03, 1.03))
    a.legend(fontsize=6, loc="lower center", bbox_to_anchor=(0.5, 1.0), frameon=False, ncol=4, handletextpad=0.2, columnspacing=1.0)
    a.text(0.01, 0.97, "(a)  bar: 5$-$95% of grains; marker: median; tick: maximum", transform=a.transAxes, fontsize=6.5, va="top")
    b.set(xlabel=r"$x_{\rm Cr}$", ylabel="grains in band (%)", ylim=(-1, None))
    if "qeaexp" in summary:
        q = summary["qeaexp"]
        b.annotate(f"decomposed Fe$-$36Cr (lever rule, QE exp. $a$): {q['lever_decomposed']:.1f}% vs random {q['random_36']:.1f}%",
                   xy=(X0, q["lever_decomposed"]), xytext=(0.47, 0.55), textcoords="axes fraction", fontsize=6, color="#0072B2",
                   arrowprops=dict(arrowstyle="-|>", color="#0072B2", lw=0.7))
        b.plot([X0], [q["lever_decomposed"]], "*", color="#0072B2", ms=8)
    b.text(0.01, 0.95, "(b)  fraction of the $10^4$ grains above 609 MHz", transform=b.transAxes, fontsize=6.5, va="top")
    fig.subplots_adjust(left=0.10, right=0.98, top=0.94, bottom=0.08)
    for ext in ("png", "pdf", "svg"):
        fig.savefig(OUT.with_suffix(f".{ext}"), dpi=600)
    (OUT.parent / "S_saw_landscape_summary.json").write_text(json.dumps(
        {"summary": summary, "per_composition": {k: v.round(2).tolist() for k, v in data.items()}}, indent=1))
    print("wrote", OUT)
    print(json.dumps({k: {kk: round(vv, 2) for kk, vv in v.items()} for k, v in summary.items()}, indent=1))


if __name__ == "__main__" and len(sys.argv) == 1:
    main()


# ---------------------------------------------------------------------------------------------------------------------
# v2 (2026-10-15, user request): the same data as split violins. Top: QE at the experimental (left half) and the
# PBE-relaxed (right half) lattice parameter, both random SQS. Middle: MD random (left) and decomposed (right), five
# MD alloys each. Bottom: fraction of grains in the band, as v1(b). v1 is kept.
# ---------------------------------------------------------------------------------------------------------------------
def all_grains(prefix):
    out = {}
    for p in sorted(glob.glob(str(SAW / f"{prefix}_x*.csv"))):
        out[float(re.search(r"_x([\d.]+)\.csv", p).group(1))] = grain_top(p)
    return out


def half_violin(ax, data, pos, side, col, width=0.034):
    v = ax.violinplot(data, positions=[pos], widths=width, showextrema=False, bw_method=0.15)
    for b in v["bodies"]:
        verts = b.get_paths()[0].vertices
        verts[:, 0] = np.clip(verts[:, 0], -np.inf, pos) if side == "left" else np.clip(verts[:, 0], pos, np.inf)
        b.set_facecolor(col); b.set_edgecolor(col); b.set_alpha(0.6); b.set_linewidth(0.4)
    m = data.max()
    ax.plot([pos - width / 2, pos] if side == "left" else [pos, pos + width / 2], [m, m], color=col, lw=1.0)


def main_v2():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    sys.path.insert(0, str(HERE / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    fig, (a1, a2, b) = plt.subplots(3, 1, figsize=(6.8, 7.4), sharex=True,
                                    gridspec_kw=dict(height_ratios=[1.3, 1.3, 0.9], hspace=0.10))
    pairs = ((a1, ("qeaexp", "#0072B2", "QE, exp. $a$ (random SQS)"), ("qe", "#8F3A00", "QE, PBE-relaxed (random SQS)")),
             (a2, ("md5r", "#009E73", r"MD random ($\phi$ = 0.5)"), ("md5d", "#0B5D3B", r"MD decomposed ($\phi$ = 1)")))
    for ax, left, right in pairs:
        ax.axhspan(EDGE, 860, color="0.90", lw=0); ax.axhline(EDGE, color="0.4", lw=0.8, ls="--")
        for (pre, col, lab), side in ((left, "left"), (right, "right")):
            for x, g in all_grains(pre).items():
                half_violin(ax, g, x, side, col)
        ax.legend(handles=[Patch(color=left[1], alpha=0.6, label=f"left half: {left[2]}"),
                           Patch(color=right[1], alpha=0.6, label=f"right half: {right[2]}")],
                  fontsize=6, loc="upper left", frameon=False, ncol=2)
        for v, t in ((XA, r"$\alpha$ (0.14)"), (X0, "Fe$-$36Cr"), (XB, r"$\alpha'$ (0.88)")):
            ax.axvline(v, color="#FF3B00" if v == X0 else "0.55", lw=0.6, ls=":")
            ax.text(v + 0.006, 760, t, ha="left", va="center", fontsize=6, color="#FF3B00" if v == X0 else "0.3")
        ax.set(ylabel="highest (P)SAW per grain (MHz)", ylim=(350, 860), xlim=(-0.04, 1.04))
        ax.text(0.62, 0.80, "unallowed in pristine CASS", transform=ax.transAxes, ha="center", fontsize=6.3, color="0.35")
    a1.text(0.99, 0.95, "(a) QE: random SQS at two lattice parameters; tick: maximum", transform=a1.transAxes, fontsize=6.5, ha="right", va="top")
    a2.text(0.99, 0.95, "(b) MD: random vs fully decomposed, same composition", transform=a2.transAxes, fontsize=6.5, ha="right", va="top")
    a2.annotate("", xy=(XA, 700), xytext=(X0, 700), arrowprops=dict(arrowstyle="-|>", color="#FF3B00", lw=1.0))
    a2.annotate("", xy=(XB, 700), xytext=(X0, 700), arrowprops=dict(arrowstyle="-|>", color="#FF3B00", lw=1.0))
    a2.text(0.62, 708, "decomposition of Fe$-$36Cr (CALPHAD, 748 K)", ha="center", va="bottom", fontsize=6, color="#FF3B00")
    summary = {}
    for pre, lab, col, mk, _ in SERIES:
        s = load(pre)
        if not len(s):
            continue
        b.plot(s[:, 0], s[:, 5], mk + "-", color=col, ms=3.0, lw=0.9, label=lab)
        fa, fb, f0 = (float(np.interp(v, s[:, 0], s[:, 5])) for v in (XA, XB, X0))
        w = (X0 - XA) / (XB - XA)
        summary[pre] = dict(random_36=f0, lever_decomposed=(1 - w) * fa + w * fb)
    for v in (XA, X0, XB):
        b.axvline(v, color="#FF3B00" if v == X0 else "0.55", lw=0.6, ls=":")
    q = summary.get("qeaexp")
    if q:
        b.plot([X0], [q["lever_decomposed"]], "*", color="#0072B2", ms=8)
        b.annotate(f"decomposed Fe$-$36Cr (lever rule, QE exp. $a$): {q['lever_decomposed']:.1f}% vs random {q['random_36']:.1f}%",
                   xy=(X0, q["lever_decomposed"]), xytext=(0.36, 0.72), textcoords="axes fraction", fontsize=5.8, color="#0072B2",
                   arrowprops=dict(arrowstyle="-|>", color="#0072B2", lw=0.7))
    b.set(xlabel=r"$x_{\rm Cr}$", ylabel="grains in band (%)", ylim=(-1, None))
    b.legend(fontsize=5.6, loc="upper left", bbox_to_anchor=(0.0, 0.88), frameon=False)
    b.text(0.01, 0.95, "(c) fraction of the $10^4$ grains above 609 MHz", transform=b.transAxes, fontsize=6.5, va="top")
    fig.subplots_adjust(left=0.10, right=0.98, top=0.98, bottom=0.07)
    out = OUT.parent / "S_saw_landscape_v2"
    for ext in ("png", "pdf", "svg"):
        fig.savefig(out.with_suffix(f".{ext}"), dpi=600)
    print("wrote", out)


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "v2":
    main_v2()


# ---------------------------------------------------------------------------------------------------------------------
# v3 (2026-10-15, user request: "also show the decomposed violins"): random (left halves) vs decomposed (right halves)
# for every series. QE decomposed = effective lever-rule tensor of alpha (0.14) and alpha' (0.88), see
# run_saw_decomposed_qe.py (validated on MD, where the decomposed cell exists). v1, v2 kept.
# ---------------------------------------------------------------------------------------------------------------------
def main_v3():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    sys.path.insert(0, str(HERE / "kdense_scripta" / "scripts"))
    from style_presets import apply_publication_style
    apply_publication_style("scripta"); plt.rcParams["svg.fonttype"] = "none"
    fig, axs = plt.subplots(4, 1, figsize=(6.8, 9.0), sharex=True, gridspec_kw=dict(height_ratios=[1.2, 1.2, 1.2, 0.95], hspace=0.10))
    rows = (("qeaexp", "qeaexpdec", "#0072B2", "#003F6B", "(a) QE, experimental lattice parameter"),
            ("qe", "qedec", "#C0712F", "#7A3A00", "(b) QE, PBE-relaxed lattice parameter"),
            ("md5r", "md5d", "#3CC29A", "#0B5D3B", "(c) MD, 25$^3$ cell, five alloys"))
    for ax, (r, d, cr, cd, title) in zip(axs[:3], rows):
        ax.axhspan(EDGE, 860, color="0.90", lw=0); ax.axhline(EDGE, color="0.4", lw=0.8, ls="--")
        for pre, side, col in ((r, "left", cr), (d, "right", cd)):
            for x, g in all_grains(pre).items():
                half_violin(ax, g, x, side, col)
        ax.legend(handles=[Patch(color=cr, alpha=0.6, label="left: random solid solution"),
                           Patch(color=cd, alpha=0.6, label="right: decomposed into $\\alpha$ + $\\alpha'$" +
                                 (" (effective lever-rule tensor)" if r != "md5r" else " ($\\phi$ = 1 cell)"))],
                  fontsize=5.8, loc="upper left", frameon=False, ncol=2)
        ax.text(0.99, 0.95, title, transform=ax.transAxes, fontsize=6.5, ha="right", va="top")
        for v, t in ((XA, r"$\alpha$ (0.14)"), (X0, "Fe$-$36Cr"), (XB, r"$\alpha'$ (0.88)")):
            ax.axvline(v, color="#FF3B00" if v == X0 else "0.55", lw=0.6, ls=":")
            ax.text(v + 0.006, 760, t, ha="left", va="center", fontsize=6, color="#FF3B00" if v == X0 else "0.3")
        ax.text(0.62, 0.80, "unallowed in pristine CASS", transform=ax.transAxes, ha="center", fontsize=6.3, color="0.35")
        ax.set(ylabel="highest (P)SAW\nper grain (MHz)", ylim=(350, 860), xlim=(-0.04, 1.04))
    b = axs[3]
    res = {}
    for r, d, cr, cd, title in rows:
        for pre, col, ls, mk in ((r, cr, "-", "o"), (d, cd, "--", "s")):
            s = load(pre)
            if len(s):
                lab = {"qeaexp": "QE exp. $a$", "qe": "QE relaxed", "md5r": "MD"}[r] + (" random" if pre == r else " decomposed")
                b.plot(s[:, 0], s[:, 5], mk + ls, color=col, ms=2.8, lw=0.9, label=lab)
                res[pre] = float(np.interp(X0, s[:, 0], s[:, 5]))
    for v in (XA, X0, XB):
        b.axvline(v, color="#FF3B00" if v == X0 else "0.55", lw=0.6, ls=":")
    if "qeaexp" in res and "qeaexpdec" in res:
        b.text(0.99, 0.95, f"Fe$-$36Cr, QE exp. $a$: random {res['qeaexp']:.1f}%, decomposed {res['qeaexpdec']:.1f}% of grains in band",
               transform=b.transAxes, ha="right", va="top", fontsize=6, color="#003F6B")
    b.set(xlabel=r"$x_{\rm Cr}$", ylabel="grains in band (%)", ylim=(-1, None))
    b.legend(fontsize=5.2, loc="upper left", bbox_to_anchor=(0.0, 0.86), frameon=False, ncol=2)
    b.text(0.01, 0.95, "(d) fraction of the $10^4$ grains above 609 MHz", transform=b.transAxes, fontsize=6.5, va="top")
    fig.subplots_adjust(left=0.11, right=0.98, top=0.99, bottom=0.06)
    out = OUT.parent / "S_saw_landscape_v3"
    for ext in ("png", "pdf", "svg"):
        fig.savefig(out.with_suffix(f".{ext}"), dpi=600)
    print("wrote", out, json.dumps({k: round(v, 2) for k, v in res.items()}))


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "v3":
    main_v3()
