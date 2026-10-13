#!/usr/bin/env python3
"""Supplementary figures on the electronic origin of the Fe-rich elastic anomaly (K-Dense Scripta style).

  S_fe_rich_elastic  B, C44 and C' of bcc Fe(1-x)Cr(x), 0 <= x <= 0.2: EMTO-CPA (Razumovskiy et al. 2011,
                     digitized), present DFT (VASP 40-atom SQS; QE 40-atom SQS where available) and MD (EAM)
                     random alloys on a 1 at.% grid. The shaded band marks the 8-12 at.% Cr range in which the
                     spin-up Fermi surface changes topology (Korzhavyi et al. 2009; Razumovskiy et al. 2011).
Output: manuscript/figures/kdense_scripta/supplementary/
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_figures_kdense as mk  # noqa: E402  (Scripta palette and checks)
from make_figures_kdense import mf, plt  # noqa: E402
from style_presets import apply_publication_style  # noqa: E402

OUT = mf.ROOT / "manuscript" / "figures" / "kdense_scripta" / "supplementary"
QE_RUNS = mf.ROOT / "dft_qe" / "engaging" / "runs"


def moduli(c11, c12, c44):
    return dict(B=(c11 + 2 * c12) / 3, C44=c44, Cp=(c11 - c12) / 2)


def qe_sqs():
    """QE 40-atom SQS results from Engaging (result.json), if present."""
    rows = []
    for p in sorted(QE_RUNS.glob("sqs40_x*_k442_e60/result.json")):
        r = json.loads(p.read_text())
        rows.append((r["x_cr"], r["C11"], r["C12"], r["C44"]))
    for p in sorted(QE_RUNS.glob("Fe_fm_k8x8x10_e60/result.json")):
        r = json.loads(p.read_text())
        rows.append((0.0, r["C11"], r["C12"], r["C44"]))
    return np.array(sorted(rows)) if rows else None


def fe_rich_elastic():
    lit = mf.read_csv("sheet_raw_razumovskiy.csv")
    lit_m = dict(B=lit["bulk_gpa"], C44=lit["c44_gpa"], Cp=lit["cprime_gpa"])
    md = mf.read_csv("md_random_alloys.csv")
    md_m = moduli(md["c11"], md["c12"], md["c44"])
    md_sd = dict(B=np.hypot(md["c11_sd"], 2 * md["c12_sd"]) / 3, C44=md["c44_sd"],
                 Cp=np.hypot(md["c11_sd"], md["c12_sd"]) / 2)
    dft = mf.read_csv("sheet_dft.csv")
    dft_m = moduli(dft["c11_gpa"], dft["c12_gpa"], dft["c44_gpa"])
    qe = qe_sqs()
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.5))
    labels = {"B": "$B$ (GPa)", "C44": "$C_{44}$ (GPa)", "Cp": "$C'$ (GPa)"}
    for ax, key, tag in zip(axs, ("B", "C44", "Cp"), ("(a)", "(b)", "(c)")):
        ax.axvspan(0.08, 0.12, color="0.9", lw=0, zorder=0)
        m = lit["x_cr"] <= 0.21
        ax.plot(lit["x_cr"][m], lit_m[key][m], "-o", color="0.15", ms=3.5, lw=0.9, label="EMTO-CPA [8]")
        mm = md["x_cr"] <= 0.21
        ax.errorbar(md["x_cr"][mm], md_m[key][mm], yerr=md_sd[key][mm], fmt="-s", color=mf.GREEN, ms=3.2,
                    lw=0.9, capsize=1.5, label="MD (EAM), random alloys")
        md_ = dft["x_cr"] <= 0.21
        ax.plot(dft["x_cr"][md_], dft_m[key][md_], "D", color=mf.RED, mec="white", mew=0.5, ms=5.5,
                label="DFT (VASP, this work)")
        if qe is not None:
            qm = qe[:, 0] <= 0.21
            q = moduli(qe[qm, 1], qe[qm, 2], qe[qm, 3])
            ax.plot(qe[qm, 0], q[key], "D", mfc="white", mec=mf.RED, mew=1.1, ms=5.5, label="DFT (QE, this work)")
        ax.set(xlabel=r"$x_{\rm Cr}$", ylabel=labels[key], xlim=(-0.01, 0.21))
        mk.style(ax, tag)
    axs[0].text(0.10, 0.97, "Fermi-surface\ntopology change", transform=axs[0].get_xaxis_transform(),
                ha="center", va="top", color="0.35")
    h, l = axs[2].get_legend_handles_labels()
    fig.legend(h, l, loc="upper center", ncol=len(l), frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.set_size_inches(mk.WIDTH_IN, 2.9)
    mk.normalize_text(fig)
    fig.tight_layout(rect=(0, 0, 1, 0.9), w_pad=1.2)
    mk.export_figure(fig, OUT / "S_fe_rich_elastic", formats=["pdf", "png", "svg"], dpi=mk.DPI, bbox_inches=None,
                     overwrite=True, mkdir=True, write_manifest=True,
                     provenance={"source": "scripts/make_si_figures.py"})
    found, sizes = mk.audit(fig)
    print("S_fe_rich_elastic", "min pt", min(sizes), "issues", found)
    plt.close(fig)


if __name__ == "__main__":
    apply_publication_style("scripta", palette_name="okabe_ito_on_white")
    fe_rich_elastic()
