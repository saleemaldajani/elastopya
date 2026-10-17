#!/usr/bin/env python3
"""Comparison versions of Figs. 2 and 3 with both QE sets in panel (a) (2026-10-17, for the user's decision; the paper
figures are unchanged): PBE-relaxed (open symbols, as now) and experimental lattice parameter (filled blue diamonds).
The experimental-lattice 80% Cr tensor is left out: it lies on the antiferromagnetic branch (discussion_items.tex).

  make_fig23_both_qe.py  -> manuscript/figures/kdense_scripta/fig2_..._qe_qeboth.*, fig3_fsaw_mddecomp_6p4um_qe_qeboth.*
"""
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import make_figures_kdense as mk  # noqa: E402,F401  (Scripta style on make_figures)
from make_figures_kdense import mf  # noqa: E402

EXCLUDE = {0.80}   # AF-branch tensor


def qe_aexp():
    rows = {}
    for d in sorted((mf.ROOT / "dft_qe" / "engaging" / "runs").glob("sqs40_x*_k442_e60_aexp")):
        j = d / "result.json"
        if j.exists():
            r = json.loads(j.read_text())
            if round(r["x_cr"], 3) not in EXCLUDE:
                rows[round(r["x_cr"], 3)] = (r["C11"], r["C12"], r["C44"])
    x = np.array(sorted(rows))
    c11, c12, c44 = (np.array([rows[v][i] for v in x]) for i in range(3))
    b, cp = (c11 + 2 * c12) / 3, (c11 - c12) / 2
    g = ((2 * cp + 3 * c44) / 5 + 5 * (2 * cp) * c44 / (4 * c44 + 3 * (2 * cp))) / 2
    e = 9 * b * g / (3 * b + g)
    return dict(x=x, c11=c11, c12=c12, c44=c44, f=mf.rayleigh_f(b, e, x * mf.RHO_CR + (1 - x) * mf.RHO_FE))


def final_legend(ax, h, l, fid):
    """Chosen layout (user, 2026-10-17): square axes kept; one column at 6.25 pt, tight spacing; 0.65 em from the axes
    so it clears the inward ticks (3.5 pt); 0.16 em line spacing leaves a gap above the QE circles at 45-50% Cr in 3a."""
    ax.legend(h, l, fontsize=6.25, loc="upper left", frameon=False, handlelength=1.2, handletextpad=0.35,
              labelspacing=0.16, borderaxespad=0.65)


def main():
    lit, dft, mdd, mdr = mf.literature(), mf.read_csv("sheet_dft.csv"), mf.md(), mf.md_random()
    qe = mf.qe_dft()
    mf.QE2, mf.STEM_SUFFIX = qe_aexp(), "_qeboth"
    mf.LEGEND_A = final_legend
    print("exp. a QE:", mf.QE2["x"], np.round(mf.QE2["c11"], 1), np.round(mf.QE2["f"], 4))
    mf.fig2_orientation_bulk(lit, dft, mdd, mdr, kind="violin", with_decomp=True, qe=qe)
    mf.fig3(lit, dft, mdd, mdr, lam=6.4, stem="fig3_fsaw_mddecomp_6p4um_qe_qeboth", md_decomp_a=True, qe=qe)
    mf.QE2, mf.STEM_SUFFIX, mf.LEGEND_A = None, "", None
    print("wrote to", mf.FIG)


if __name__ == "__main__":
    main()


def _split(labels):
    dft = [i for i, l in enumerate(labels) if "MD" not in l]
    md = [i for i, l in enumerate(labels) if "MD" in l]
    return dft, md


def legend_options():
    """Four layouts for the panel-(a) legend with both QE sets (2026-10-17): A taller y axis; B two columns above the
    panel; C DFT entries upper left, MD entries lower right; D legend outside on the right."""
    lit, dft, mdd, mdr = mf.literature(), mf.read_csv("sheet_dft.csv"), mf.md(), mf.md_random()
    qe = mf.qe_dft()
    mf.QE2 = qe_aexp()

    def A(ax, h, l, fid):
        lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + (0.33 if fid == "fig2" else 0.30) * (hi - lo))
        ax.legend(h, l, frameon=False, fontsize=6.3, loc="upper left", ncol=2, columnspacing=0.8, handletextpad=0.4)

    def B(ax, h, l, fid):
        ax.legend(h, l, frameon=False, fontsize=6.2, loc="lower left", bbox_to_anchor=(0.0, 1.06), ncol=2,
                  columnspacing=0.8, handletextpad=0.4, borderaxespad=0.0)

    def C(ax, h, l, fid):
        d, m = _split(l)
        lg = ax.legend([h[i] for i in d], [l[i] for i in d], frameon=False, fontsize=6.3, loc="upper left")
        ax.add_artist(lg)
        if fid == "fig2":
            ax.legend([h[i] for i in m], [l[i] for i in m], frameon=False, fontsize=6.3, loc="lower right")
        else:   # the alpha/alpha' arrows occupy the bottom; MD entries go to the right, mid-height
            ax.legend([h[i] for i in m], [l[i] for i in m], frameon=False, fontsize=6.3, loc="center right",
                      bbox_to_anchor=(1.0, 0.36))

    def D(ax, h, l, fid):
        ax.legend(h, l, frameon=False, fontsize=6.2, loc="upper left", bbox_to_anchor=(1.01, 1.0), borderaxespad=0.0)

    for tag, fn in (("A", A), ("B", B), ("C", C), ("D", D)):
        mf.LEGEND_A, mf.STEM_SUFFIX = fn, f"_qeboth_leg{tag}"
        mf.fig2_orientation_bulk(lit, dft, mdd, mdr, kind="violin", with_decomp=True, qe=qe)
        mf.fig3(lit, dft, mdd, mdr, lam=6.4, stem=f"fig3_fsaw_mddecomp_6p4um_qe_qeboth_leg{tag}", md_decomp_a=True, qe=qe)
    mf.LEGEND_A, mf.QE2, mf.STEM_SUFFIX = None, None, ""


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "legends":
    legend_options()


def legend_options_compact():
    """Square axes kept (user preference, 2026-10-17); legend made compact instead: E one column 5.6 pt, tight spacing;
    F two columns 5.8 pt across the top."""
    lit, dft, mdd, mdr = mf.literature(), mf.read_csv("sheet_dft.csv"), mf.md(), mf.md_random()
    qe = mf.qe_dft()
    mf.QE2 = qe_aexp()
    kw = dict(frameon=False, handlelength=1.2, handletextpad=0.35, labelspacing=0.22, borderaxespad=0.25)

    def E(ax, h, l, fid):
        ax.legend(h, l, fontsize=5.6, loc="upper left", **kw)

    def F(ax, h, l, fid):
        ax.legend(h, l, fontsize=5.8, loc="upper left", ncol=2, columnspacing=0.7, **kw)

    for tag, fn in (("E", E), ("F", F)):
        mf.LEGEND_A, mf.STEM_SUFFIX = fn, f"_qeboth_leg{tag}"
        mf.fig2_orientation_bulk(lit, dft, mdd, mdr, kind="violin", with_decomp=True, qe=qe)
        mf.fig3(lit, dft, mdd, mdr, lam=6.4, stem=f"fig3_fsaw_mddecomp_6p4um_qe_qeboth_leg{tag}", md_decomp_a=True, qe=qe)
    mf.LEGEND_A, mf.QE2, mf.STEM_SUFFIX = None, None, ""


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "legends2":
    legend_options_compact()
