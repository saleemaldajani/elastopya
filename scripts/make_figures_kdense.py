#!/usr/bin/env python3
"""Scripta Materialia versions of every paper figure, built with the K-Dense helpers.

Same data and panels as scripts/make_figures.py; only the presentation changes:
  - drawn at the final Elsevier full width (190 mm), exported without bbox cropping, so the
    PDF/TIFF is placed at 100 % with no rescaling;
  - K-Dense 'scripta' preset (scripts/kdense_scripta): Arial 7 pt lettering including mathtext,
    7.5 pt axis labels, 9 pt bold panel labels, closed frame with inward ticks;
  - Okabe-Ito colour-blind-safe palette in the same hue families as the draft figures;
  - every text item is checked at final size against other text and against plotted data; the
    result is written to <stem>.audit.json and printed;
  - strip versions (panel c: atomic snapshots, red Fe / blue Cr as in Acta Fig. 7) are drawn
    natively in the same figure, so they stay vector.

Output: manuscript/figures/kdense_scripta/{<stem>.pdf,.png,(.tiff),.export.json,.audit.json}
        manuscript/figures/kdense_scripta/with_strips/<stem>_strip.*
"""
import json
import sys
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import PathCollection  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.text import Text  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "kdense_scripta" / "scripts"))
import make_figures as mf  # noqa: E402
import figures_with_strips as fws  # noqa: E402
from figure_export import check_figure_size, export_figure  # noqa: E402
from style_presets import apply_publication_style, load_publisher_profiles  # noqa: E402
from visualize_decomposition import CR, FE, series  # noqa: E402

OUT = mf.ROOT / "manuscript" / "figures" / "kdense_scripta"
PROFILE = load_publisher_profiles()["profiles"]["scripta"]
WIDTH_IN = PROFILE["widths_mm"]["full"] / 25.4
MIN_PT = PROFILE["lettering_pt"]["normal"]
DPI = 600
TIFF_STEMS = {"fig1_free_energy", "fig2_elastic_mapping", "fig3_fsaw"}   # manuscript figures

# Okabe-Ito (Wong 2011) replacements, same hue family as the draft's Google-Sheets colours
mf.BLUE, mf.RED, mf.YELLOW, mf.GREEN, mf.ORANGE = "#0072B2", "#D55E00", "#E69F00", "#009E73", "#CC79A7"
mf.PURPLE = "#D55E00"   # Fe-36Cr in Fig. 3b (vermillion, distinct from Fe-20Cr blue)
mf.FIG = OUT
AUDIT = {}
STRIPS = {}   # stem -> (title, panels); filled in main()


def style(ax, label):
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True)
    ax.set_title(label, loc="left", fontsize=9, fontweight="bold", pad=3)


def normalize_text(fig):
    """Final-size lettering: Arial, nothing under 7 pt, axis labels 7.5 pt, panel labels 9 pt."""
    for ax in fig.axes:
        ax.tick_params(axis="both", which="both", labelsize=MIN_PT)
        ax.xaxis.label.set_fontsize(7.5); ax.yaxis.label.set_fontsize(7.5)
    for t in fig.findobj(Text):
        if t.get_text().strip():
            t.set_text(t.get_text().replace(r"$\mu$m", "\u00b5m"))   # upright micro sign
            t.set_fontfamily("Arial")
            t.set_fontsize(max(MIN_PT, t.get_fontsize()))


def _data_points(ax, r):
    """Display coordinates of every marker, and of plotted lines sampled densely."""
    pts = []
    inside = ax.get_window_extent(r)
    for ln in ax.get_lines():
        if not ln.get_visible() or ln.get_gid() == "guide":
            continue
        xy = ln.get_transform().transform(np.asarray(ln.get_xydata(), float))
        xy = xy[np.all(np.isfinite(xy), axis=1)]
        if len(xy) == 0:
            continue
        if ln.get_linestyle() not in ("None", "", " ", "none") and len(xy) > 1:
            seg = [np.linspace(a, b, 40) for a, b in zip(xy[:-1], xy[1:])]
            xy = np.vstack([xy] + seg)
        pts.append(xy)
    for c in ax.collections:
        if isinstance(c, PathCollection) and c.get_visible() and len(c.get_offsets()):
            pts.append(c.get_offset_transform().transform(c.get_offsets()))
    if not pts:
        return np.empty((0, 2))
    p = np.vstack(pts)
    m = (p[:, 0] >= inside.x0) & (p[:, 0] <= inside.x1) & (p[:, 1] >= inside.y0) & (p[:, 1] <= inside.y1)
    return p[m]


def audit(fig):
    """Text-text overlaps anywhere, and text/legend boxes covering plotted data (mm, at final size)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    mm = 25.4 / fig.dpi
    hidden = set()   # labels of ticks outside the view limits exist as Text objects but are not drawn
    for ax in fig.axes:
        for axis in (ax.xaxis, ax.yaxis):
            lo, hi = sorted(axis.get_view_interval())
            eps = 1e-9 * max(1.0, abs(hi - lo))
            for tk in axis.get_major_ticks() + axis.get_minor_ticks():
                if not lo - eps <= tk.get_loc() <= hi + eps:
                    hidden.update((id(tk.label1), id(tk.label2)))
    texts = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip() and id(t) not in hidden
             and t.get_window_extent(r).width > 0]
    boxes = [(t, t.get_window_extent(r).padded(-1.0)) for t in texts]
    found = []
    for (t1, b1), (t2, b2) in combinations(boxes, 2):
        if b1.overlaps(b2) and not (t1.get_gid() == "strip" and t2.get_gid() == "strip"):
            found.append(dict(kind="text_text", a=t1.get_text(), b=t2.get_text()))
    for ax in fig.axes:
        if getattr(ax, "_strip", False):
            continue
        pts = _data_points(ax, r)
        cands = [(t.get_text(), t.get_window_extent(r)) for t in ax.texts if t.get_visible() and t.get_text().strip()]
        leg = ax.get_legend()
        if leg is not None:
            cands.append(("<legend>", leg.get_window_extent(r)))
        for name, bb in cands:
            bb = bb.padded(-0.5)
            n = int(np.sum((pts[:, 0] > bb.x0) & (pts[:, 0] < bb.x1) & (pts[:, 1] > bb.y0) & (pts[:, 1] < bb.y1)))
            if n:
                found.append(dict(kind="text_over_data", a=name, points=n,
                                  at_mm=[round(bb.x0 * mm, 1), round(bb.y0 * mm, 1)]))
    sizes = sorted({round(t.get_fontsize(), 2) for t in texts})
    return found, sizes


def _export(fig, base, stem):
    fmts = ["pdf", "png"] + (["tiff"] if stem in TIFF_STEMS else [])
    rep = export_figure(fig, base, formats=fmts, dpi=DPI, tiff_rgb="tiff" in fmts, bbox_inches=None,
                        overwrite=True, mkdir=True, write_manifest=True,
                        metadata={"Title": stem, "Author": "S. A. Al Dajani"},
                        provenance={"profile": "scripta", "profile_accessed": PROFILE["accessed"],
                                    "source": "scripts/make_figures_kdense.py"})
    found, sizes = audit(fig)
    size = check_figure_size(fig, "elsevier")
    rec = dict(width_mm=round(size["figure"]["width_mm"], 1), height_mm=round(size["figure"]["height_mm"], 1),
               font_sizes_pt=sizes, min_font_pt=min(sizes), issues=found,
               outputs=[Path(o["path"]).name for o in rep["outputs"]])
    Path(f"{base}.audit.json").write_text(json.dumps(rec, indent=1, ensure_ascii=False) + "\n")
    AUDIT[Path(base).name] = rec


def save(fig, stem):
    w, h = fig.get_size_inches()
    fig.set_size_inches(WIDTH_IN, h * WIDTH_IN / w)
    normalize_text(fig)
    fig.tight_layout(w_pad=1.4, pad=0.3)
    _export(fig, OUT / stem, stem)
    if stem in STRIPS:
        add_strip(fig, *STRIPS[stem])
        _export(fig, OUT / "with_strips" / f"{stem}_strip", stem + "_strip")
    plt.close(fig)


# ----------------------------------------------------------------------------- panel (c) strip
CMAP = ListedColormap([FE, CR])


def _grid_ax(ax, grid):
    ax.imshow(grid, cmap=CMAP, vmin=0, vmax=1, interpolation="none", origin="lower")
    ax.set_xticks([]); ax.set_yticks([]); ax.minorticks_off()
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)
    ax._strip = True


def add_strip(fig, title, panels):
    """Append a row of atomic-plane snapshots as panel (c), keeping (a)/(b) at their size."""
    W, H = fig.get_size_inches()
    n = len(panels)
    gap, side = 0.05 * W / 6, None
    left = fig.axes[0].get_position().x0 * W
    right = max(a.get_position().x1 for a in fig.axes) * W
    cell = min(0.95, (right - left - (n - 1) * 0.12) / n)
    head, lab, foot = 0.20, 0.14, 0.04
    extra = head + lab + cell + foot
    pos = [(a, a.get_position()) for a in fig.axes]
    fig.set_size_inches(W, H + extra)
    Hn = H + extra
    for a, p in pos:
        a.set_position([p.x0, (p.y0 * H + extra) / Hn, p.width, p.height * H / Hn])
    yrow = foot / Hn
    total = n * cell + (n - 1) * 0.12
    x0 = (left + (right - left - total) / 2) / W
    for k, (label, draw) in enumerate(panels):
        ax = fig.add_axes([x0 + k * (cell + 0.12) / W, yrow, cell / W, cell / Hn])
        draw(ax)
        ax.set_title(label, fontsize=MIN_PT, pad=1.5).set_gid("strip")
    ytxt = (foot + cell + lab + 0.02) / Hn
    fig.text(left / W, ytxt, "(c)", fontsize=9, fontweight="bold", fontfamily="Arial", va="bottom", ha="left")
    fig.text((left + 0.32) / W, ytxt + 0.012 / Hn, title, fontsize=MIN_PT, fontfamily="Arial", va="bottom", ha="left")
    # colour key at the right end of the title row
    kx = right / W
    for name, col in (("Cr", CR), ("Fe", FE)):
        t = fig.text(kx, ytxt + 0.012 / Hn, name, fontsize=MIN_PT, fontfamily="Arial", va="bottom", ha="right")
        fig.canvas.draw()
        bb = t.get_window_extent(fig.canvas.get_renderer())
        sq = 0.075
        fig.patches.append(matplotlib.patches.Rectangle(
            (bb.x0 / fig.dpi / W - (sq + 0.03) / W, ytxt + 0.025 / Hn), sq / W, sq / Hn,
            transform=fig.transFigure, facecolor=col, edgecolor="none", figure=fig))
        kx = bb.x0 / fig.dpi / W - (sq + 0.14) / W


def random_grid(x, n=25, seed=7):
    return (np.random.default_rng(seed).random((n, n)) < x).astype(int)


def md_grid(path):
    t, p, box = fws.read_data(path)
    L = box[:, 1] - box[:, 0]
    n = round(L[0] / 2.866); a = L[0] / n
    rel = (p - box[:, 0]) / a
    on = np.abs(rel[:, 2] - np.round(rel[:, 2])) < 0.1
    i, j = np.round(rel[on, 0]).astype(int) % n, np.round(rel[on, 1]).astype(int) % n
    g = np.zeros((n, n), int); g[j, i] = (t[on] == 2).astype(int)
    return g


def strips():
    ser = series()
    snap = lambda f: (lambda ax: _grid_ax(ax, md_grid(f)))
    xs = [x for x in (0.1, 0.2, 0.36, 0.5, 0.7, 0.95) if x in ser]
    STRIPS["fig1_free_energy"] = (
        r"Fully decomposed MD structures ($\phi$ = 1), one (001) plane of a 25$^3$ bcc cell: "
        r"two phases form inside the spinodal window",
        [(f"{round(100 * x)}% Cr", snap(ser[x][-1])) for x in xs])
    rnd = [(f"{round(100 * x)}% Cr", (lambda x: lambda ax: _grid_ax(ax, random_grid(x)))(x))
           for x in (0.0, 0.2, 0.36, 0.5, 0.8, 1.0)]
    for stem in ("fig2_elastic_mapping", "fig2_elastic_mapping_orientation_violin_bulk"):
        STRIPS[stem] = ("Homogeneous (random) solid solutions, one (001) plane: each tensor in (a) and (b) "
                        "describes one of these", rnd)
    if 0.36 in ser:
        p36 = [(rf"$\phi$ = {float(f.stem[3:]):.1f}", snap(f)) for f in ser[0.36]]
        for stem in ("fig3_fsaw", "fig3_fsaw_relative", "fig3_fsaw_mddecomp", "fig3_fsaw_alpha",
                     "fig3_fsaw_crcolor", "fig3_fsaw_relative_crcolor"):
            STRIPS[stem] = (r"Fe$_{64}$Cr$_{36}$ decomposing (MD, one (001) plane): Fe-rich $\alpha$ and Cr-rich "
                            r"$\alpha'$ domains form as $\phi$ goes 0.5 $\to$ 1", p36)


def main():
    apply_publication_style("scripta", palette_name="okabe_ito_on_white")
    import matplotlib as mpl
    keep = dict(mpl.rcParams)
    mf.style, mf.save = style, save
    # make_figures.main() resets the font family; restore the preset right after it does
    orig_update = plt.rcParams.update
    plt.rcParams.update = lambda *a, **k: (orig_update(*a, **k), orig_update(
        {k2: keep[k2] for k2 in ("font.family", "font.sans-serif", "font.size", "axes.linewidth")}))
    strips()
    sys.argv = sys.argv[:1]
    mf.main()
    (OUT / "audit_summary.json").write_text(json.dumps(AUDIT, indent=1, ensure_ascii=False) + "\n")
    bad = {k: v["issues"] for k, v in AUDIT.items() if v["issues"]}
    print(f"\n{len(AUDIT)} figures -> {OUT}")
    print("min lettering (pt):", min(v["min_font_pt"] for v in AUDIT.values()))
    for k, v in bad.items():
        print(f"  {k}:")
        for i in v:
            print("    ", i)


if __name__ == "__main__":
    main()
