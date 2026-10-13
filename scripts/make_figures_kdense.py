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
TIFF_STEMS = {"fig1_free_energy", "fig2_elastic_mapping_orientation_violin_bulk_mddecomp",
              "fig3_fsaw_mddecomp_6p4um"}   # manuscript figures

# Okabe-Ito (Wong 2011) replacements, same hue family as the draft's Google-Sheets colours
mf.BLUE, mf.RED, mf.YELLOW, mf.GREEN, mf.ORANGE = "#0072B2", "#D55E00", "#E69F00", "#009E73", "#CC79A7"
mf.PURPLE = "#D55E00"   # Fe-36Cr in Fig. 3b (vermillion, distinct from Fe-20Cr blue)
mf.FIG = OUT
AUDIT = {}
STRIPS = {}   # stem -> (title, panels); filled in main()
TRENDS = {}   # stem -> list of trend lanes for the _strip_trend versions


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
            # legends need a visible gap (2 pt) to the data; free text only must not cover it
            bb = bb.padded(2 * fig.dpi / 72) if name == "<legend>" else bb.padded(-0.5)
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
        undo = add_strip(fig, *STRIPS[stem])
        _export(fig, OUT / "with_strips" / f"{stem}_strip", stem + "_strip")
        if TRENDS.get(stem):   # same strip with trend lanes under the blocks (separate files)
            undo()
            add_strip(fig, *STRIPS[stem], lanes=TRENDS[stem])
            _export(fig, OUT / "with_strips" / f"{stem}_strip_trend", stem + "_strip_trend")
    plt.close(fig)


# ----------------------------------------------------------------------------- panel (c) strip
CMAP = ListedColormap([FE, CR])


def _grid_ax(ax, grid):
    ax.imshow(grid, cmap=CMAP, vmin=0, vmax=1, interpolation="none", origin="lower")
    ax.set_xticks([]); ax.set_yticks([]); ax.minorticks_off()
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)
    ax._strip = True


LANE_H, LABEL_W, GAP = 0.46, 0.80, 0.12   # inches: trend-lane height, lane-label column, block spacing


def add_strip(fig, title, panels, lanes=()):
    """Append a row of atomic-plane snapshots as panel (c), keeping (a)/(b) at their size.

    lanes: optional trend rows drawn under the blocks, each
      dict(label, values (one per block), fmt, mode='seq' | 'pairs', zero=False)
    'seq' draws arrows between consecutive blocks; 'pairs' draws an arrow inside each (2i, 2i+1)
    pair with its % change. Returns a function that removes everything added."""
    W, H = fig.get_size_inches()
    n = len(panels)
    old_axes, old_texts, old_patches = list(fig.axes), list(fig.texts), list(fig.patches)
    pos = [(ax, ax.get_position()) for ax in old_axes]
    left = old_axes[0].get_position().x0 * W
    right = max(ax.get_position().x1 for ax in old_axes) * W
    lab_w = LABEL_W if lanes else 0.0
    cell = min(0.95, (right - left - lab_w - (n - 1) * GAP) / n)
    head, lab, foot = 0.20, 0.14, (0.12 if lanes else 0.04)
    lanes_h = len(lanes) * LANE_H + (0.04 if lanes else 0.0)
    extra = head + lab + cell + lanes_h + foot
    fig.set_size_inches(W, H + extra)
    Hn = H + extra
    for ax, q in pos:
        ax.set_position([q.x0, (q.y0 * H + extra) / Hn, q.width, q.height * H / Hn])
    total = n * cell + (n - 1) * GAP
    x0_in = left + lab_w + (right - left - lab_w - total) / 2
    yblock = foot + lanes_h
    centers = [k * (cell + GAP) + cell / 2 for k in range(n)]
    for k, (label, draw) in enumerate(panels):
        ax = fig.add_axes([(x0_in + k * (cell + GAP)) / W, yblock / Hn, cell / W, cell / Hn])
        draw(ax)
        ax.set_title(label, fontsize=MIN_PT, pad=1.5).set_gid("strip")
    for i, lane in enumerate(lanes):
        y0 = foot + (len(lanes) - 1 - i) * LANE_H
        ax = fig.add_axes([x0_in / W, y0 / Hn, total / W, (LANE_H - 0.04) / Hn])
        ax._strip = True
        _draw_lane(ax, lane, centers, total)
        fig.text(left / W, (y0 + (LANE_H - 0.04) / 2) / Hn, lane["label"], fontsize=MIN_PT, fontfamily="Arial",
                 va="center", ha="left", linespacing=1.1)
    ytxt = (yblock + cell + lab + 0.02) / Hn
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

    def undo():
        for ax in list(fig.axes):
            if ax not in old_axes:
                ax.remove()
        for t in list(fig.texts):
            if t not in old_texts:
                t.remove()
        fig.patches[:] = old_patches
        fig.set_size_inches(W, H)
        for ax, q in pos:
            ax.set_position(q)
    return undo


def _draw_lane(ax, lane, centers, total):
    """Values under each block, joined by trend arrows."""
    v = np.asarray(lane["values"], float)
    fmt, mode = lane["fmt"], lane.get("mode", "seq")
    lo, hi = np.nanmin(v), np.nanmax(v)
    if lane.get("zero"):
        lo, hi = min(lo, 0.0), max(hi, 0.0)
    pad = 0.45 * (hi - lo or abs(hi) or 1.0)
    ax.set_xlim(0, total); ax.set_ylim(lo - pad, hi + pad)
    ax.axis("off")
    if lane.get("zero") and lo < 0 < hi:
        ax.axhline(0, color="0.6", lw=0.5, ls=(0, (2, 2)), zorder=0)
        ax.text(total, 0, " 0", fontsize=MIN_PT, color="0.45", va="center", ha="left")
    ax.plot(centers, v, "o", color="0.15", ms=3.2, zorder=3, clip_on=False)
    arrow = dict(arrowstyle="-|>", color="0.25", lw=0.9, mutation_scale=7, shrinkA=3.5, shrinkB=3.5)
    pairs = [(k, k + 1) for k in range(len(v) - 1)] if mode == "seq" else [(k, k + 1) for k in range(0, len(v) - 1, 2)]
    for i, j in pairs:
        ax.annotate("", xy=(centers[j], v[j]), xytext=(centers[i], v[i]), arrowprops=arrow, zorder=2)
        if mode == "pairs":
            ax.text((centers[i] + centers[j]) / 2, max(v[i], v[j]) + 0.12 * pad, f"{100 * (v[j] / v[i] - 1):+.1f}%".replace("-", "\u2212"),
                    fontsize=MIN_PT, fontweight="bold", ha="center", va="bottom", color="0.1")
    for c, val in zip(centers, v):
        below = mode == "pairs" or lane.get("below", False)
        ax.annotate(fmt.format(val).replace("-", "\u2212"), (c, val), xytext=(0, -4.5 if below else 4.5), textcoords="offset points",
                    ha="center", va="top" if below else "bottom", fontsize=MIN_PT, color="0.15")


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
    # Fig. 2 with decomposed-MD tensors: the actual MD cells, random (phi = 0.5) next to decomposed (phi = 1)
    pairs = []
    for x in (0.2, 0.36, 0.5):
        if x in ser:
            pairs += [(f"{round(100 * x)}% Cr, random", snap(ser[x][0])),
                      (rf"{round(100 * x)}% Cr, $\phi$ = 1", snap(ser[x][-1]))]
    if pairs:
        STRIPS["fig2_elastic_mapping_orientation_violin_bulk_mddecomp"] = (
            "MD cells behind the green tensors in (b), one (001) plane: random alloy and fully decomposed "
            "structure at three compositions", pairs)
    if 0.36 in ser:
        p36 = [(rf"$\phi$ = {float(f.stem[3:]):.1f}", snap(f)) for f in ser[0.36]]
        for stem in ("fig3_fsaw", "fig3_fsaw_relative", "fig3_fsaw_mddecomp", "fig3_fsaw_alpha",
                     "fig3_fsaw_crcolor", "fig3_fsaw_relative_crcolor", "fig3_fsaw_mddecomp_6p4um"):
            STRIPS[stem] = (r"Fe$_{64}$Cr$_{36}$ decomposing (MD, one (001) plane): Fe-rich $\alpha$ and Cr-rich "
                            r"$\alpha'$ domains form as $\phi$ goes 0.5 $\to$ 1", p36)


def trends():
    """Values shown in the trend lanes of the _strip_trend versions (same data as the figures)."""
    lev = mf.read_csv("sheet_levesque_parameters.csv")
    a0, b0, c0, d0 = (lev[k][0] for k in ("alpha", "beta", "gamma", "delta"))     # PWSCF-PAW, as Fig. 1b
    dg_mix = lambda x: -1e3 * float(mf.levesque_model(np.array([x]), a0, b0, c0, d0)[0][0])   # meV/atom
    xs1 = [float(lab.split("%")[0]) / 100 for lab, _ in STRIPS["fig1_free_energy"][1]]
    TRENDS["fig1_free_energy"] = [dict(label="$\\Delta G_{\\rm mix}$\n(meV/atom)", values=[dg_mix(x) for x in xs1],
                                       fmt="{:+.0f}", zero=True)]
    mdr = mf.md_random()
    rnd_x = (0.0, 0.2, 0.36, 0.5, 0.8, 1.0)
    c11r = [float(np.interp(x, mdr["x"], mdr["c11"])) for x in rnd_x]
    for stem, lam in (("fig2_elastic_mapping", 5.5), ("fig2_elastic_mapping_orientation_violin_bulk", 6.4)):
        fr = [float(np.interp(x, mdr["x"], mdr["f"])) * 6.4 / lam for x in rnd_x]
        TRENDS[stem] = [dict(label="MD $C_{11}$\n(GPa)", values=c11r, fmt="{:.0f}"),
                        dict(label=f"MD $f_{{\\rm SAW}}$\n(GHz, {lam:g} µm)", values=fr, fmt="{:.3f}")]
    dd = mf.read_csv("md_decomposition.csv")
    dec = mf.md_decomposition()
    c11p, fp = [], []
    for x in (0.2, 0.36, 0.5):
        for ph in (0.5, 1.0):
            m = (np.abs(dd["x_cr"] - x) < 1e-9) & (np.abs(dd["phi"] - ph) < 1e-9)
            c11p.append(float(dd["c11"][m][0]))
            fp.append(float(np.interp(ph, dec[x]["phi"], dec[x]["f"])))
    TRENDS["fig2_elastic_mapping_orientation_violin_bulk_mddecomp"] = [
        dict(label="MD $C_{11}$\n(GPa)", values=c11p, fmt="{:.0f}", mode="pairs"),
        dict(label="MD $f_{\\rm SAW}$\n(GHz, 6.4 µm)", values=fp, fmt="{:.3f}", mode="pairs")]
    ph36, f36 = dec[0.36]["phi"], dec[0.36]["f"] * 6.4 / 5.5
    rel36 = 100 * (f36 / f36[np.argmin(ph36)] - 1)
    for stem in ("fig3_fsaw", "fig3_fsaw_mddecomp", "fig3_fsaw_alpha", "fig3_fsaw_crcolor"):
        TRENDS[stem] = [dict(label="$f_{\\rm SAW}$\n(GHz, 5.5 µm)", values=list(f36), fmt="{:.3f}")]
    TRENDS["fig3_fsaw_mddecomp_6p4um"] = [dict(label="$f_{\\rm SAW}$\n(GHz, 6.4 µm)", values=list(dec[0.36]["f"]),
                                               fmt="{:.3f}")]
    for stem in ("fig3_fsaw_relative", "fig3_fsaw_relative_crcolor"):
        TRENDS[stem] = [dict(label="$\\Delta f_{\\rm SAW}$\n(%)", values=list(rel36), fmt="{:+.1f}")]


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
    trends()
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
