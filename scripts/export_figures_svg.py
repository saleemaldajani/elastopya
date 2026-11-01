#!/usr/bin/env python3
"""SVG copies of every figure in the paper, named by figure number (2026-11-01, user request).

For each \\includegraphics of the main text and the Supplementary Material (in order, numbered from the combined
paper's .aux), the SVG is
  - the native matplotlib SVG next to the source figure, when the figure scripts wrote one (editable text);
  - otherwise converted from the PDF with pdftocairo (vector; text as outlines);
  - for raster-only figures (PNG), an SVG wrapper that embeds the original PNG (still raster inside).
Panels of one figure get a, b, c. The other main-text figure styles (trend lanes, cubes) go to variants/.

  export_figures_svg.py <combined .tex> <combined .aux> <out dir>
Writes <out dir>/{main,supplementary,variants}/*.svg and <out dir>/index.csv (figure, source, how).
"""
import base64
import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

MS = Path(__file__).resolve().parent.parent / "manuscript"


def to_svg(src, dst):
    src = MS / src
    native = src.with_suffix(".svg")
    if native.exists() and native.stat().st_mtime >= src.stat().st_mtime - 120:   # native SVG from the same export
        shutil.copy2(native, dst)
        return "native matplotlib SVG"
    if src.suffix.lower() == ".pdf":
        subprocess.run(["pdftocairo", "-svg", str(src), str(dst)], check=True)
        return "converted from PDF (pdftocairo)"
    from PIL import Image
    w, h = Image.open(src).size
    b64 = base64.b64encode(src.read_bytes()).decode()
    mime = "image/png" if src.suffix.lower() == ".png" else "image/jpeg"
    dst.write_text(f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
                   f'width="{w}" height="{h}" viewBox="0 0 {w} {h}">'
                   f'<image width="{w}" height="{h}" xlink:href="data:{mime};base64,{b64}"/></svg>\n')
    return "raster (PNG embedded in SVG)"


def main():
    tex, aux, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    labels = dict(re.findall(r"\\newlabel\{([^}]*)\}\{\{([^}]*)\}", aux.read_text()))
    t = tex.read_text()
    i_sm = t.index(r"\section*{Supplementary Material}")
    rows = []
    for env in re.finditer(r"\\begin\{(figure\*?|sidewaysfigure)\}(.*?)\\end\{\1\}", t, re.S):
        body = env.group(2)
        lab = re.findall(r"\\label\{([^}]*)\}", body)
        num = labels.get(lab[-1], "") if lab else ""
        if not num:
            continue
        part = "main" if env.start() < i_sm else "supplementary"
        figs = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}", body)
        for k, f in enumerate(figs):
            name = f"Fig{num}" + (chr(ord("a") + k) if len(figs) > 1 else "")
            d = out / part; d.mkdir(parents=True, exist_ok=True)
            how = to_svg(f, d / f"{name}.svg")
            rows.append((part, name, f, how))
            if part == "main":                       # every other figure style of Figs. 1-3 (cubes, trend lanes)
                stem = Path(f); base = stem.stem.split("_strip")[0]
                for v in sorted((MS / stem.parent).glob(f"{base}_strip*.pdf")):
                    if v.name == stem.name:
                        continue
                    suf = v.stem[len(base) + len("_strip"):] or "_without_dft_row"
                    dv = out / "variants"; dv.mkdir(parents=True, exist_ok=True)
                    rel = v.relative_to(MS)
                    rows.append(("variants", f"{name}{suf}", str(rel), to_svg(rel, dv / f"{name}{suf}.svg")))
    with open(out / "index.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["part", "figure", "source", "how"]); w.writerows(rows)
    from collections import Counter
    print(len(rows), "SVGs;", dict(Counter(r[3] for r in rows)), dict(Counter(r[0] for r in rows)))


if __name__ == "__main__":
    main()
