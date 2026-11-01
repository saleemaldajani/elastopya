#!/usr/bin/env python3
"""All graphical abstract versions in one PDF to flip through (one version per page, with its notes from VERSIONS.md),
plus an overview page with all versions at once (2026-11-01).

Output: manuscript/figures/graphical_abstract_versions/all_versions.pdf
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

V = Path(__file__).resolve().parent.parent / "manuscript" / "figures" / "graphical_abstract_versions"


def main():
    notes = dict(re.findall(r"^\| (v\d+) \| [^|]* \| (.*?) \|$", (V / "VERSIONS.md").read_text(), re.M))
    vs = sorted((d for d in V.glob("v*") if d.is_dir()), key=lambda d: int(d.name[1:]))
    imgs = {d.name: mpimg.imread(d / "graphical_abstract.png")[::3, ::3] for d in vs}
    with PdfPages(V / "all_versions.pdf") as pdf:
        n = len(vs); cols = 3; rows = -(-n // cols)
        fig, axs = plt.subplots(rows, cols, figsize=(16, 2.6 * rows + 0.6))
        for ax in axs.flat:
            ax.axis("off")
        for ax, d in zip(axs.flat, vs):
            ax.imshow(imgs[d.name]); ax.set_title(d.name, fontsize=13, fontweight="bold", loc="left")
        fig.suptitle("Graphical abstract: all versions (newest last)", fontsize=15, fontweight="bold")
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        for d in vs:
            fig = plt.figure(figsize=(13, 6.2))
            ax = fig.add_axes([0.02, 0.12, 0.96, 0.80]); ax.imshow(imgs[d.name]); ax.axis("off")
            fig.text(0.02, 0.95, d.name, fontsize=18, fontweight="bold", va="center")
            fig.text(0.02, 0.06, notes.get(d.name, ""), fontsize=11, va="center", wrap=True)
            pdf.savefig(fig); plt.close(fig)
    print(V / "all_versions.pdf")


if __name__ == "__main__":
    main()
