#!/usr/bin/env python3
"""Shrink a paper PDF for e-mail (2026-10-17): raster images wider or taller than MAXPX are downsampled to MAXPX and
re-encoded as JPEG (quality Q); vector graphics and text are untouched; identical objects are merged. The full-resolution
PDFs stay in to_send/full_resolution/. Needs pypdf + Pillow (the scratch venv used for annotation extraction has both).

  compress_pdf.py in.pdf out.pdf [maxpx] [quality]
"""
import sys

from pypdf import PdfReader, PdfWriter


def main():
    src, dst = sys.argv[1], sys.argv[2]
    maxpx = int(sys.argv[3]) if len(sys.argv) > 3 else 1800
    q = int(sys.argv[4]) if len(sys.argv) > 4 else 88
    w = PdfWriter(clone_from=PdfReader(src))
    n = 0
    for page in w.pages:
        for im in page.images:
            pil = im.image
            if max(pil.size) <= maxpx and len(im.data) < 150_000:
                continue
            s = min(1.0, maxpx / max(pil.size))
            if s < 1:
                pil = pil.resize((round(pil.size[0] * s), round(pil.size[1] * s)))
            if pil.mode not in ("RGB", "L"):
                pil = pil.convert("RGB")
            im.replace(pil, quality=q)
            n += 1
    w.compress_identical_objects(remove_identicals=True, remove_orphans=True)
    with open(dst, "wb") as fh:
        w.write(fh)
    print(f"{src} -> {dst}: {n} images re-encoded")


if __name__ == "__main__":
    main()
