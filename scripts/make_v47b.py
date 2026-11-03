#!/usr/bin/env python3
"""Paper v47b (2026-11-03): v47 (the submitted v44 without section headings) with one leftover fixed: the Methods
sentence ended "(see Data Availability)", a heading that no longer exists; the phrase is dropped (the repository URL
stays in the closing paragraph). Nothing else changes. v47 is kept.

  make_v47b.py    writes finalpaper_qe_mechanism_v47b_validation[_cubes|_trend|_trend_cubes].tex from v47
"""
from pathlib import Path

MS = Path(__file__).resolve().parent.parent / "manuscript"
OLD = "provided in the accompanying repository (see Data Availability)."
NEW = "provided in the accompanying repository."

for v in ("", "_cubes", "_trend", "_trend_cubes"):
    src = MS / f"finalpaper_qe_mechanism_v47_validation{v}.tex"
    dst = MS / f"finalpaper_qe_mechanism_v47b_validation{v}.tex"
    assert not dst.exists(), dst
    t = src.read_text()
    assert t.count(OLD) == 1, (v, t.count(OLD))
    t = t.replace(OLD, NEW)
    assert "Data Availability" not in t[:t.index(r"\section*{Supplementary Material}")]
    dst.write_text(t)
    print("wrote", dst.name)
