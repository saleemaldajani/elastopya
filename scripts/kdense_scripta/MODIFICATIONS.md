# K-Dense scientific-visualization helpers, adapted for Scripta Materialia

Vendored from K-Dense (MIT License, see `LICENSE.md`) and modified on 2026-10-12:

- `scripts/style_presets.py`: new `scripta` preset (7 pt Arial lettering incl. mathtext, closed
  frame with inward major+minor ticks, 0.6 pt axes, 600 dpi raster).
- `assets/publisher_profiles.json`: new `scripta` profile (Elsevier general widths 90/140/190 mm,
  500 dpi combination artwork, 7 pt / 6 pt lettering), marked as inherited because the journal
  Guide for Authors could not be fetched.
- `journal_requirements.md`: Scripta section.

`scripts/make_figures_kdense.py` (one level up) uses these to rebuild every figure of the paper at
the final 190 mm width into `manuscript/figures/kdense_scripta/`, checks lettering size, text–text
and text–data collisions, and writes an export manifest per figure.
