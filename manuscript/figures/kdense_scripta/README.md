# Scripta Materialia figure set (K-Dense)

Built by `scripts/make_figures_kdense.py` with the K-Dense helpers adapted for Scripta
(`scripts/kdense_scripta/`). Same data, series and panels as `manuscript/figures/`; presentation only:

- drawn at the final Elsevier full width, 190 mm, and exported without cropping (place at 100 %);
- Arial lettering at final size: 7 pt text, ticks and legends, 7.5 pt axis labels, 9 pt panel labels;
  mathtext in Arial and an upright µ;
- Okabe–Ito colour-blind-safe colours in the same hue families as the draft (blue, vermillion for
  red, amber for yellow, bluish green, reddish purple for orange);
- PDF (vector, TrueType fonts embedded) and 600 dpi PNG for every figure, plus 600 dpi RGB TIFF for
  the three manuscript figures (Elsevier: >= 500 dpi for combination artwork);
- `<stem>.audit.json`: size in mm, font sizes used, and an automatic check for text overlapping
  other text or plotted data (all figures currently pass); `<stem>.export.json`: export manifest;
- `with_strips/`: versions with panel (c), atomic-plane snapshots drawn natively (red Fe, blue Cr as
  in Acta Mater. 246 (2023) Fig. 7). Captions in `with_strips/CAPTIONS.md`.

Colour words to update in the captions if this set is used: DFT diamonds in Fig. 2a and literature
circles in Fig. 2b are amber (was yellow); red series are vermillion.
