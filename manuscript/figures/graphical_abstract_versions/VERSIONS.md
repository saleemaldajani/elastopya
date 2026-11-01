# Graphical abstract versions

Each folder holds that version's outputs (PDF, PNG, TIFF at 600 dpi; `_screen.png` at journal display size) and the
script that made it (`make_graphical_abstract.py`; run it from `mike_final_touches/` to rebuild). To revert, copy a
version's script over `scripts/make_graphical_abstract.py` and run it. The current version is always
`manuscript/figures/graphical_abstract/`.

| Version | Date | What changed |
|---|---|---|
| v1 | 2026-11-01 | First draft: four numbered steps (unstable window, MD cubes, SAW stiffening, indicator box) |
| v2 | 2026-11-01 | DFT added: QE demixing energy (1) and QE C' (3); conclusion as banner |
| v3 | 2026-11-01 | MD and DFT side by side; DFT top (P)SAW entering the unallowed range (Fig. 2b) |
| v4 | 2026-11-01 | Two halves: key result plotted above its structures (MD cubes along phi; QE cells alpha, Fe-36Cr, alpha') |
| v5 | 2026-11-01 | v4 restyled: tinted cards, MD/DFT tags, open axes, highlighted alpha', banner |
| v6 | 2026-11-01 | Both cards on one x_Cr axis: all twelve MD compositions (phi = 0.7-1), six MD cells, five QE cells with C' |
| v7 | 2026-11-01 | "Spinodal decomposition" in title and banner; alpha / alpha' labels back |
| v8 | 2026-11-01 | Mechanism line as coloured chips including K/G down (embrittlement) |
| v9 | 2026-11-01 | Decluttered v8: phi = 1 curve only, QE cells without moment arrows, no window shading or extra alpha label in the DFT plot, mechanism line as plain coloured text |
| v10 | 2026-11-01 | One spinodal band: the 748 K CALPHAD spinodal (30-71% Cr, Sec. 3.1) instead of the two 298 K PWSCF-PAW windows of Fig. 1b, whose narrow second window near 95% Cr was confusing; +5.7% under the peak |
| v11 | 2026-11-01 | New design after the Acta 2023 graphical abstract: three panels with bold claims beneath (MD decomposed Fe-36Cr cell; (P)SAW spectrum over grain orientations with inset on the unallowed range, DFT alpha' reaching 0.69 GHz; f_SAW up and K/G down at all twelve MD compositions); no spinodal shading |
| v12 | 2026-11-01 | v11 plus the QE 40-atom cells of alpha (5% Cr) and alpha' (95% Cr) under the MD cube in panel 1, with arrows from the cube (MD separates, DFT models the two phases) |
| v13 | 2026-11-01 | v12 beautified, same content: rounded cards and claim strips, smoothed spectra, open axes, soft shadow under the cube, DFT cells symmetric about the cube, shorter legend |
| v14 | 2026-11-01 | Fe(1-x)Cr(x) made explicit (panel 1 claim 'of Fe1-xCrx (shown: x = 0.36)', panel 3 axis 'x in Fe1-xCrx'); Fe-36Cr marked '(TGS)', the alloy in which TGS detected spinodal decomposition (Rae et al.), not a CASS composition |
| v15 | 2026-11-01 | Centre legend ordered MD (green), DFT Fe-36Cr (blue), DFT alpha' (orange); panel 1 claim without Fe-36Cr ('Simulated spinodal decomposition of Fe1-xCrx'); cube smaller and clear of the key; Fe/Cr atom key under the DFT cells |
| v16 | 2026-11-01 | Headroom so the alpha' peak clears the centre legend; Fe-36Cr marker removed from panel 3 (Fe-36Cr named only in the centre legend) |
| v17 | 2026-11-01 | Wider frequency axis (0.27-0.80 GHz) so the centre legend ends before the band; smaller inset off the curves; panel-3 labels moved off their curves; panel-1 key inset from the card edge |
| v18 | 2026-11-01 | alpha / alpha' cell labels without percentages (the MD cube does not show a composition); SVG output added to the script |
