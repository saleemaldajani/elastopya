# Homogenizing decomposed Fe–Cr for TGS (agent report, 2026-10-15; saved by the main session)

**Conclusion.** No homogenization of the QE or EMTO phase tensors puts decomposed Fe–36Cr in the band (>609 MHz, PSAW
intensity ≥ 2×10⁻³). The Voigt tensor is a rigorous upper bound in the positive-semidefinite order for any
microstructure, and the per-grain SAW is monotone in C, so the decomposed SAW cannot exceed the Voigt value (checked grain
by grain over 10⁴ grains: never more than 1 MHz above, solver resolution). The weak link is not the homogenization; it is
the phase tensors themselves (coherency strain / TOECs, α′ magnetism, true phase compositions) or the PSAW intensity
criterion (at 5×10⁻⁴, 1–3 % of decomposed grains are in band, 0.2–1.4 % of random).

## Validation against MD φ = 1 (12 compositions; RMS vs 5-seed mean)
| Scheme | f RMS (MHz) | shift recovered |
|---|---|---|
| random alloy | 17.1 | 0 % |
| Voigt | 3.27 | 80 % |
| Reuss | 5.07 | 73 % |
| Hill | 4.09 | 77 % |
| laminate, 3 ⟨100⟩ variants | 3.51 | 79 % |
| self-consistent ⟨100⟩ | 3.50 | 79 % |
| FFT on MD cells, two-phase | 3.72 | 79 % |
| FFT on MD cells, continuous composition (27-atom) | 2.84 | 84 % |
MD cannot rank schemes (both MD phases have A > 1; scheme differences ~1 MHz); the residual is the phase-composition representation.

## Fe–36Cr (α 0.14, α′ 0.88, w = 0.2973): highest per-grain mode / % > 609 MHz
| | QE (exp. a) | EMTO |
|---|---|---|
| random alloy | 600.3 / 0 | 551.3 / 0 |
| Voigt | 569.1 / 0 | 561.5 / 0 |
| Reuss | 569.1 / 0 | 581.8 / 0 |
| laminate, tetragonal | 569.7 / 0 | 571.4 / 0 |
| laminate, 3 variants | 568.5 / 0 | 570.3 / 0 |
| FFT ⟨100⟩ | 569.7 / 0 | 573.8 / 0 |
Phase tensors (GPa): QE α (237.7, 121.6, 106.7) A = 1.84; α′ (360.6, 93.5, 91.8) A = 0.69. EMTO α (263.7, 109.5, 111.0)
A = 1.44; α′ (385.4, 94.1, 101.9) A = 0.70. Medians: QE 507 → ~502 MHz; EMTO 489 → ~526 MHz.

## Methods (summary)
Voigt/Reuss/Hill; Backus laminate (layers ⊥ [001], welded), three-variant iso-strain cubic projection; Willis/HS and
self-consistent estimates; Moulinec–Suquet FFT on the MD cells and on synthetic ⟨100⟩ Cahn structures. General-anisotropy
per-grain SAW (saw_general.py) reproduces the cubic code bit-for-bit on 300 grains. Checks in test_homog.py all pass.
Recommendation: FFT ⟨100⟩ (or laminate V3 as analytic stand-in); quote Voigt as the rigorous upper bound.

Caveats: linear elasticity only (no 0.39 % coherency misfit / TOECs); QE α′ interpolated between 0.80 (AF-like) and 0.95;
EMTO values digitized; Hill and spherical SC not run through SAW (between Voigt and Reuss).
Files: fecr36_saw_summary.txt, md_validation_summary.txt, fecr36_tensors.csv; scripts homog.py, validate_md.py,
fecr36_tensors.py, saw_general.py.
