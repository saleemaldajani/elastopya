# Experiment and measured-tensor route (agent report, 2026-10-15; saved by the main session)

## Acta (al2023detecting) facts
- CF8 and CF8M, aged 290/330/360/400 °C for 10 kh and 30 kh plus unaged; 20 TGS spots each (10 air, 10 vacuum), Λ = 6.4 µm.
- Peaks by peak-finding on the filtered power spectrum; no amplitude criterion, no per-spot phase, no ferrite composition given.
- Table 3 mean second peak: only CF8 10 kh/360 °C (622.7 MHz) is above 609 MHz.
- Fig. 9 (digitized, ±1–2 counts): ~10–20 individual peaks per aged temperature in 627–833 MHz; unaged none above 627 (~3 in 597–631).
- Pristine band from an austenite (304) tensor 209/133/121 GPa, ρ 7950: 609.0 MHz reproduced only with a PSAW intensity cut ≥ 1–2×10⁻³ of the SAW (not stated in the paper); at ≤ 5×10⁻⁴ the edge is 665 MHz (cutoff-free ceiling). Fig. 9 grey band starts ~627 MHz.

## Measured data
- Single crystals only for Fe (231.4/134.6/116.4 GPa, A = 2.41; Ledbetter & Reed 1973, doi:10.1063/1.3253127) and Cr (350/67.8/100.8, A = 0.71; Bolef & de Klerk 1963 via Ledbetter & Kim 2001; not re-read).
- Fe–36Cr (Rae et al., arXiv:2601.06659): C′ = 56.7 GPa, C44 = 113 GPa, A = 1.99; C′ 31 % below Vegard. Aging 500 °C/10 h: C′ +5.4 %, soft-direction SAW +2.1 %, C44 unchanged.
- Duplex-steel bulk ultrasound: few-% changes of inconsistent sign (e.g. COBEM 2003 ≈ −3.5 %).

## SAW (10⁴ grains, cut 2×10⁻³)
Random Fe–36Cr (measured) 416–600 MHz, 0 % > 609; decomposed (Voigt/Hill/Reuss, anchored) 466–596, 0 %; APT-partial max 599; measured aged max 596.
The top of the band of an A > 1 tensor is shear-controlled: f_max ≈ √(C44/ρ)/Λ. Needed: +3 % in C44/ρ for 609 MHz, +9 % for 627, +90 % for 830.
At 5×10⁻⁴, weak high-velocity PSAWs reach ~840 MHz: 0.10 % of grains (random), 0.45 % (decomposed Voigt), 0.84 % (Vegard); pure Cr 28 %.

## Assessment
Effective tensor, anisotropy and α′ magnetism cannot lift the shear-controlled top; measured ~5 % stiffening (in C′, not C44) explains peaks to ~610–620 MHz, not to 830. The 650–830 MHz peaks need another mode: weak high-velocity (quasi-longitudinal-type) PSAWs, more frequent as decomposition lowers the anisotropy (or possibly a 2f₁ harmonic, unverified). Depends on the unreported intensity threshold and peak amplitudes.
Files: results.csv, tensors.csv, acta_figs/, runs/, digitize_fig9.py, effective_tensor.py.
