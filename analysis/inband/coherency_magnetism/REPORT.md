# Coherency stress and α′ magnetism (agent report, 2026-10-15; saved by the main session)

**Answer:** at the paper's PSAW threshold (2×10⁻³) neither ingredient moves decomposed Fe–36Cr into the band; all variants
top out at 566.7–578.1 MHz (edge 609). The top mode is a shear-type PSAW, f ∝ C44^0.51, nearly independent of C11, C12;
609 MHz needs ~15 % higher effective C44 (~+50 GPa in α′).

## A. Coherency stresses (TOECs)
Measured Fe TOECs (Powell & Skove 1968; Rotter & Smith 1966; Guinan & Beshers 1968 via Ledbetter & Reed 1973); no
measured Cr TOECs, so Cr from pressure data (Rodière et al. 2026) and QE volume derivatives. Each phase changes by up to
~10 GPa; stresses balance, so the Fe–36Cr effective tensor changes by < 0.9 GPa. SAW on the 40 top grains: −0.8 to
+0.7 MHz; 10× the real misfit reaches only 594 MHz.

## B. α′ magnetic state
At room temperature Cr with 12 % Fe is paramagnetic (T_N 190 K at 10 % Fe, gone by 15 %; Benka et al.). Paramagnetic ≠
non-magnetic: measured RT Cr (350/67.8/100.8 GPa) has B = 162 GPa, close to our magnetic DFT states (153–186) and far from
non-magnetic (218–264); the softening persists above T_N (Katahara et al. 1979). Non-magnetic α′ is the least realistic.
α′ alone ranges 549–684 MHz across choices, but Fe–36Cr moves only −2.4 to +9.0 MHz (best 578 MHz, measured tensors).

## PSAW threshold
| cut | Fe–36Cr grains > 609 MHz |
|---|---|
| 2×10⁻³, 1×10⁻³ | none |
| 1×10⁻⁴, decomposed | 1.2–7.0 % (weak modes ~750–810 MHz) |
| 1×10⁻⁴, pristine SS304 | 0.01 % (one grain, 665 MHz) |
At 1×10⁻⁴ decomposition raises the fraction: QE 0.37 → 2.1–3.0 %, EMTO 2.3 → 7.0 %, MD 0.02 → 0.09 %.
Caveats: "measured" α′ is pure Cr; several sources read second-hand or abstract only.
Files: tensors_magnetic.csv, thresholds.csv, coherency.py, coherency_saw.py, magnetic_tensors.py, saw/.
