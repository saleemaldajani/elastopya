# Fig. 4 — mechanism summary: notes (draft, 13 Oct 2026)

Goal: one schematic, full width (190 mm), explaining (i) the anomalous composition dependence of the
elastic/SAW response of bcc Fe–Cr and (ii) why spinodal decomposition stiffens the alloy, linking atoms,
moments, electrons and the TGS measurement. Every claim below is tagged with the evidence it rests on and
its status: **[shown]** (in our figures), **[running]** (Engaging jobs), **[lit]** (literature only), or
**[hypothesis]**.

## The causal chain in one line
Positive mixing enthalpy from magnetic frustration → spinodal α/α′ separation → each phase sits on a stiffer
part of a *convex* f(x) curve than the random alloy → stiffening, amplified where electronic and magnetic
anomalies make the random alloy unusually soft (near 36–40 % Cr), and enough to push the highest
(P)SAW modes toward the 609–833 MHz band.

## Proposed panels

### (a) Atomic scale: magnetic frustration and the driving force (ΔH, ΔG)
- Three bcc snippets with moment arrows, styled like Fig. 1(d):
  1. **Fe-rich (≤ 7–10 % Cr):** isolated Cr carries a moment antiparallel to Fe, about −1.8 μB (Klaver 2006) **[lit]**. This makes ΔH_mix slightly negative (Olsson 2003/2006, Levesque, Korzhavyi) **[lit]**.
  2. **Intermediate (20–50 % Cr):** Cr–Cr neighbours are frustrated: each Cr wants to be antiparallel to Fe and also to its Cr neighbour. Cr moments collapse.
     - QE SQS ⟨m_Cr⟩ = −0.58 / −0.21 / −0.03 μB at 20 / 35 / 50 % Cr, while ⟨m_Fe⟩ falls from 2.22 to 1.41 μB toward 95 % Cr **[shown]**, Fig. 1(d) lanes.
     - ΔH_mix > 0, peaking at about 0.08–0.09 eV/atom near 50–65 % Cr: QE 0.083, VASP 0.089, PAW literature 0.073–0.095 **[shown]**, Fig. 1a (entropy-sign-fixed version).
  3. **Decomposed:** a ferromagnetic Fe-rich α domain next to a Cr-rich α′ domain, where Cr recovers antiferromagnetic / SDW order **[hypothesis; SDW supercell running]**.
- Inset: ΔG_mix(x) and its curvature → spinodal window 0.21–0.75 (Fig. 1b) **[shown]**. Optionally ΔG(T) including vibrational entropy, about +0.2 k_B (Fultz 1995) **[running: phonopy SQS]**.

### (b) Electronic scale: Fermi surfaces and density waves
- **Fe-rich:** the majority-spin Fermi surface changes topology between 5 and 10 % Cr, an electronic topological transition (Korzhavyi 2009 Figs. 6–7; Olsson 2003 √-type B anomaly near 6 %) **[lit]**.
  - Our Fe Fermi surface in Korzhavyi's layout **[shown]**.
  - Unfolded SQS Fermi surfaces at 0 / 5 / 10 % Cr **[running]**.
  - Fe-rich elastic SM figure with the shaded 0.08–0.12 band **[shown]**.
- **Cr-rich:** the electron "jack" at Γ nests onto the hole octahedron at H, with Q = (2π/a)(1−δ), δ ≈ 0.05 → incommensurate SDW plus a CDW and strain wave at 2Q (Fawcett 1988) **[lit]**. Ours:
  - Cr Fermi-surface sections **[shown]**.
  - Nesting function ξ(q) along Γ–H, with rigid-band Fe filling to show the nesting being detuned **[running]**.
  - 40-atom SDW/CDW supercell: m(z), site-charge modulation at 2Q, forces along Q, and E(SDW) vs E(AFM) **[running]**.
  - The NM-Cr phonon instability at coarse k is the same nesting showing up as a Kohn anomaly (Fig. R8 / SM) **[shown]**.
- **Intermediate:** local Cr moments disappear near 40 % Cr, which is Razumovskiy's explanation of the non-monotonic C44, C′ and B **[lit]**. Our QE moments show the collapse **[shown]**. Whether the dip in f follows at the experimental lattice parameter is the aexp test **[running]**.
- Spin- and charge-density maps of the SQS on a (110) plane (where the Cr moment collapse is visible) **[running]**.

### (c) Elastic scale: why decomposition stiffens (the key quantitative panel)
- A schematic f(x) curve with the random alloy at x₀ = 0.36 and the decomposed alloy as a lever-rule chord between α (x ≈ 0.12) and α′ (x ≈ 0.74). The chord lies above the curve because f(x) is convex: decomposition stiffens.
- **Numbers (new, from data we have):** lever rule (1 − w′) f(x_α) + w′ f(x_α′), with the MD-measured phase compositions, against the actual stiffening on full decomposition:

| x₀ | Lever rule on MD f(x) | MD decomposition (φ = 1) | Lever rule on literature f(x) |
|---|---|---|---|
| 0.20 | +1.9 % | +2.6 % | −2.0 % |
| 0.36 | +4.0 % | **+5.1 %** | **+7.7 %** |
| 0.50 | +4.5 % | +5.5 % | +5.3 % |
| 0.70 | +3.3 % | +4.1 % | 0.0 % |

  - The composition (convexity) effect explains about 75–80 % of the MD stiffening at every x. The remaining roughly 1 % comes from what composition alone misses: α/α′ interfaces, coherency strain, and chemical short-range order within each phase **[shown; quantify in SM]**.
  - With the electronic (EMTO-CPA) f(x), the same lever rule gives **+7.7 % at Fe–36Cr**. The random alloy sits in the anomalous dip, so separation into stiffer α and α′ gains more than the classical potential predicts. This is the *quantum-elastodynamic* amplification. It is negative on the Fe-rich side, where the literature f(x) has its maximum near 20 % Cr.
  - With QE at the relaxed volume, f(x) is nearly linear, so the lever rule gives only a small gain. Whether the amplification survives depends on the aexp test **[running]**. Report both.
- Caveat: averaging frequencies linearly is a crude mixing rule. For the SM, redo it with Voigt/Reuss/Hill averages of the phase tensors and the anisotropic SAW.

### (d) Measurement scale: what TGS sees
- (P)SAW and bulk-wave spectra of the random alloy vs the decomposed alloy at λ = 6.4 μm. Decomposition shifts the distribution up, and the highest (P)SAW modes move toward the 609–833 MHz band of pristine CASS (Figs. 2b, 3) **[shown]**. Link to Rae et al. peak splitting **[lit]**.

## Electronic–mechanical interplay (the thread that ties (a)–(c) together)
An elastic constant is the strain curvature of the total energy. For a metal it splits into a "stiff" ionic and
electrostatic part and a band or Fermi-surface part, with C ≈ C_ion + C_band:
- **C_band ≈ −Σ_k (∂ε_k/∂η)² δ(ε_k − E_F) + ….** Strain η moves bands through E_F, electrons refill, and the
  lattice softens in proportion to the Fermi-level DOS and the deformation potentials. So whatever reshapes
  the Fermi surface reshapes the shear moduli:
  - **Electronic topological transition (Fe-rich, 5–10 % Cr):** a Fermi-surface neck opens or closes and gives a √-type anomaly in B
    and C′ (Lifshitz; Olsson 2003 near 6 %) **[lit]**. Our SM Fe-rich figure **[shown]**; unfolded Fermi surfaces **[running]**.
  - **Nesting (Cr-rich):** the same jack–octahedron nesting that drives the SDW produces Kohn anomalies. We see it directly
    as the coarse-k NM-Cr phonon instability and the Γ–H softening (Fig. R8) **[shown]**. The magnetoelastic anomalies of Cr
    at T_N (Fawcett 1988) **[lit]**. Fe detunes the nesting **[running: ξ(q) rigid band]**.
  - **Fermi-level DOS of the SQS (QE):** N(E_F) = 0.79 / 0.86 / 0.88 / 0.98 / 0.81 states/eV/atom at
    20 / 35 / 50 / 65 / 95 % Cr. Between 35 and 65 % Cr, C44 falls (130 → 116 GPa) while N(E_F) rises; C′ keeps rising
    (65 → 161 GPa), so the DOS is not the whole story **[shown, partial]**.
- **Magnetoelastic and magnetovolume coupling:** forming or losing local moments changes the volume and the
  bonding. Cr-moment collapse near 35–50 % Cr is Razumovskiy's origin of the C44 / C′ / B anomalies **[lit]**. Our moments
  collapse **[shown]**; the volume dependence is under test **[running: aexp]**.
- **Mechanical → electronic (the reverse arrow):** decomposition puts each phase at a composition, and the coherency strain puts it at a
  lattice parameter, where the Fermi surface is *not* at a topological or nesting singularity. The soft
  electronic term is switched off in both α and α′, which is why the electronic lever-rule gain (+7.7 %) exceeds the
  classical one (+4.0 %) at Fe–36Cr.
- **One-sentence version for the caption or abstract:** "Spinodal decomposition stiffens Fe–Cr because it moves
  matter out of a composition range where Fermi-surface and magnetic instabilities soften the shear moduli."

## Design direction ("make it hot")
- **Look:** one continuous left-to-right story band across the full width. A diverging Fe-red → Cr-blue colour field
  runs under all four panels, with glowing (soft-edged) Fermi-surface line art and ball-and-stick lattices with
  moment arrows (as Fig. 1d).
- **Hero element:** the convex f(x) curve with the lever-rule chord at Fe–36Cr and the gap labelled "+5 % (classical) →
  +8 % (electronic)", the money shot.
- **Typography:** big, bold panel verbs: FRUSTRATE → NEST → STIFFEN → DETECT.
- **Data:** real data in every inset, no cartoon numbers: moments, ΔG, ξ(q), f(x), spectra. Each inset regenerates from
  the scripts so it updates as the Engaging results land.

## Layout idea
Four columns left → right: atoms (a) → electrons (b) → elasticity (c) → measurement (d). Arrows labelled
"frustration", "nesting / ETT", "convexity (lever rule)", "SAW". Reuse the KDense style, Fe red / Cr blue,
moment arrows as in Fig. 1(d), and Fermi-surface line art from S_fermi_pure (Cr NM) with the Q arrow.

## Open questions before drawing
1. Does the aexp set reproduce the dip? This decides whether panel (b)'s "intermediate" story rests on our DFT or only on EMTO-CPA.
2. Does PBE stabilise the SDW relative to commensurate AFM in our supercell? Literature: nearly degenerate, sensitive to volume.
3. VASP ΔH at 95 % Cr (−0.042 eV/atom) is an outlier → keep it out of the mechanism panel.
4. Should panel (c) use the anisotropic SAW of the phase tensors rather than the isotropic f? The MD φ series gives this directly.
5. Main text vs SM: (a)–(d) as Fig. 4, with the density-wave figure and the Fultz phonon comparison in the SM.

## Sources to cite in the caption
Korzhavyi 2009; Olsson 2003, 2006; Klaver 2006; Levesque 2011; Razumovskiy 2011; Fawcett 1988; Fultz 1995;
Andersson & Sundman 1987; Dench 1963; Al Dajani et al. 2023 (Acta); Rae et al. (TGS).
