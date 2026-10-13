# Fawcett (1988): SDW antiferromagnetism in Cr. Facts extracted for the DFT comparison

Source: E. Fawcett, "Spin-density-wave antiferromagnetism in chromium", Rev. Mod. Phys. 60, 209–283 (1988).
Page numbers below are journal pages (PDF page = journal page − 208). The text was read from the scan's OCR layer, and every number quoted here was checked against the page image.
Values marked "≈ (read from fig.)" were read by eye from a figure. Allow about ±3–5 % for these, or ± the last digit given.

Companion CSVs (digitized):
- `cr_sdw_delta_vs_T_fawcett1988.csv`: Q(T) and δ(T) for pure Cr, from Fig. 4(b).
- `cr_fe_neel_fawcett1988.csv`: T_N vs at.% Fe (Fig. 5) and Q at 100 K vs at.% Fe (Fig. 6).

---

## (a) SDW wave vector, δ(T), Néel and spin-flip temperatures

- **Wave vector.** Magnetic satellites appear at (0,0,1±δ), with **Q = (2π/a)(1 ± δ)** along a cube axis (p. 211; Bykov et al. 1959, Corliss et al. 1959). In the review's convention Q = (0,0,1−δ)a*, with a* = 2π/a and δ = 1 − Q (Fig. 3 caption, p. 213; Fig. 4 caption, p. 214).
- **δ at low T.** **δ = 0.048** at low temperature and **δ = 0.037 near T_N** (p. 246, text; also Fig. 4(b), p. 214). From Fig. 4(b), Q ≈ 0.9516 (δ ≈ 0.0484) below ~80 K and Q ≈ 0.963 (δ ≈ 0.037) at ~310 K. Other quoted values: δ = 0.047 at 130 K (p. 220), and δ ≈ 0.05 (p. 214).
- **T dependence.** "Both Q and δ appear to become constant and independent of temperature below about 80 K". However, synchrotron measurements of the 2Q CDW by Gibbs et al. (1987) show dQ/dT approaching zero only below 40 K (p. 214). The modulation wavelength changes continuously from ≈21 unit cells at 78 K to ≈28 unit cells near T_N (Shirane & Takei 1962; p. 211).
- **Néel temperature.** **T_N ≈ 311 K** (pp. 210–211). Fig. 80 gives T_N = 311.4 K. The transition is **weakly first order** (pp. 212, 263–264). It is superimposed on a continuous transition about 3 K wide, with latent heat ≈ 1.2 ± 0.2 J mol⁻¹ (p. 264) and total magnetic entropy ≈ 35 mJ mol⁻¹ K⁻¹ (p. 264, Table XVI).
- **Spin-flip transition.** **T_SF ≈ 123 K** (p. 211). Above T_SF the SDW is **transverse** (S ⊥ Q, "TSDW", 123–311 K). Below it the SDW is **longitudinal** (S ∥ Q, "LSDW") (pp. 211, 214; Fig. 4(a)). Hastings (1960) first located the transition near 110 K (p. 211).
- **Pressure dependence.** dT_N/dp = −0.51 K GPa⁻¹ (Mitsui & Tomizuka 1965, p. 266). d ln T_N / d ln V = +26.5, taking B = 162 GPa (McWhan & Rice 1967, p. 267).

## (b) SDW amplitude and moment at low T

- At 4.2 K the rms moment is **⟨μ(0)⟩ = 0.43 μ_B**, and the **peak SDW amplitude is μ₀ = √2⟨μ(0)⟩ = 0.62 μ_B** per Cr atom (Arrott et al. 1967, p. 217).
- Shirane & Takei (1962) gave μ = 0.59 μ_B at 78 K, which is consistent with the 4.2 K value (p. 217). In the 1950s, powder Cr gave 0.40 μ_B at 4.2 K for the (strain-induced) commensurate (100) peak (p. 211).
- In commensurate AFM **CrMn** the moment rises to about **0.8 μ_B per atom** (Koehler et al. 1966; p. 216).
- Theory: Asano & Yamashita obtained μ₀ = 0.6 μ_B with Xα, α = 0.5. Full Slater exchange (α = 1) gives μ₀ > 3 μ_B, which is far too large (p. 248).
- Magnetovolume: ω_M(T → 0) ≈ 0.15 % (read from Fig. 40, p. 236). It scales with ⟨μ²⟩, with C/B = 2.4 % μ_B⁻² (Fawcett, Kaiser & White 1986; p. 236, Eq. 37).

## (c) Charge-density wave and strain wave (second harmonic)

- The incommensurate SDW of wave vector Q is accompanied by a **second-harmonic CDW and strain wave (periodic lattice distortion) with wave vector 2Q** (period λ_SDW/2). They appear as satellites of nuclear Bragg peaks at (h,k,l ± 2δ)a* (pp. 228–229, Eq. 20a, Fig. 27).
  - The strain wave was first seen by x rays (Tsunoda et al. 1974). The CDW was predicted via electron–phonon coupling (Young & Sokoloff 1974).
  - The displacement A₂ is **longitudinal (∥ Q)** (Eagen & Werner 1975; p. 229).
  - A₂ ∝ S₁² holds over the whole T range, with no change at T_SF (Fig. 28, p. 229).
  - A third-harmonic SDW at 3Q (S₃ ∝ S₁³) is also observed.
- **Strain-wave amplitude A₂/a.** The recommended value is **A₂/a = (1.45 ± 0.4) × 10⁻³ at 153 K** (from Pynn et al. 1976, rescaled; p. 232). Other determinations:

  | Value of A₂/a | Temperature | Method | Page |
  |---|---|---|---|
  | (1.3 ± 0.3) × 10⁻³ | 200 K | neutrons | p. 232 |
  | (1.7 ± 0.2) × 10⁻³ | 130 K | Eagen & Werner, neutrons | p. 231 |
  | (2.5 ± 0.5) × 10⁻³ | 153 K | x rays, Tsunoda et al. 1974 | p. 232 |
  | 3.5 × 10⁻³ | not stated | indirect, Tsunoda et al. 1975 | p. 231 |

- **Phase of the strain wave.** A₂/S₁ is negative, as shown in Fig. 32(b) on p. 231. This is the converse of case (a), in which the minimum lattice spacing coincides with the SDW antinode.
- **CDW amplitude.** Mori et al. (1975) measured the relative charge-density modulation ρ = ρ⁰[1 + σ cos(2δ·r)] by x rays at 150 K. They found **σ ≈ 8 × 10⁻³**, with a positive sign: "the electron density is large at the antinode of the SDW". The CDW amplitude is of the same order as the strain wave (p. 232). Theoretical estimates are σ ≳ 10⁻² (Young & Sokoloff), σ ≈ (5 ± 2) × 10⁻³ (Kotani) and σ ≈ 1.4 × 10⁻³ (Nakajima & Kurihara) (p. 232).
- **Third harmonic.** |S₃/S₁| = (1.65 ± 0.05) × 10⁻² at 200 K (Pynn et al. 1976). S₃/S₁ is **negative**, which gives a "rectangular" (squared-up) waveform (p. 232–233; Figs. 31, 33).
- Tetragonality of single-Q Cr: just below T_N c/a > 1, and the sign changes to c/a < 1 at 230 K (Steinitz et al. 1969; p. 212).

## (d) Fermi-surface nesting

- The SDW originates from nesting between the **electron "jack" centred at Γ** and the **hole octahedron centred at H**. The bodies of both are approximately octahedral with nearly flat faces. They are separated by a nesting vector Q̄ that is close to the SDW wave vector Q (Lomer 1962; pp. 211, 246–247).
- Other Fermi-surface pieces are the **hole pockets (ellipsoids) at N** and the **electron "balls"/lens along ΓH**. The lens touches the jack (pp. 246–247; Figs. 57–59).
- **Nesting vector** Q̄ = (0,0,1−δ̃)a*:
  - Loucks (APW) gave Q̄ = 0.90a*, against the experimental Q = 0.95a* (p. 246).
  - In Rath & Callaway's (1973) calculation, δ̃ varies between 0.024 and 0.055 over the (001) plane through Γ, and between 0.026 and 0.037 in a plane displaced by a*/24 (Fig. 54, pp. 246–247). Experiment gives δ = 0.037 near T_N and 0.048 at low T.
  - Rigid-band estimate: Q̄ should change by about 1 % per at.% of V or Mn. Experiment shows about twice this rate (p. 246).
- **Kohn-anomaly calipers** (Shaw & Muhlestein 1971 phonons, room T; Table VII, p. 247):

  | Nesting feature | Experiment | KSG theory | VBH theory |
  |---|---|---|---|
  | [111] electron-jack body to hole octahedron at H | 0.94 ΓH | 0.98 | 0.99 |
  | [011] electron jack to hole pocket at N | 0.90 ΓN | 0.84 | 0.84 |
  | [111] jack ball to hole N | 0.54 ΓH | 0.54 | 0.54 |
  | [111] extremum of ball on jack | 0.50 ΓP | 0.54 | 0.52 |
  | [001] anomaly | 0.85 ΓH | 0.82 | 0.83 |

  The [001] anomaly at 0.85 ΓH was attributed by Shaw & Muhlestein to jack–octahedron nesting. Fawcett notes that 0.85a* is smaller than Q = 0.95a* (p. 247).
- Møller & Mackintosh (1965) found phonon energies unchanged between 100 K (AFM) and 400 K (paramagnetic). Only the linewidths grew with T (p. 233). This supports using the 300 K phonon data for paramagnetic-Cr comparisons.

## (e) Fe in Cr (CrFe alloys)

The review gives Fe data **only in Figs. 5 and 6** (p. 215; CrFe data from Burke & Rainford 1978). Digitized values are in `cr_fe_neel_fawcett1988.csv`.

- **T_N decreases with Fe content.** Starting from 311 K for pure Cr, Fig. 5 gives approximately:

  | at.% Fe | T_N (K) |
  |---|---|
  | 1.5 | 271 |
  | 1.8 | 268 |
  | 2.3 | 258 |
  | 3.3 | 251 |
  | 4.9 | 232 |

  This is roughly −15 to −20 K per at.% Fe. Uncertainty is ±10 K: Fig. 5's pure-Cr square reads ≈319 K, against 311 K in the text.
- **Commensurate/incommensurate.** Fe raises Q, so δ decreases. At 100 K (Fig. 6), δ goes from 0.048 (Cr) to ≈0.045 (0.5 % Fe), ≈0.041 (0.75–0.9 % Fe) and ≈0.036–0.040 (1.7 % Fe). At **≈2.3 at.% Fe the SDW is commensurate (Q = 1.00)**. Fig. 5 shows, schematically as a dotted line, a transition from the incommensurate SDW to a **commensurate SDW (CSDW) in CrFe for x ≳ 2.3 at.%** below ≈255 K.
- **Moment vs Fe content.** No values are given for CrFe. The only alloy moment quoted is about 0.8 μ_B in commensurate CrMn (p. 216).
- **Local moments on Fe.** The review says **nothing** about local moments on Fe impurities in Cr, and nothing about spin-glass behaviour of Cr–Fe. Its "local moment" discussions concern only the absence of local moments in paramagnetic pure Cr (pp. 211, 217, 228), and a band calculation in which Fe, Co and Ni (but not Cr) develop disordered local moments (p. 249).
- Impurity effects on the transition: 0.5 at.% V turns the first-order Néel transition into an apparently continuous one, and Mo does the same (p. 267). Fe is not discussed there.

## Elastic-constant (magnetoelastic) anomalies at T_N in Cr

Sec. III.D (pp. 236–237) and Sec. VI.A.3 (pp. 265–266) cover these anomalies.

- **Young's modulus.** Fine et al. (1951) found a ΔE softening in polycrystalline Cr. It begins well above T_N and has a sharp minimum at T_N; Young's-modulus data also showed the spin-flip transition near 121 K (pp. 210, 236).
  - Munday & Street (1971) measured E along [110] (Fig. 41, p. 236). Reading the figure:
    - ≈2.83 × 10¹¹ Pa in the LSDW phase at 100 K.
    - A drop at T_SF to 2.72–2.78 × 10¹¹ Pa.
    - Broad minima of ≈2.57 × 10¹¹ Pa (single-Q) to ≈2.64 × 10¹¹ Pa (field-cooled along [1̄10]) near 250 K.
    - A jump at T_N to ≈2.70–2.72 × 10¹¹ Pa in the paramagnetic phase.
  - So the TSDW phase is up to about 5 % softer than the paramagnetic value.
  - In the single-S single-Q state the intrinsic change is "relatively small" (dotted curve). The large effects come from rotation of the polarization S, i.e. technical antiferromagnetism.
- **Longitudinal moduli C₁₁, C₂₂, C₃₃** (Muir et al. 1987b; Fig. 43, p. 237). These show the strong anomaly. Reading the figure:
  - ≈3.87–3.92 × 10¹¹ Pa at 60–70 K.
  - Decreasing to ≈3.5–3.6 × 10¹¹ Pa at about 290 K.
  - A sharp cusp at T_N reaching ≈3.27 × 10¹¹ Pa (C₂₂, C₃₃) and ≈3.48 × 10¹¹ Pa (C₁₁, read on its offset right-hand axis).
  - Recovery to ≈3.5–3.55 × 10¹¹ Pa just above T_N.
- **Shear and quasishear moduli.** Muir et al. find ΔC_ii ≫ ΔC_Mi ≈ ΔC_αα (Eq. 81, p. 265). That is, **C₄₄ and C′ = ½(C₁₁ − C₁₂) (≈ C_M, Eq. 84) have only small anomalies** at T_N. C₄₄ is ≈1.00–1.035 × 10¹¹ Pa over 50–350 K, with changes of about 1–2 % (Fig. 44, p. 237).
  - The exception is **C₆₆**, which has "a fairly strong magnetic anomaly below T_N" (p. 265).
  - This happens because T_N depends mainly on **volume strain** and only weakly on shear or tetragonal strain (Eqs. 87–88, p. 265).
- **Bulk modulus.** The magnetic anomaly ΔB(t) (Cr vs paramagnetic Cr₉₅V₅) is proportional to the anomaly in volume expansivity Δβ(t), both below and above T_N (Fig. 82, p. 266).
  - ΔB reaches ≈4 × 10¹⁰ N m⁻² near T_N, at 304 K (read from Fig. 82). The reference value is B̄(Cr₉₅V₅) = 2.07 × 10¹¹ N m⁻². That makes the softening about 20 % at T_N (my ratio, not stated in the paper).
  - ΔB is still ≈2 × 10¹⁰ N m⁻² at 130 K and at 380 K, and ≈0.3 × 10¹⁰ N m⁻² at 500 K (read from Fig. 82).
  - The resulting Grüneisen parameters are γ_{T_N−} = −37 (T < T_N) and γ_{T_N+} = −144 (T > T_N; paramagnetic spin fluctuations) (p. 266).
- **Under pressure** at T < T_N, Katahara et al. (1979) find a strong anomaly in C₁₁ at the Néel pressure p_N. C₄₄ and C_M show no discontinuity, only a change in dC/dp (p. 266).
- **Volume/strain at T_N.** The first-order strain discontinuity Δε ≤ 11 × 10⁻⁶ (an upper limit) and is isotropic to about 5 % (p. 264, Fig. 81). The anomaly in linear expansivity is Δα ≈ 4.0 × 10⁻⁶ K⁻¹ (Lee & Asgar 1969; p. 265).
- **Ultrasonic attenuation** at T_N in pure Cr is very large (> 40 dB cm⁻¹). It is strongly suppressed by 0.67 at.% V (p. 267).
- Fultz et al. (1995) cite this review for the statement that in AFM Cr and Cr-rich Cr–Fe the elastic constants are strongly affected by the Néel transition. This review itself contains no elastic-constant data for Cr–Fe alloys.
