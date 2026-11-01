#!/usr/bin/env python3
"""Paper v44 (2026-11-01): v43 with AI-style phrasing removed from the main text (user request): joining semicolons split,
empty "rather than"/"not X" contrasts and summary sentences without content deleted, Discussion paragraphs 1-3 cut to
their content, generic Conclusions replaced by the paper's findings. Deletion first; the user's wording kept elsewhere.
Only the main text changes (the SM, incl. the full text of Secs. 2-3 in S11, is untouched); the Highlights block is
taken out of the manuscript because Elsevier takes highlights as a separate file. v43 is kept.

  make_v44.py     writes finalpaper_qe_mechanism_v44_validation[_cubes|_trend|_trend_cubes].tex from v43
"""
from pathlib import Path

MS = Path(__file__).resolve().parent.parent / "manuscript"

EDITS = [
    # Highlights are a separate submission file (upload/2_highlights), not part of the manuscript
    (r"""\textbf{Highlights}

\begin{itemize}
\item Ab initio simulations predict stiffening during spinodal decomposition in $\delta$--ferrite
\item Composition-dependent elasticity maps to surface acoustic wave frequency shifts
\item Theory and experiment agree: stiffening non-destructively indicates embrittlement
\end{itemize}
""", ""),
    # Sec. 2.3: semicolon, "not VASP"
    (r"""These simulations were performed with LAMMPS \cite{plimpton1995}, not VASP, and include no electronic structure.""",
     r"""These classical simulations were performed with LAMMPS \cite{plimpton1995} and include no electronic structure."""),
    (r"""phase-separated by Monte Carlo atom swaps; intermediate states were produced""",
     r"""phase-separated by Monte Carlo atom swaps. Intermediate states were produced"""),
    # Sec. 2.4: semicolon
    (r"""the isotropic VRH moduli used in Eq.~\ref{eqn:rayleigh} are calculated as follows; the same moduli enter the Pugh embrittlement criterion and the thermal estimates reported in the Supplementary Material.""",
     r"""the isotropic VRH moduli used in Eq.~\ref{eqn:rayleigh}, in the Pugh embrittlement criterion and in the thermal estimates of the Supplementary Material are calculated as follows:"""),
    # Sec. 3.1: summary sentence without content
    (r""" Instead, any observable elastic response must be understood as emerging from a system that is intrinsically driven toward compositional modulation, providing the foundation for linking free-energy curvature to stiffness evolution and ultimately to ultrasonic observables.""",
     ""),
    # Sec. 3.2
    (r"""so we interpret relative trends rather than absolute frequencies; Quantum ESPRESSO calculations""",
     r"""so we interpret relative trends. Quantum ESPRESSO calculations"""),
    (r"""in which no surface mode of the pristine material exists; (P)SAW peaks observed there""",
     r"""in which no surface mode of the pristine material exists. (P)SAW peaks observed there"""),
    (r"""in step with the SAW shift; meanwhile the shear-modulus misfit between $\alpha$ and $\alpha'$ (+16 to +23~GPa from QE and EMTO-CPA) hardens the ferrite""",
     r"""in step with the SAW shift. The shear-modulus misfit between $\alpha$ and $\alpha'$ (+16 to +23~GPa from QE and EMTO-CPA) also hardens the ferrite"""),
    (r""" Within the spinodally unstable composition window, elastic stiffening therefore translates into measurable acoustic signatures, enabling stiffness changes to be inferred non-destructively.""",
     ""),
    # Sec. 3.3
    (r"""We therefore report two complementary relationships that isolate these contributions within a unified elastodynamic framework.""",
     r"""We therefore treat the two separately."""),
    (r"""compositions (arrows); in both directions the homogeneous-alloy frequency increases.""",
     r"""compositions (arrows), and in both directions the homogeneous-alloy frequency increases."""),
    (r"""nearly vanishes for the dilute alloys at either end; the fully separated state is an idealized limit""",
     r"""nearly vanishes for the dilute alloys at either end. The fully separated state is an idealized limit"""),
    (r""" Together, $f_{\mathrm{SAW}}(x_{\mathrm{Cr}})$ captures the baseline elastic response of a homogeneous alloy, and $f_{\mathrm{SAW}}(\phi)$ reflects the additional stiffening associated with progressive spinodal conversion.""",
     ""),
    # Discussion: paragraphs 1-3 cut to their content
    (r"""Taken together, the results establish a closed thermodynamics-to-observable pathway for interpreting picosecond ultrasonic measurements in $\delta$--ferrite--containing stainless steels. The free-energy curvature criterion identifies the composition window in which the homogeneous solid solution is intrinsically unstable, defining the regime where elastic and acoustic properties must reflect ongoing or incipient spinodal decomposition rather than equilibrium alloy behavior.

Within this thermodynamically unstable regime, first-principles elasticity reveals systematic stiffening trends that are propagated through elastodynamic constraints into Rayleigh SAW frequency shifts. The direct mapping between elastic constants and SAW observables provides a quantitative basis for interpreting peak shifts and splitting as signatures of stiffness evolution, rather than purely phenomenological indicators.

Crucially, separating the solid-solution trend $f_{\mathrm{SAW}}(x_{\mathrm{Cr}})$ from the spinodal conversion trend $f_{\mathrm{SAW}}(\phi)$ resolves a key ambiguity in experimental interpretation. While composition alone can shift SAW frequencies, progressive spinodal decomposition introduces additional stiffening that cannot be captured by composition-dependent elasticity alone. The combined framework therefore enables SAW peak splitting to be interpreted as a mechanistically grounded indicator of microstructural evolution toward embrittlement, rather than as a purely compositional effect.""",
     r"""While composition alone can shift SAW frequencies, progressive spinodal decomposition introduces additional stiffening that cannot be captured by composition-dependent elasticity alone. Separating the solid-solution trend $f_{\mathrm{SAW}}(x_{\mathrm{Cr}})$ from the spinodal conversion trend $f_{\mathrm{SAW}}(\phi)$ therefore allows SAW peak splitting to be read as a sign of microstructural evolution toward embrittlement."""),
    (r"""This separation clarifies why SAW peak splitting emerges as a robust, non-destructive diagnostic of spinodal decomposition in CASS and related alloys, and suggests a general strategy for extending ultrasonic embrittlement detection to other phase-separating structural materials.""",
     r"""The same approach may extend ultrasonic embrittlement detection to other phase-separating structural materials."""),
    # Conclusions: the findings instead of generic statements
    (r"""Free-energy curvature provides a direct stability criterion that identifies spinodally unstable composition windows. DFT-derived cubic elasticity quantifies composition-dependent stiffening trends and supplies inputs for elastodynamic prediction. Combining these elements yields experimentally actionable relationships for interpreting picosecond ultrasonics across both composition and decomposition-state variables.""",
     r"""Spinodal decomposition raises the Rayleigh SAW frequency by up to 5.7\% in MD, most strongly within the spinodally unstable window that contains Fe--36Cr, by stiffening the tetragonal shear modulus $C'$ through the Cr-rich $\alpha'$ phase. Only the first-principles $\alpha'$ tensors reach the frequency range unallowed for pristine CASS. The same shear stiffening lowers Pugh's ratio $K/G$, linking the acoustic signal to embrittlement."""),
]


def main():
    for v in ("", "_cubes", "_trend", "_trend_cubes"):
        src = MS / f"finalpaper_qe_mechanism_v43_validation{v}.tex"
        dst = MS / f"finalpaper_qe_mechanism_v44_validation{v}.tex"
        assert not dst.exists(), dst
        t = src.read_text()
        i = t.index(r"\section*{Supplementary Material}")
        main_, sm = t[:i], t[i:]
        for a, b in EDITS:
            assert main_.count(a) == 1, (v, a[:70])
            main_ = main_.replace(a, b)
        dst.write_text(main_ + sm)
        print("wrote", dst.name)


if __name__ == "__main__":
    main()
