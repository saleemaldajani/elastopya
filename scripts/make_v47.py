#!/usr/bin/env python3
"""Paper v47 (2026-11-03): the submitted v44 reformatted to Scripta Materialia style after the journal sent it back
(SMM-26-1622): no section headings in the main text except Abstract, Acknowledgements and References. The text is the
submitted v44 unchanged apart from (i) the removed headings (the Data Availability sentence stays as a plain closing
paragraph) and (ii) references to main-text section numbers, which no longer exist, reworded to equations, figures or
"the main text" (in the main text and in the SM). v44 and v46 are kept.

  make_v47.py     writes finalpaper_qe_mechanism_v47_validation[_cubes|_trend|_trend_cubes].tex from v44
"""
import re
from pathlib import Path

MS = Path(__file__).resolve().parent.parent / "manuscript"

HEADINGS = [r"\section{Introduction and Background}", r"\section{Methods}", r"\subsection{Elastodynamic Methods}",
            r"\subsection{DFT Methods}", "\\subsection{Atomistic simulations of decomposed Fe--Cr}\n\\label{sec:md}",
            "\\subsection{Analysis of DFT Results}\n\\label{sec:analysis}", r"\section{Results}",
            r"\subsection{Spinodal stability from demixing free-energy curvature}",
            r"\subsection{Elastic stiffening and elastodynamic mapping to Rayleigh SAW observables}",
            r"\subsection{Separating composition trends from spinodal conversion trends}",
            r"\section{Discussion}", r"\section{Conclusions}", r"\section*{Data Availability}"]

EDITS = [   # (old, new, expected count in the whole file: main text + SM)
    # VRH average: point to the moduli equation (main text and its copy in SM S11)
    (r"""average of each elastic tensor (Sec.~\ref{sec:analysis}).""", r"""average of each elastic tensor (Eq.~\ref{eqn:moduli}).""", 2),
    (r"""and the full text of Sections 2--3 before condensation are given in SM Section S11.""",
     r"""and a longer version of the methods and results are given in SM Section S11.""", 1),
    # figure captions
    (r"""($25^3$ \textit{bcc} cells; Sec.~\ref{sec:md})""", r"""($25^3$ \textit{bcc} cells)""", None),          # 2D-map styles
    (r"""(three outer faces, corner sublattice; Sec.~\ref{sec:md})""", r"""(three outer faces, corner sublattice)""", None),  # cube styles
    (r"""three elsewhere; Sec.~\ref{sec:md})""", r"""three elsewhere)""", 1),
    (r"""from the MD (EAM) decomposition simulations (Sec.~\ref{sec:md}):""", r"""from the MD (EAM) decomposition simulations:""", 1),
    # the surface boundary-value problem (main text and S11)
    (r"""through the surface boundary-value problem of Sec.~2.1.""", r"""through the surface boundary-value problem of Eq.~\ref{eqn:rayleigh}.""", 2),
    # SM references to main-text sections
    (r"""Three checks constrain the mechanism stated in Sec.~3.2.""", r"""Three checks constrain the mechanism stated in the main text.""", 1),
    (r"""The conversion procedure (Sec.~\ref{sec:md}) is""", r"""The conversion procedure (main text) is""", 1),
    (r"""to 0.95 (Sec.~\ref{sec:md}; Fig.~\ref{fig:fsaw_combined}b)""", r"""to 0.95 (Fig.~\ref{fig:fsaw_combined}b)""", 1),
    # SM S11: section numbers no longer exist in the main text
    (r"""\subsection*{S11. Full text of Sections 2--3}""", r"""\subsection*{S11. Full text of the methods and results}""", 1),
    (r"""The main text condenses Sections 2--3 to the journal's word limit. Their full text before condensation follows; equation, figure and section numbers refer to the main text.""",
     r"""The main text condenses the methods and results to the journal's word limit. Their full text before condensation follows; equation and figure numbers refer to the main text.""", 1),
]
S11_TITLES = [("2.1 Elastodynamic Methods.", "Elastodynamic methods."), ("2.2 DFT Methods.", "DFT methods."),
              ("2.3 Atomistic simulations of decomposed Fe--Cr.", "Atomistic simulations of decomposed Fe--Cr."),
              ("3.1 Spinodal stability from demixing free-energy curvature.", "Spinodal stability from demixing free-energy curvature."),
              ("3.2 Elastic stiffening and elastodynamic mapping to Rayleigh SAW observables.", "Elastic stiffening and elastodynamic mapping to Rayleigh SAW observables."),
              ("3.3 Separating composition trends from spinodal conversion trends.", "Separating composition trends from spinodal conversion trends.")]


def main():
    for v in ("", "_cubes", "_trend", "_trend_cubes"):
        src = MS / f"finalpaper_qe_mechanism_v44_validation{v}.tex"
        dst = MS / f"finalpaper_qe_mechanism_v47_validation{v}.tex"
        assert not dst.exists(), dst
        t = src.read_text()
        for a, b, n in EDITS:
            assert (t.count(a) <= 1) if n is None else (t.count(a) == n), (v, a[:60], t.count(a))
            t = t.replace(a, b)
        for a, b in S11_TITLES:
            a2, b2 = rf"\paragraph{{{a}}}", rf"\paragraph{{{b}}}"
            assert t.count(a2) == 1, (v, a2)
            t = t.replace(a2, b2)
        i = t.index(r"\section*{Supplementary Material}")
        main_, sm = t[:i], t[i:]
        for h in HEADINGS:
            assert main_.count(h) == 1, (v, h)
            main_ = main_.replace(h + "\n", "").replace(h, "")
        t = main_ + sm
        left = re.findall(r"\\(?:sub)*section\*?\{[^}]*\}", t[:t.index(r"\section*{Supplementary Material}")])
        assert left == [r"\section*{Acknowledgements}"], left
        assert "sec:md" not in t and "sec:analysis" not in t
        own = re.sub(r"\(Sec\.~3\.7 and Fig\.~8 therein\)|\(Sec\.~3\.8, Table~4\)", "", t)   # sections of the Acta paper
        assert not re.search(r"Sec\.~[23]\.[0-9]", own), re.findall(r".{60}Sec\.~[23]\.[0-9].{20}", own)
        dst.write_text(t)
        print("wrote", dst.name)


if __name__ == "__main__":
    main()
