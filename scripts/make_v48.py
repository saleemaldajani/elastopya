#!/usr/bin/env python3
"""Paper v48 (2026-11-03): v47 (submitted v44 in Scripta format) with two main-text edits (user request):
  - the repository is cited as a numbered [software] reference (elastopya2026) at its first mention (Methods), as the
    Scripta guide for authors asks (Research data, Option B: cite and link the dataset); the standalone closing sentence
    with the URL and the dangling "(see Data Availability)" left from v47 go; the data statement goes in Editorial Manager;
  - full postal address with the country (guide for authors);
  - the one-sentence paragraph "The same approach may extend ..." joins the end of the conclusions paragraph.
Main text only; the SM is unchanged. v47 is kept.

  make_v48.py     writes finalpaper_qe_mechanism_v48_validation[_cubes|_trend|_trend_cubes].tex from v47
"""
from pathlib import Path

MS = Path(__file__).resolve().parent.parent / "manuscript"
URL = r"\url{https://github.com/saleemaldajani/elastopya}"

EDITS = [   # applied to the main text only; each must occur exactly once there
    (r"""All LAMMPS inputs, random seeds, and potential files are provided in the accompanying repository.""",
     r"""All LAMMPS inputs, random seeds, and potential files are provided in the accompanying repository \cite{elastopya2026}."""),
    # guide for authors: full postal address including the country
    ("Cambridge, MA 02139}", "Cambridge, MA 02139, USA}"),
    (r"""All figures are generated from tabulated inputs by scripts provided in the accompanying repository (see Data Availability).""",
     r"""All figures are generated from tabulated inputs by scripts provided in the same repository, which also holds the first-principles inputs and outputs."""),
    ("\n\nThe same approach may extend ultrasonic embrittlement detection to other phase-separating structural materials.\n", "\n"),
    (r"""The same shear stiffening lowers Pugh's ratio $K/G$, linking the acoustic signal to embrittlement.""",
     r"""The same shear stiffening lowers Pugh's ratio $K/G$, linking the acoustic signal to embrittlement. This approach may extend ultrasonic embrittlement detection to other phase-separating structural materials."""),
    (rf"""The tabulated data, the atomistic and first-principles simulation inputs and outputs, and the scripts that generate all figures are publicly available at {URL}.""", ""),
]


def main():
    for v in ("", "_cubes", "_trend", "_trend_cubes"):
        src = MS / f"finalpaper_qe_mechanism_v47_validation{v}.tex"
        dst = MS / f"finalpaper_qe_mechanism_v48_validation{v}.tex"
        assert not dst.exists(), dst
        t = src.read_text()
        i = t.index(r"\section*{Supplementary Material}")
        main_, sm = t[:i], t[i:]
        for a, b in EDITS:
            assert main_.count(a) == 1, (v, a[:60], main_.count(a))
            main_ = main_.replace(a, b)
        assert "github.com/saleemaldajani/elastopya" not in main_ and "Data Availability" not in main_
        dst.write_text(main_ + sm)
        print("wrote", dst.name)


if __name__ == "__main__":
    main()
