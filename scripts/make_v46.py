#!/usr/bin/env python3
"""Paper v46 (2026-11-01): post-submission working copy of v44 (the version submitted to Scripta Materialia), for the
revision (user request):
  - an empty Funding section with a visible placeholder, to be filled during the revision;
  - additional affiliations for S. A. Al Dajani (KAUST Applied Physics Program; MIT Civil and Environmental Engineering),
    taken from the cover letter's signature block, marked TODO in the source to confirm wording and author assignment;
  - corresponding author as an elsarticle \\cortext (same e-mail).
(v45 is the parked draft without the Rae citation, see manuscript/v45_not_used/.)

  make_v46.py     writes finalpaper_qe_mechanism_v46_validation[_cubes|_trend|_trend_cubes].tex from v44
"""
from pathlib import Path

MS = Path(__file__).resolve().parent.parent / "manuscript"

OLD_AUTH = r"""\author[1]{Saleem A. Al Dajani}
\author[1]{Michael P. Short}

\address{
Department of Nuclear Science and Engineering \\
Massachusetts Institute of Technology\\
77 Massachusetts Avenue \\
Cambridge, MA 02139}

\address[1]
{
Corresponding Author: \href{mailto:sdajani@mit.edu}{\texttt{sdajani@mit.edu}}
}"""

NEW_AUTH = r"""% TODO (revision): confirm the additional affiliations (2, 3), their wording, and which authors they belong to.
\author[1,2,3]{Saleem A. Al Dajani\corref{cor1}}
\author[1]{Michael P. Short}

\address[1]{Department of Nuclear Science and Engineering, Massachusetts Institute of Technology, 77 Massachusetts Avenue, Cambridge, MA 02139, USA}
\address[2]{Applied Physics Program, King Abdullah University of Science and Technology (KAUST), 4700 KAUST, Thuwal 23955-6900, Kingdom of Saudi Arabia}
\address[3]{Department of Civil and Environmental Engineering, Massachusetts Institute of Technology, 77 Massachusetts Avenue, Cambridge, MA 02139, USA}

\cortext[cor1]{Corresponding author: \href{mailto:sdajani@mit.edu}{\texttt{sdajani@mit.edu}}}"""

OLD_ACK = r"""\section*{Acknowledgements}"""
NEW_ACK = r"""\section*{Funding}
% TODO (revision): funding sources of both authors (grant numbers).
\textit{[Funding statement to be added.]}

\section*{Acknowledgements}"""


def main():
    for v in ("", "_cubes", "_trend", "_trend_cubes"):
        src = MS / f"finalpaper_qe_mechanism_v44_validation{v}.tex"
        dst = MS / f"finalpaper_qe_mechanism_v46_validation{v}.tex"
        assert not dst.exists(), dst
        t = src.read_text()
        for a, b in ((OLD_AUTH, NEW_AUTH), (OLD_ACK, NEW_ACK)):
            assert t.count(a) == 1, (v, a[:40])
            t = t.replace(a, b)
        dst.write_text(t)
        print("wrote", dst.name)


if __name__ == "__main__":
    main()
