#!/usr/bin/env python3
"""Paper v45 (2026-11-01): v44 without the citation of Rae et al. (arXiv 2026), at the user's request: it was cited only
for why Fe-36Cr is singled out ("the alloy in which spinodal decomposition was recently detected through TGS
stiffening"); no data or method comes from it. Fe-36Cr is justified in Sec. 3.1 by the 748 K CALPHAD spinodal of Xiong
et al. Changes: Sec. 3.1 sentence (main text and its copy in SM S11), Fig. 3 caption. v44 is kept.

  make_v45.py     writes finalpaper_qe_mechanism_v45_validation[_cubes|_trend|_trend_cubes].tex from v44
"""
from pathlib import Path

MS = Path(__file__).resolve().parent.parent / "manuscript"

EDITS = [   # (old, new, expected count in the whole file)
    (r"""Both windows contain the Fe--36Cr alloy in which spinodal decomposition was recently detected through TGS stiffening \cite{rae2026tgs}.""",
     r"""Both windows contain Fe--36Cr.""", 2),                    # Sec. 3.1 and its copy in SM S11
    (r"""the star marks Fe--36Cr \cite{rae2026tgs} and the arrows""",
     r"""the star marks Fe--36Cr and the arrows""", 1),
]


def main():
    for v in ("", "_cubes", "_trend", "_trend_cubes"):
        src = MS / f"finalpaper_qe_mechanism_v44_validation{v}.tex"
        dst = MS / f"finalpaper_qe_mechanism_v45_validation{v}.tex"
        assert not dst.exists(), dst
        t = src.read_text()
        for a, b, n in EDITS:
            assert t.count(a) == n, (v, a[:60], t.count(a))
            t = t.replace(a, b)
        assert "rae2026tgs" not in t
        dst.write_text(t)
        print("wrote", dst.name)


if __name__ == "__main__":
    main()
