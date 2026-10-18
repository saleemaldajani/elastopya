#!/usr/bin/env python3
"""Build the journal-submission folder (2026-10-18): main text and Supplementary Material as separate LaTeX documents,
each self-contained (own bibliography), from the combined paper, with the figures, tables and bibliography they use.

Cross-document references are resolved to the numbers of the combined paper (read from its .aux): main text -> SM
figures (S1, S13), SM -> main-text Fig. 1 and Sec. 2.3. SM equations are numbered S1, S2, ... as its figures and
tables already are.

  make_submission.py <combined .tex> <combined .aux> <out dir>
"""
import re
import shutil
import sys
from pathlib import Path

MS = Path(__file__).resolve().parent.parent / "manuscript"


def main():
    src, aux, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    labels = dict(re.findall(r"\\newlabel\{([^}]*)\}\{\{([^}]*)\}", aux.read_text()))
    t = src.read_text().split("\n")
    i_doc = next(i for i, l in enumerate(t) if l.startswith("\\begin{document}"))
    i_sm = next(i for i, l in enumerate(t) if l.startswith("\\section*{Supplementary Material}"))
    i_bib = next(i for i, l in enumerate(t) if l.startswith("\\bibliographystyle"))
    i_end = next(i for i, l in enumerate(t) if l.startswith("\\end{document}"))
    pre = t[: i_doc + 1]
    main_body = t[i_doc + 1: i_bib + 2]                      # through \bibliography{references}
    sm_body = t[i_sm + 1: i_end]
    main_txt, sm_txt = "\n".join(main_body), "\n".join(sm_body)
    main_labels = set(re.findall(r"\\label\{([^}]*)\}", main_txt))
    sm_labels = set(re.findall(r"\\label\{([^}]*)\}", sm_txt))
    for tabfile in re.findall(r"\\input\{([^}]*)\}", sm_txt):
        sm_labels |= set(re.findall(r"\\label\{([^}]*)\}", (MS / f"{tabfile}.tex").read_text()))

    def resolve(txt, foreign):
        return re.sub(r"\\ref\{([^}]*)\}", lambda m: labels[m.group(1)] if m.group(1) in foreign else m.group(0), txt)
    main_txt = resolve(main_txt, sm_labels)
    sm_txt = resolve(sm_txt, main_labels)
    sm_txt = sm_txt.replace("Eq.~(2)", "Eq.~(S2)")
    if "\\renewcommand{\\theequation}" in sm_txt:
        sm_head_eq = ""
    else:
        sm_head_eq = "\\renewcommand{\\theequation}{S\\arabic{equation}}\n\\setcounter{equation}{0}\n"
    sm_head = r"""
\begin{center}
{\Large\textbf{Supplementary Material}}\\[6pt]
{\large First principles evidence of stiffening during spinodal decomposition in FeCr}\\[6pt]
Saleem A. Al Dajani and Michael P. Short\\
Department of Nuclear Science and Engineering, Massachusetts Institute of Technology, Cambridge, MA 02139, USA
\end{center}
"""
    out.mkdir(parents=True, exist_ok=True)
    (out / "main.tex").write_text("\n".join(pre) + "\n" + main_txt + "\n\n\\end{document}\n")
    (out / "supplement.tex").write_text("\n".join(pre) + "\n" + sm_head + sm_head_eq + sm_txt
                                        + "\n\n\\bibliographystyle{elsarticle-num}\n\\bibliography{references}\n\n\\end{document}\n")
    # assets: figures, tables, bibliography
    used = set()
    for txt in (main_txt, sm_txt):
        used |= set(re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}", txt))
    for f in sorted(used):
        p = MS / f
        if not p.exists():
            sys.exit(f"missing figure {f}")
        (out / f).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(p, out / f)
    for tabfile in re.findall(r"\\input\{([^}]*)\}", sm_txt + main_txt):
        (out / tabfile).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(MS / f"{tabfile}.tex", out / f"{tabfile}.tex")
    shutil.copy2(MS / "references.bib", out / "references.bib")
    print("main refs resolved:", sorted(set(re.findall(r"\\ref\{([^}]*)\}", "\n".join(main_body))) & sm_labels))
    print("SM refs resolved:", sorted(set(re.findall(r"\\ref\{([^}]*)\}", "\n".join(sm_body))) & main_labels))
    print("figures:", len(used))


if __name__ == "__main__":
    main()
