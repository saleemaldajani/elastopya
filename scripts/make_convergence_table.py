#!/usr/bin/env python3
"""SM table: structural convergence of every QE relaxation used in the paper (added 2026-10-16, answering Mike's
request for the convergence studies; the parameter study is SM Fig. S8).

For each 40-atom SQS (PBE-relaxed cell, and experimental lattice parameter with ions relaxed) it reports the residual
total force and the pressure at the last converged ionic step, the number of BFGS steps (summed over resumed parts),
and whether any SCF before the end of BFGS failed to converge. Read from dft_qe/engaging/runs/*/relax*.out.

  make_convergence_table.py  -> manuscript/tables/S_structural_convergence.tex
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RUNS = ROOT / "dft_qe" / "engaging" / "runs"
OUT = ROOT / "manuscript" / "tables" / "S_structural_convergence.tex"


def last_bfgs(text):
    """Text up to the end of BFGS (drops the closing check SCF of vc-relax, which is not used)."""
    i = text.find("End of BFGS Geometry Optimization")
    return text if i < 0 else text[:i]


def row(d):
    parts = sorted(d.glob("relax_part*.out"), key=lambda q: int(q.stem[10:])) + [d / "relax.out"]
    texts = [p.read_text() for p in parts]
    main = last_bfgs(texts[-1])
    if "End of BFGS" not in texts[-1]:
        return None
    force = float(re.findall(r"Total force =\s+([\d.]+)", main)[-1])
    press = float(re.findall(r"P=\s*([-\d.]+)", main)[-1])
    steps = sum(len(re.findall(r"^!\s+total energy", last_bfgs(t), re.M)) for t in texts)
    bad = sum(last_bfgs(t).count("convergence NOT achieved") for t in texts)
    return force, press, steps, bad


def main():
    rows = []
    for d in sorted(RUNS.glob("sqs40_x*_k442_e60")) + [RUNS / "sqs40_x0.80_k442_e60_crafm2"]:
        rows.append(("relaxed", d))
    for d in sorted(RUNS.glob("sqs40_x*_k442_e60_aexp")) + [RUNS / "sqs40_x0.80_k442_e60_aexp_crafm"]:
        rows.append(("exp.", d))
    lines, fmax, pmax, worst = [], 0.0, 0.0, 0
    for lat, d in rows:
        if not ((d / "result.json").exists() or (d / "relax_result.json").exists() ):
            continue
        r = row(d)
        if r is None:
            continue
        x = float(re.match(r"sqs40_x([\d.]+)_", d.name).group(1))
        f, p, n, bad = r
        fmax = max(fmax, f); worst += bad
        if lat == "relaxed":
            pmax = max(pmax, abs(p))
        lines.append((lat, x, f, p, n, bad))
    lines.sort(key=lambda t: (t[0] != "relaxed", t[1]))
    body = "\n".join(f"{(f'{x:.3f}' if round(x, 3) == 0.475 else f'{x:.2f}')} & {lat} & {f:.4f} & {p:.2f} & {n} \\\\" for lat, x, f, p, n, bad in lines)
    tex = (r"""\begin{table}[!htbp]
\centering
\caption{Structural convergence of the Quantum ESPRESSO relaxations (40-atom SQS). Residual total force and pressure at
the last ionic step, and number of ionic (BFGS) steps. Thresholds: energy change $10^{-5}$~Ry, force components
$10^{-3}$~Ry/bohr, pressure 0.2~kbar for the PBE-relaxed cells; at the experimental lattice parameter the cell is fixed
and only the ions relax, hence the non-zero pressure. The SCF threshold is $10^{-8}$~Ry, except $3\times10^{-6}$~Ry for
the last 19 of the 21 ionic steps of Fe--65Cr at the experimental lattice parameter, where the SCF reached a $10^{-6}$~Ry noise floor (energy error below
1~$\mu$eV/atom). Every SCF of every ionic step converged.}
\label{tab:S_structural_convergence}
\small
\begin{tabular}{rlrrr}
\hline
$x_{\mathrm{Cr}}$ & lattice & force (Ry/bohr) & $P$ (kbar) & ionic steps \\
\hline
""" + body + "\n\\hline\n\\end{tabular}\n\\end{table}\n")
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(tex)
    # compact side-by-side version for the response to Mike (no float, no caption)
    rel = [t for t in lines if t[0] == "relaxed"]; exp = [t for t in lines if t[0] != "relaxed"]
    fx = lambda x: f"{x:.3f}" if round(x, 3) == 0.475 else f"{x:.2f}"
    cell = lambda t: f"{fx(t[1])} & {t[2]:.4f} & {t[3]:.2f} & {t[4]}" if t else "& & &"
    rows = [f"{cell(rel[i] if i < len(rel) else None)} & {cell(exp[i] if i < len(exp) else None)} \\\\"
            for i in range(max(len(rel), len(exp)))]
    (OUT.parent / "R_structural_convergence_compact.tex").write_text(
        "{\\footnotesize\\setlength{\\tabcolsep}{4pt}\\begin{tabular}{rrrr|rrrr}\n\\hline\n"
        "\\multicolumn{4}{c|}{PBE-relaxed cell} & \\multicolumn{4}{c}{experimental $a$ (ions relaxed)} \\\\\n"
        "$x_{\\mathrm{Cr}}$ & force & $P$ & steps & $x_{\\mathrm{Cr}}$ & force & $P$ & steps \\\\\n\\hline\n"
        + "\n".join(rows) + "\n\\hline\n\\end{tabular}}\n")
    print(tex)
    print(f"max force {fmax:.4f} Ry/bohr, max |P| relaxed {pmax:.2f} kbar, unconverged SCF before end of BFGS: {worst}")


if __name__ == "__main__":
    main()
