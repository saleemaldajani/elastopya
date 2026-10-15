#!/usr/bin/env python3
"""Rebuild the Supplementary Material of the validation paper in topic order (2026-10-15).

Every figure and table is moved verbatim (caption and label unchanged). Each new section starts on a new page with its
text first and its figures after, pinned to the section with \\FloatBarrier, so figures no longer split unrelated text.
The original-draft sections S1-S2 stay first and unchanged.

  reorganize_sm.py <in.tex> <out.tex>
"""
import re
import sys
from pathlib import Path

src, dst = Path(sys.argv[1]), Path(sys.argv[2])
s = src.read_text()


def block(label):
    """The float environment (figure/table) that carries \\label{label}."""
    i = s.index("\\label{" + label + "}")
    a = max(s.rindex("\\begin{figure}", 0, i) if "\\begin{figure}" in s[:i] else -1,
            s.rindex("\\begin{table}", 0, i) if "\\begin{table}" in s[:i] else -1)
    env = "table" if s.startswith("\\begin{table}", a) else "figure"
    b = s.index("\\end{" + env + "}", i) + len("\\end{" + env + "}")
    out = s[a:b]
    return re.sub(r"\\begin\{(figure|table)\}\[[^\]]*\]", r"\\begin{\1}[!htbp]", out, count=1)


def between(start, end):
    a = s.index(start); b = s.index(end, a)
    return s[a:b]


sm0 = s.index("\\subsection*{S3. First-principles calculations with Quantum ESPRESSO}")
head = s[:sm0].rstrip()
head = head.replace("Sections S3 and S4 contain the figures added in this revision.",
                    "Sections S3--S9 contain the material added in this revision.")
methods = between("\\paragraph{Methods.}", "\\begin{figure}").strip()
val_text = between("Every quantity shown in Section S3 was", "\\begin{table}").strip().replace(
    "Every quantity shown in Section S3 was", "Every quantity shown in Sections S4--S7 was")
val_table = block("tab:validation")
val_concl = between("The electronic structure and magnetism agree well with the literature", "\\begin{figure}").strip()
ett_text = between("The Fe-rich anomaly of the elastic constants (Fig.~\\ref{fig:S_fe_rich})", "\\begin{figure}").strip()
mech = s[s.index("\\subsection*{S4. Proposed mechanism figures}"):s.index("\\end{document}")]
mech = mech.replace("\\subsection*{S4. Proposed mechanism figures}", "\\subsection*{S9. Proposed mechanism figures}")
tail = s[s.index("\\end{document}"):]

dos = block("fig:S_dos").replace("supplementary/S_dos_sqs.pdf", "supplementary/S_dos_sqs_wide_v2.pdf")
dos = re.sub(r"\\upd\{the Fe-80Cr DOS and the wide-window DOS above the Fermi level will be added\.\}",
             r"The energy window extends to 4 eV above $E_{\\mathrm F}$ (1.6 times the occupied bands).", dos)
bulk = block("fig:V_bulk").replace("additional QE series of Section S3", "additional QE series of Section S4")

FB = "\n\\FloatBarrier\n"


def sec(title, *parts):
    return "\\clearpage\n\\subsection*{" + title + "}\n" + "\n\n".join(p.strip() for p in parts if p) + FB


S3 = sec("S3. First-principles calculations: methods and overview", methods,
         r"""\paragraph{Overview.} The additional first-principles material is organised by topic: elastic constants at the PBE-relaxed and the room-temperature experimental lattice parameter (Section S4); electronic structure, including the unfolded Fermi surfaces and band structures of the Fe-rich alloys (S5); the electronic topological transition that links them to the Fe-rich elastic anomaly (S6); lattice dynamics (S7); a validation of every computed quantity against published data and experiment, including the thermal-transport estimates of Section S1 (S8); and the proposed mechanism figures (S9). Each section gives its text first and its figures after.""",
         r"""\upd{still running: unfolded Fermi surfaces and band structures at 5 and 10 at.-pct Cr; Fe-80Cr elastic constants at the PBE-relaxed lattice parameter; the Cr spin- and charge-density-wave supercell, the Cr nesting function and spin-density maps of the SQS; SQS phonon DOS with the vibrational entropy of mixing.}""")

S4 = sec("S4. Elastic constants",
         r"""Figure~\ref{fig:S_fe_rich} compares the Fe-rich elastic constants with EMTO-CPA and MD. Figure~\ref{fig:S_qe_aexp} repeats the QE calculations at the room-temperature experimental lattice parameter (Vegard between Fe 2.8665 and Cr 2.8846~\AA; cubic cell fixed, ions relaxed), the convention of the EMTO-CPA reference data \cite{razumovskiy2011first}. PBE relaxes all compositions to about 11.4~\AA$^3$/atom, about 4\% below experiment; the bulk modulus is the property most sensitive to this (Section S8).""",
         r"""At the experimental lattice parameter $f_{\mathrm{SAW}}$ rises from 446~MHz for Fe (EMTO-CPA: 449~MHz) to 485 at 20\%, flattens at 497, 500 and 499~MHz between 30 and 40\%~Cr (EMTO-CPA: 503, 488, 490~MHz), and reaches 513 at 50\% (513), 542 at 80\% and 580~MHz for AFM Cr (602). The bulk modulus has a local minimum at 40\%~Cr (163~GPa) where the Cr moments collapse, and the shear modulus flattens there. The experimental lattice parameter thus reproduces the magnitude of $f_{\mathrm{SAW}}$ and a magnetic anomaly near Fe--36Cr, but not the Fe-rich maximum near 20\%~Cr of the EMTO-CPA curve, so it gives a plateau rather than a minimum. At the PBE-relaxed volume $f_{\mathrm{SAW}}$ rises steadily from 504 to 591~MHz. The bulk-modulus drops at 80\% and 100\%~Cr reflect antiferromagnetic (Cr) and ferrimagnetic (Fe--80Cr, Cr antiparallel to Fe, net 0.31~$\mu_{\mathrm B}$/atom) order, which collapses $C_{12}$ as in pure Cr.""",
         block("fig:S_fe_rich"),
         r"""\begin{figure}[!htbp]
\centering
\includegraphics[width=\textwidth]{figures/kdense_scripta/supplementary/S_qe_aexp_progress.pdf}
\caption{QE 40-atom SQS and pure elements at the PBE-relaxed (open, dashed) and the room-temperature experimental (filled) lattice parameter, against EMTO-CPA \cite{razumovskiy2011first}: (a) $f_{\mathrm{SAW}}$ at $\lambda=6.4\,\mu$m, (b) bulk modulus, (c) VRH shear modulus. The shaded band marks 30--40\%~Cr. Pure Cr is AFM; Fe--80Cr at the experimental lattice parameter converged to a ferrimagnetic state.}
\label{fig:S_qe_aexp}
\end{figure}""")

S5 = sec("S5. Electronic structure",
         r"""Fermi surfaces and band structures of the pure metals (Figs.~\ref{fig:S_fermi_Fe}--\ref{fig:S_bands_aexp}) set the reference for the alloys. The unfolded Fermi surfaces and band structures of the Fe-rich SQS (Figs.~\ref{fig:S_fermi_unfolded} and \ref{fig:S_bands_unfolded}) follow how Cr changes the spin-up Fermi surface; the projected densities of states of the SQS (Fig.~\ref{fig:S_dos}) show the composition dependence of the Fe and Cr states and of the Fermi-level DOS.""",
         block("fig:S_fermi_Fe"), block("fig:S_fermi_pure"), block("fig:S_bands_pure"), block("fig:S_bands_aexp"),
         block("fig:S_fermi_unfolded"), block("fig:S_bands_unfolded"), dos)

S6 = sec("S6. Electronic topological transition and the Fe-rich elastic anomaly", ett_text, block("fig:S_ett"))

S7 = sec("S7. Lattice dynamics",
         r"""Phonon dispersions of Fe, non-magnetic Cr and antiferromagnetic Cr are shown in Fig.~\ref{fig:S_phonons}, the $k$-point convergence that removes spurious imaginary modes of non-magnetic Cr in Fig.~\ref{fig:S_phonon_conv}, and the comparison with neutron and nuclear inelastic scattering in Fig.~\ref{fig:S_phonons_exp}. Cr is computed with density-functional perturbation theory; Fe by finite displacements, for the reason given in the captions.""",
         block("fig:S_phonons"), block("fig:S_phonon_conv"), block("fig:S_phonons_exp"))

THERMAL = r"""\paragraph{Thermal transport (Section S1).} The estimates of Section S1 depend only on $E$, $\rho$ and the molar mass and were evaluated for every tensor of this work (Fig.~\ref{fig:S_thermal}). The elastic Debye temperature of Fe is 450~K (QE, experimental lattice parameter), 453~K (EMTO-CPA) and 472~K (MD), against 470~K from calorimetry \cite{kittel2005introduction}; for Cr, 585--605~K against 630~K. Debye temperatures from the second moment of the phonon DOS, a definition that runs lower, agree with experiment on the same definition: 415~K for our Fe against 417--431~K from NRIXS, and 466~K for our AFM Cr against 460~K from the neutron-based DOS of Fultz \textit{et al.} The minimum thermal conductivity of Eq.~(2), $\kappa_{\min}=1.15$--1.6~W/(m\,K), gives $D_{\min}=\kappa_{\min}/\rho c_p=0.32$--$0.48\times10^{-6}$~m$^2$/s, an order of magnitude below the thermal diffusivity measured by TGS in CF8 ($3.5\times10^{-6}$~m$^2$/s, Fig.~\ref{fig5}d). This is expected, since $\kappa_{\min}$ is a lower bound on the lattice contribution and electrons carry most of the heat in a metal. The relative changes do not depend on the prefactor. On full decomposition of MD Fe--36Cr, $E$ rises by 9.5\% and the S1 scaling $D\propto E^{1/2}$ gives $\Delta D/D=+4.7\%$, close to the $+4.0\%$ measured after aging at 290~$^\circ$C. The larger increases measured after aging at 330--400~$^\circ$C ($+20$ to $+29\%$) exceed what lattice stiffening alone can produce; they point to the electronic contribution, which the S1 estimate does not include. \upd{Fig. S7e,f of the original draft evaluated kappa with Planck's constant (times 1e11) in place of k\_B and m = 2, which scales kappa and D by 7.6 relative to Eq. (2); the trends are unchanged but the magnitudes need correcting.}"""

THERMAL_FIG = r"""\begin{figure}[!htbp]
\centering
\includegraphics[width=\textwidth]{figures/kdense_scripta/supplementary/S_thermal.pdf}
\caption{Validation of the thermal-transport estimates of Section S1. (a) Debye temperature: elastic (Anderson) from every tensor of this work and EMTO-CPA; from the second moment of the phonon DOS (Fultz \textit{et al.}\ \cite{fultz1995phonon}, NRIXS \cite{mauger2014nonharmonic}, ours; this definition runs about 10\% below the elastic one); calorimetric values of Fe and Cr \cite{kittel2005introduction}. (b) $D_{\min}=\kappa_{\min}/\rho c_p$ from Eq.~(2) of Section S1 for each tensor set; it is a lower bound on the lattice contribution. (c) Relative change of the thermal diffusivity: measured by TGS in CF8 after aging (circles, Fig.~\ref{fig5}d) and predicted by the S1 scaling for full decomposition of MD Fe--36Cr and Fe--50Cr (dashed). Script: \texttt{dft\_qe/thermal\_check.py}.}
\label{fig:S_thermal}
\end{figure}"""

V_FIGS = [block(l) if l != "fig:V_bulk" else bulk for l in
          ("fig:V_bulk", "fig:V_moments", "fig:V_nef", "fig:V_bands", "fig:V_phdisp_Cr", "fig:V_phdos_Cr", "fig:V_phdisp_Fe", "fig:V_phdos_Fe")]
S8 = sec("S8. Validation against literature and experiment", val_text, val_table, val_concl, THERMAL, *V_FIGS, THERMAL_FIG)

out = head + "\n\n" + FB + S3 + S4 + S5 + S6 + S7 + S8 + mech.rstrip() + "\n\n" + tail
out = out.replace("SM Section S3.}", "SM Section S4.}")
if "\\usepackage{placeins}" not in out:
    out = out.replace("\\usepackage{booktabs}", "\\usepackage{booktabs}\n\\usepackage{placeins}", 1)
for lab in ("fig:S_fe_rich", "fig:S_dos", "fig:S_ett", "fig:S_phonons", "fig:V_bulk", "fig:S_mechanism", "fig:S_embrittlement", "fig5"):
    assert out.count("\\label{" + lab + "}") == 1, lab
dst.write_text(out)
print("wrote", dst, "| figures:", out.count("\\begin{figure}") + out.count("\\begin{sidewaysfigure}"), "(in:", s.count("\\begin{figure}") + s.count("\\begin{sidewaysfigure}"), ")")
