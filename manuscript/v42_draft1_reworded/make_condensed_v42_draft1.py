#!/usr/bin/env python3
"""Paper v42 (2026-11-01): main text condensed to Scripta's 2,500-word limit (user request), figures and captions unchanged.

Every condensed paragraph of Secs. 1-3 is replaced by a shorter one, and the full original text of those sections is
appended, verbatim, as SM Section S11, so no content is lost (labels removed, equations unnumbered). The Discussion and
Conclusions (Mike's wording) are unchanged. v41 stays as the long version.

  make_condensed_v42.py      writes finalpaper_qe_mechanism_v42_validation[_cubes|_trend|_trend_cubes].tex from v41
"""
import re
from pathlib import Path

MS = Path(__file__).resolve().parent.parent / "manuscript"
BASE = MS / "finalpaper_qe_mechanism_v41_validation.tex"

# (first line, last line) of each condensed block in the v41 base file, 1-indexed, inclusive
BLOCKS = {
    "intro": (159, 169), "elastodyn": (175, 239), "dft": (243, 251), "md": (256, 258),
    "stab1": (291, 296), "stab2": (306, 306), "elastic": (318, 320), "sep0": (324, 324), "sep": (334, 336),
}
FULL_TITLES = {"intro": "1. Introduction and Background", "elastodyn": "2.1 Elastodynamic Methods",
               "dft": "2.2 DFT Methods", "md": "2.3 Atomistic simulations of decomposed Fe--Cr",
               "stab1": "3.1 Spinodal stability from demixing free-energy curvature", "stab2": None,
               "elastic": "3.2 Elastic stiffening and elastodynamic mapping to Rayleigh SAW observables",
               "sep0": "3.3 Separating composition trends from spinodal conversion trends", "sep": None}

NEW = {}
NEW["intro"] = r"""Spinodal decomposition stiffens the \textit{bcc} $\delta$--ferrite of cast austenitic stainless steels (CASS), and this stiffening can be detected through picosecond ultrasonic surface acoustic wave (SAW) peak splitting \cite{al2023detecting, al2020detecting,al2020non}. This paper predicts that SAW signature from first principles, combining elastodynamics with density functional theory (DFT) to connect composition, elastic response, and Rayleigh-wave observables.

For high--melting-point metals such as $\delta$--ferrite, elasticity is closely tied to plastic response: dislocation energetics are largely stored in the Hookean lattice outside the dislocation core, and thermally activated processes during deformation can modify elastic properties \cite{pugh1954xcii}. Lattice elasticity also correlates with melting point through the shared dependence on interatomic bonding and electronic structure \cite{chang1966correlation,fine1984elastic, korzhavyi2009electronic, razumovskiy2011first, hosemann2018materials}. In ferromagnetic Fe--Cr these quantities vary strongly with concentration \cite{korzhavyi2009electronic}, and the observed peak splitting is consistent with an anomalous concentration dependence of the elastic constants associated with electronic-structure effects, including Fermi-surface topology changes \cite{razumovskiy2011first}.

During decomposition, chromium partitions into Fe-rich $\alpha$ and Cr-rich $\alpha'$ regions of differing stiffness, so elastic stiffening can serve as a mechanistically grounded proxy for embrittlement in stainless steels undergoing thermal and/or irradiation-assisted aging \cite{du2017facile, sdajani2019cassspinodal,al2020non, hofmann2019transient, al2020detecting,byun2016thermal,lach2020precipitation}. Because the Rayleigh SAW frequency depends directly on the elastic moduli, picosecond ultrasonic pump--probe spectroscopy can track this progression non-destructively \cite{hofmann2019transient}.

We therefore develop a thermodynamics-to-observable workflow for \textit{bcc} Fe$_{1-x}$Cr$_x$ as a model for $\delta$--ferrite: (i) spinodally unstable composition windows from the curvature of the demixing free energy, (ii) the composition dependence of the cubic elastic constants from DFT, and (iii) their mapping onto Rayleigh SAW frequencies. We then separate a solid-solution trend $f_{\mathrm{SAW}}(x_{\mathrm{Cr}})$ from an evolution trend $f_{\mathrm{SAW}}(\phi)$ in the spinodal conversion fraction $\phi$ of atomistic simulations of decomposed Fe--Cr \cite{al2023detecting}, so that SAW peak shifts and splitting can be read in terms of composition and decomposition state. The full text of Sections 1--3, before condensation to the journal's length, is given in SM Section S11."""

NEW["elastodyn"] = r"""Elastic waves in each homogeneous alloy obey the elastodynamic wave equation \cite{du2017facile}:
\begin{equation}
    \rho \frac{\partial^2 u_i}{\partial t^2} = C_{ijkl} \frac{\partial^2 u_k}{\partial x_j \partial x_l}
    \label{eqn:ewe}
\end{equation}
where $C_{ijkl}$ is the stiffness tensor, $\rho$ the density, and $u$ the displacement field. The SAW measured by transient grating spectroscopy (TGS) is the Rayleigh wave of a traction-free half-space. Because crystal orientation is not resolved in this step, we solve the surface problem for the isotropic Voigt--Reuss--Hill (VRH) average of each tensor (Sec.~\ref{sec:analysis}): $v_R = v_s\sqrt{q}$, where $q$ is the root in $0<q<1$ of
\begin{equation}
q^3 - 8q^2 + (24-16\beta)q - 16(1-\beta) = 0, \qquad \beta = v_s^2/v_p^2,
\label{eqn:rayleigh}
\end{equation}
with $v_s=\sqrt{G/\rho}$ and $v_p=\sqrt{(B+4G/3)/\rho}$. The SAW frequency is $f_{\mathrm{SAW}} = v_R/\lambda$ with $\lambda = 6.4\,\mu$m, the grating period of our TGS measurements \cite{al2023detecting}; density is interpolated linearly between \textit{bcc} Fe (7870 kg m$^{-3}$) and Cr (7200 kg m$^{-3}$). Orientation-resolved (pseudo-)SAW and bulk-wave frequencies are computed for 10{,}000 random grain orientations with the surface Green's-function method of Du and Zhao \cite{du2017facile}. The boundary conditions, the closed-form root used, and the inversion of stress--strain data for the elastic tensor are given in SM Section S11."""

NEW["dft"] = r"""Total energies, mixing enthalpies and elastic tensors were calculated with the Vienna \textit{ab initio} Simulation Package (VASP) \cite{kresse1996vasp} on 40-atom special quasirandom structure (SQS) supercells of \textit{bcc} Fe--Cr \cite{VandeWalleLiu2013}, generated with \texttt{mcsqs} (ATAT), with spin polarization, the projector-augmented wave method \cite{kresse1999paw}, the PBE functional \cite{pbe1996}, a 500 eV cutoff and $\Gamma$-point sampling, at $x_{\mathrm{Cr}} = 0$, 0.20, 0.50, 0.95, and 1.00. An independent Quantum ESPRESSO (QE) set on 40-atom SQS at twelve compositions, at the PBE-relaxed and the experimental lattice parameters, is described with its convergence and validation in SM Sections S3--S8; the self-consistent field procedure is outlined in SM Section S11."""

NEW["md"] = r"""The spinodal conversion fraction $\phi$ and the corresponding moduli come from our earlier LAMMPS study \cite{al2023detecting,plimpton1995}, which includes no electronic structure: a $50^3$ \textit{bcc} cell of 250,000 atoms (Fe$_{50}$Cr$_{50}$) was phase-separated by Monte Carlo swaps driven by a repulsive nearest-neighbour Fe--Cr interaction, intermediate states were produced by converting a fraction of the atoms back, and the elastic constants of each relaxed state were calculated by finite deformation with the Fe--Cr two-band embedded-atom-method (EAM) potential of Bonny \textit{et al.} \cite{bonny2011iron}, for ten random seeds. $\phi=0.5$ is the random alloy and $\phi=1$ the fully separated structure; the density is held at 7535 kg m$^{-3}$. With the same potential and stiffness procedure we computed random alloys across composition and repeated the decomposition at twelve compositions, $x_{\mathrm{Cr}}=0.05$--0.95, in $25^3$ cells (31,250 atoms; five alloys each, with a composition-preserving conversion), which reproduce the original Fe$_{50}$Cr$_{50}$ series (SM Fig.~\ref{fig:S_md_validation}). The full procedure is given in SM Section S11; all LAMMPS inputs, random seeds and potential files are in the repository."""

NEW["stab1"] = r"""The local stability of a binary solid solution is governed by the curvature of its free energy of mixing: where $\partial^2 \Delta G_{\mathrm{mix}}/\partial x_{\mathrm{Cr}}^2 < 0$, infinitesimal composition fluctuations lower the free energy and demixing proceeds without a nucleation barrier. Figure~\ref{fig:free_energy} is evaluated at 298.15 K from the concentration-dependent interaction models $\Omega(x)=(x-a)(bx^2+cx+d)$ that Levesque \textit{et al.} fitted to PWSCF--PAW and SIESTA--NC calculations \cite{levesque2011simple}, with the ideal configurational entropy $S_{\mathrm{mix}}$. Following their sign convention, we plot the demixing free energy
\begin{equation}
\Delta G_{\mathrm{demix}} = \Omega(x)\,x(1-x) + T S_{\mathrm{mix}} = -\Delta G_{\mathrm{mix}},
\label{eqn:demix}
\end{equation}
so the spinodal criterion appears as \emph{positive} curvature of $\Delta G_{\mathrm{demix}}$ in Fig.~\ref{fig:free_energy}b. The exact second derivative of the PWSCF--PAW model and its three-point numerical second differences coincide, and a cruder two-point scheme agrees on where the instability begins (SM Section S11)."""

NEW["stab2"] = r"""The homogeneous alloy is unstable for $0.21 \lesssim x_{\mathrm{Cr}} \lesssim 0.75$ in the PWSCF--PAW model (with a second, narrow window at $0.87$--$0.98$) and for $0.25 \lesssim x_{\mathrm{Cr}} \lesssim 0.99$ in the SIESTA--NC model. Both windows contain the Fe--36Cr alloy in which spinodal decomposition was recently detected through TGS stiffening \cite{rae2026tgs}. The entropic term narrows these windows at aging temperatures: at 748 K the CALPHAD assessment of Xiong \textit{et al.}, validated by heat capacity, Curie temperatures and atom-probe tomography, places the chemical spinodal at $x_{\mathrm{Cr}}\approx0.30$--0.71 within a miscibility gap of $\approx0.14$--0.88 \cite{xiong2011improved}. Fe--36Cr therefore lies inside the chemical spinodal at the aging temperature, close to its Fe-rich edge, consistent with their observation at 773 K that 37.8 at.\% Cr decomposes spinodally whereas 26.7 at.\% Cr decomposes by nucleation \cite{xiong2011improved}. Within these windows, elastic and acoustic properties must be read as those of a system driven toward compositional modulation, not of a homogeneous solid solution."""

NEW["elastic"] = r"""We compute the elastic tensor across $x_{\mathrm{Cr}}$ and propagate it through the surface problem of Sec.~2.1, using the longitudinal stiffness $C_{11}$ as a compact index of stiffening (Fig.~\ref{fig:elastic_elastodynamics}a). Fully decomposed MD structures are stiffer than the random alloy at every composition, most at intermediate compositions (257.4 versus 247.3 GPa at $x_{\mathrm{Cr}}=0.36$). The VASP tensors are systematically stiffer than the reference data ($C_{11}=324$ versus 260 GPa for \textit{bcc} Fe), possibly because of the finite SQS size and the magnetic configuration, so we interpret relative trends rather than absolute frequencies. The QE tensors at the PBE-relaxed lattice parameter match the literature at 20 and 95\% Cr ($f_{\mathrm{SAW}}$ within 3 MHz) and, at the experimental lattice parameter, follow it between 30 and 50\% Cr ($C_{11}$ within 11 GPa), levelling off where the literature $f_{\mathrm{SAW}}$ has its minimum.

Each homogeneous composition has a single isotropic $f_{\mathrm{SAW}}$ but a distribution of (P)SAW frequencies over grain orientations. Our earlier study defined an ``unallowed'' range for pristine CASS, above its highest (pseudo-)SAW frequency and below its bulk-wave frequency, and attributed (P)SAW peaks observed there after aging to stiffened $\delta$-ferrite \cite{al2023detecting}. Repeating that calculation for every Fe--Cr tensor (Fig.~\ref{fig:elastic_elastodynamics}b), all bulk waves lie above the range ($\geq0.87$ GHz), while in the (P)SAW branch no literature tensor up to $x_{\mathrm{Cr}}=0.7$ and no MD tensor reaches it, even after full decomposition (highest MD mode 585 to 596 MHz at $x_{\mathrm{Cr}}=0.36$). Only the stiffest tensors enter it: Cr-rich literature alloys (9\% and 29\% of grains at $x_{\mathrm{Cr}}=0.9$ and 1.0) and the present DFT tensors (30\% and 46\% at $x_{\mathrm{Cr}}=0.95$ and 1.0). Frequencies in this range therefore indicate stiffness beyond that of a homogeneous $\delta$-ferrite or Fe--36Cr solid solution, consistent with the stiffer Cr-rich $\alpha'$ formed by decomposition (Fig.~\ref{fig:fsaw_combined}a) and with the measured peak splitting (60--100\% of aged spots against 0--25\% unaged; Table~3 of Ref.~\cite{al2023detecting}).

Decomposition stiffens mainly the tetragonal shear modulus $C'$ (by 17--29\% for Fe--36Cr, through $\alpha'$) rather than $C_{44}$. It therefore cannot lift the top of the (P)SAW branch above about 585~MHz, but it lowers the elastic anisotropy and makes weak high-velocity pseudo-SAWs above 665~MHz, beyond any surface mode of pristine CASS, three to five times more frequent (QE: 0.4 to 1.9\% of grains; EMTO-CPA: 2.6 to 6.2\%), as in the measured tail (0.2\% of peaks unaged against 2--15\% after aging; Fig.~9 of Ref.~\cite{al2023detecting}; SM Section S9). Electronically, the stiffer $C'$ of $\alpha'$ follows the filling of the \textit{bcc} $d$ band \cite{soderlind1993theory} (QE: $C'$ from 56~GPa at 5\%~Cr to 145~GPa at 95\%~Cr), the electronic topological transition at $x_c\approx0.08$ sets the anomalous moduli of $\alpha$ (SM Section S6), and the spin- and charge-density waves of $\alpha'$ couple to the volume rather than to shear, lowering $K/G$ without affecting the acoustic signature (SM Section S9).

The same shear stiffening links the acoustic signal to embrittlement. Decomposition raises $G$ (+11\% for MD Fe--36Cr) at a nearly unchanged $K$ ($-1.2\%$), so Pugh's ratio $K/G$ falls by 11\% in step with the SAW shift, $\Delta(K/G)/(K/G)\approx-2\,\Delta f/f$ \cite{pugh1954xcii}; the $\alpha$/$\alpha'$ shear-modulus misfit (+16 to +23~GPa, QE and EMTO-CPA) hardens the ferrite against dislocation motion \cite{kato1981hardening}, while the cleavage resistance, which scales with $K$, stays fixed. Stiffening therefore shifts the balance from plastic flow towards cleavage, as in the Charpy energies of CF8, which fall from 88.5~J to 14--24~J after 10{,}000~h at 330--400~$^\circ$C as the SAW splitting grows from 12 to 128--229~MHz \cite{al2023detecting}. Fe--Cr, melting above 1500~$^\circ$C, lies in the athermal regime where Pugh's relations hold (SM Section S11)."""

NEW["sep0"] = r"""Composition and spinodal conversion both shift $f_{\mathrm{SAW}}$, so we report the two relationships separately (Fig.~\ref{fig:fsaw_combined})."""

NEW["sep"] = r"""Figure~\ref{fig:fsaw_combined}a shows how chromium partitioning produces stiffening. Fe--36Cr lies near a local minimum of the homogeneous $f_{\mathrm{SAW}}(x_{\mathrm{Cr}})$ curve, and decomposition moves material toward Fe-rich $\alpha$ and Cr-rich $\alpha'$ compositions (arrows), in both directions to higher frequency. The five-composition VASP set is consistent with this non-monotonic dependence but cannot resolve the minimum. The MD random alloys, an independent classical model, also have their minimum at $x_{\mathrm{Cr}}\approx0.36$--0.40 although their $C_{11}$ rises monotonically: the minimum comes from the shear moduli, which control the Rayleigh velocity ($C_{44}$ falls from 116 GPa for Fe to 106 GPa at $x_{\mathrm{Cr}}=0.5$). Figure~\ref{fig:fsaw_combined}a describes the constituents; the response of the two-phase microstructure also depends on phase fractions and on the domain size relative to the SAW wavelength.

Figure~\ref{fig:fsaw_combined}b isolates decomposition at fixed overall composition. For MD Fe$_{50}$Cr$_{50}$, $f_{\mathrm{SAW}}$ rises monotonically from 0.469 GHz for the random alloy ($\phi=0.5$) to 0.490 GHz for the fully separated structure ($\phi=1$), about 4.6\% at fixed density, the same sign as the first-principles results. Across the twelve compositions the rise is nearly symmetric about 50/50, from 0.5\% at $x_{\mathrm{Cr}}=0.05$ and 0.4\% at 0.95 to 5.0\% at Fe--36Cr and 5.7\% at 50/50 (every composition in SM Section S11; standard deviations at most 0.1 percentage points): stiffening on decomposition is largest within the spinodally unstable window of Fig.~\ref{fig:free_energy}b and nearly vanishes at the dilute ends. The fully separated state is an idealized limit, so these values bound the stiffening from above. Together, $f_{\mathrm{SAW}}(x_{\mathrm{Cr}})$ gives the baseline response of a homogeneous alloy and $f_{\mathrm{SAW}}(\phi)$ the additional stiffening of progressive spinodal conversion, so SAW measurements can distinguish early-stage spinodal evolution from equilibrium compositional effects."""


def unlabel(txt):
    """Copy of main-text paragraphs for the SM: no labels (they stay in the main text), unnumbered equations."""
    txt = re.sub(r"\s*\\label\{[^}]*\}", "", txt)
    txt = txt.replace(r"\begin{equation}", r"\begin{equation*}").replace(r"\end{equation}", r"\end{equation*}")
    return txt


def main():
    base = BASE.read_text().split("\n")
    old = {k: "\n".join(base[a - 1:b]) for k, (a, b) in BLOCKS.items()}
    full = [r"\subsection*{S11. Full text of Sections 1--3}",
            "The main text condenses Sections 1--3 to the journal's word limit. Their full text before condensation "
            "follows; equation, figure and section numbers refer to the main text."]
    for k in BLOCKS:
        if FULL_TITLES[k]:
            full.append(rf"\paragraph{{{FULL_TITLES[k]}.}}")
        full.append(unlabel(old[k]))
    sm_add = "\n\n".join(full) + "\n"
    for v in ("", "_cubes", "_trend", "_trend_cubes"):
        src = MS / f"finalpaper_qe_mechanism_v41_validation{v}.tex"
        dst = MS / f"finalpaper_qe_mechanism_v42_validation{v}.tex"
        assert not dst.exists(), dst
        t = src.read_text()
        for k in BLOCKS:
            assert t.count(old[k]) == 1, (v, k)
            t = t.replace(old[k], NEW[k])
        i = t.rindex(r"\end{document}")
        t = t[:i] + "\n" + sm_add + "\n" + t[i:]
        dst.write_text(t)
        print("wrote", dst.name)


if __name__ == "__main__":
    main()
