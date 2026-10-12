# elastopya — Fe–Cr spinodal stiffening and surface acoustic waves

Code, inputs, data and figures for

> S. A. Al Dajani and M. P. Short, *First principles evidence of stiffening during spinodal
> decomposition in FeCr* (in preparation, Scripta Materialia).

Everything that produces a number or a figure in the paper is here, together with the inputs
needed to regenerate it. Large regenerable files (LAMMPS restart files, coordinate dumps, Quantum
ESPRESSO scratch/wavefunctions) are not stored.

## Layout

| Path | Contents |
|---|---|
| `manuscript/` | LaTeX source, bibliography, figures and compiled PDF of the manuscript |
| `scripts/make_figures.py` | Regenerates every figure from `data/` (`python scripts/make_figures.py`) |
| `data/sheet_*.csv` | Bounded snapshots of the working spreadsheet used for Figs. 1–3 (literature tensors, author DFT tensors, Acta MD moduli, thermodynamic model, band edges) |
| `scripts/saw_du.py` | Python port of the Du & Zhao (2017) surface Green's-function (P)SAW solver |
| `scripts/saw_random_orientations.py`, `scripts/run_saw_*.sh` | Orientation-resolved (P)SAW spectra over 10,000 random grains |
| `data/saw/` | Grain orientations, (P)SAW speeds and intensities for every tensor; validation run for SS304 |
| `lammps/Spinodal_Generation_Care_Package/` | Original LAMMPS workflow of L. Waldron used for Al Dajani *et al.*, Acta Mater. 246 (2023) 118552 (Sec. 3.7) |
| `lammps/random_alloys/` | Elastic constants of random bcc Fe–Cr alloys (same stiffness procedure), inputs and logs |
| `lammps/decomposition/` | Spinodal-decomposition series at several compositions (25³ cells); inputs, logs, gzipped structures |
| `dft_qe/` | Quantum ESPRESSO elastic-constant drivers, SQS cells, inputs/outputs, convergence study |
| `visualizations/` | Snapshots of the decomposed structures and local-composition histograms |

## Reproducing the figures

```bash
pip install numpy scipy matplotlib openpyxl
python scripts/collect_md_random.py          # MD random-alloy constants -> data/md_random_alloys.csv
python scripts/collect_md_decomposition.py   # MD decomposition series   -> data/md_decomposition.csv
python scripts/phase_compositions.py         # alpha / alpha' compositions of decomposed structures
python scripts/visualize_decomposition.py    # structure snapshots
python scripts/make_figures.py               # all figures in manuscript/figures/
```

## Methods and validation

**Rayleigh/SAW frequencies (Figs. 2–3).** Isotropic estimates use the Voigt–Reuss–Hill moduli
of each cubic tensor. Orientation-resolved spectra use `saw_du.py`, a line-by-line port of the
MATLAB `getSAW.m` of Du & Zhao (npj Comput. Mater. 3, 17 (2017)); the MATLAB code itself is not
redistributed here. The port was validated against the original MATLAB run for 304 stainless steel
(10,000 orientations, `data/saw/SS304_speeds_10k.csv`): SAW speeds agree exactly, and with the
original PSAW selection (intensity ≥ 2×10⁻³ of the SAW) the maximum frequency is 609.0 MHz at a
6.4 µm grating, as reported in Acta Mater. 246 (2023) 118552.

**MD (LAMMPS).** Bonny *et al.* (2011) two-band Fe–Cr EAM potential, NIST Interatomic Potentials
Repository release `ipr3` (numerically identical to the original files apart from format fixes
required by current LAMMPS). Random alloys: 50×50×50 bcc cells, ten random seeds at the DFT
compositions and three elsewhere; the Fe₅₀Cr₅₀ value reproduces Acta (C₁₁ = 252.9 vs 252.7 GPa).
Decomposition: Waldron's Monte Carlo + conversion procedure in 25³ cells with the same swap
attempts per atom, with a composition-preserving conversion away from 50/50; validated against the
Acta 50³ series (C₁₁ within 0.8 GPa at every conversion fraction).

**DFT (Quantum ESPRESSO 7.5).** PBE-PAW (pslibrary 0.2.1 Fe / 1.0.0 Cr, downloaded from the QE
pseudopotential library; not stored here). Convergence study for bcc Fe in
`data/qe_fe_convergence.csv`. 40-atom SQS cells generated with `sqsgenerator`.

## Software

LAMMPS 22 Jul 2025 (conda-forge, Open MPI), Quantum ESPRESSO 7.5 (built from source), Python 3.12
with NumPy/SciPy/Matplotlib, sqsgenerator 0.5.6, ASE 3.29.

## License

MIT (see `LICENSE`). Third-party inputs (LAMMPS potential files, the Waldron workflow) remain
under their original terms and should be cited as indicated above.
