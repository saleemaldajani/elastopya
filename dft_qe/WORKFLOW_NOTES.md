# First-principles workflow: final scripts and what had to be fixed

These notes cover the Quantum ESPRESSO (QE) calculations behind the revised Scripta manuscript and its supplement. Part 1 lists the final scripts. Part 2 lists the problems found along the way and how each was fixed, so the history behind every option is on record. The scripts are the final versions.

All calculations ran on MIT Engaging (ORCD) in `/orcd/pool/004/sdajani/fecr_qe`. `pull_engaging.sh` mirrors every input, output, log and job script into `dft_qe/engaging/`, without wavefunction scratch or pseudopotentials.

## 1. Final scripts

### Elastic constants
| Script | Purpose |
|---|---|
| `elastic_sqs_qe.py` | 40-atom SQS: relax, then ±0.5% strains (xx, yz) with ionic relaxation; C11, C12, C44 from central differences. Environment options are listed in §2. |
| `elastic_qe.py` | Pure bcc Fe (FM) and Cr (AFM or NM), 2-atom cell; same procedure. |
| `engaging/sqs40.sbatch`, `engaging/sqs40_variant.sbatch` | SQS jobs. The variant version writes `runs/sqs40_x<X>_k442_e60_<TAG>` for alternative settings (experimental lattice, magnetic start, mixing). |
| `engaging/bcc2.sbatch`, `bcc2_variant.sbatch` | Pure-element jobs; the variant writes `runs/<El>_<mag>_k<K>_e60_<SUFFIX>`. |

Settings (all runs):
- PBE with pslibrary PAW (Fe 3s3p3d4s, Cr 3s3p3d4s), 60/480 Ry, Marzari–Vanderbilt smearing 0.02 Ry.
- SQS: 2×2×5 conventional cells (40 atoms) with a 4×4×2 k-mesh. Pure-element references for mixing energies use 8×8×10, the same k-point density.

### Electronic structure
| Script | Purpose |
|---|---|
| `fermi_slices.py` | Spin-resolved Fermi-surface sections of pure Fe and Cr (Korzhavyi 2009, Fig. 6 layout). |
| `bands_path.py` | Band structures along Γ–H–N–Γ–P–H (pure elements), and unfolded SQS bands (BandUPpy). |
| `unfold_fs.py` | Unfolded (Bloch spectral function) Fermi surfaces of Fe-rich SQS. |
| `dos_sqs.py`, `dos_plot.py`, `plot_dos.py` | Spin- and atom-projected DOS of the SQS (projwfc.x). |
| `nesting.py`, `sdw_cr.py`, `density_maps.py` | Cr nesting function, SDW/CDW supercell, spin- and deformation-charge maps. |
| `pure_dos.sbatch` | Total DOS of the pure metals (tetrahedra, dense mesh), for comparison with the literature. |

### Lattice dynamics
| Script | Purpose |
|---|---|
| `phonons.py` | DFPT (ph.x, q2r.x, matdyn.x) dispersions and DOS for Fe FM, Cr NM and Cr AFM. Its plot step refuses imaginary modes in phases expected to be stable. |
| `phonopy_cubic.py` + `phcubic.sbatch` | Finite-displacement Fe phonons in a cubic 3×3×3 supercell. **This replaces Fe DFPT; see §2.** |
| `phonopy_bcc_bands.py` | Dispersion and DOS from phonopy force constants, mapped to the primitive bcc cell. |
| `phonopy_sqs.py` | SQS phonon DOS and vibrational entropy of mixing (finite displacements). |
| `phonon_ref_plot.py`, `phonon_dos_plot.py` | Phonons against neutron, NRIXS and Fultz 1995 data. |
| `fe_phconv.sbatch`, `engaging/phonon_test/*`, `engaging/frozen_Fe/*` | Convergence and frozen-phonon checks that diagnosed the Fe DFPT problem. |

### Validation and figures
| Script | Purpose |
|---|---|
| `validation_overlays.py` | Every variant overlaid on literature and experiment, with deviations written to `validation_summary.csv`. |
| `refdata/` | Reference data, tabulated or digitized. Each CSV header states the source, figure, calibration and uncertainty. |
| `refdata/digitize_scripts/` | Scripts used to digitize the raster figures of Korzhavyi 2009 and Klaver 2006: tick calibration, colour/ridge tracing, marker finding, overlay checks. `*_final.py` / `*_out.py` produced the CSVs; the rest are the iterations that led to them. |
| `../scripts/make_figures*.py`, `make_si_figures.py`, `make_fig4_mechanism.py`, `make_fig5_embrittlement.py`, `make_qe_aexp_progress.py`, `visualize_*.py` | Paper and supplementary figures (K-Dense Scripta style). |

### Cluster tooling
| Script | Purpose |
|---|---|
| `pull_engaging.sh` | Mirror Engaging into `engaging/`, including every `*.py`, `*.sbatch` and `*.sh`, subfolders included. |
| `watch_connection.sh` | Keep the shared `ssh engaging` connection alive (ping every 2 min) and the Mac awake (`caffeinate`); one notification if the connection is lost. Logging in needs Duo, so it cannot reconnect by itself. |
| `wait_key.sh` | Background waiter: wakes the assistant when key results appear, a job fails, or the pool is more than 90% full. |
| `watch_engaging.sh`, `keepalive_engaging.sh`, `local_guardian.sh` | Earlier versions of the same tooling, kept for the record. |

## 2. Problems found and how they were fixed

**QE version.**
- The Engaging QE 7.3.1 module gave Fe C44 = 118.7 GPa, against 114.1 GPa from our local QE 7.5.
- Fix: QE 7.5 was built on Engaging with the same options (`build_qe75*.sbatch`). It reproduces the local results to within 0.1 GPa (`validate_fe_qe75.sbatch`).

**Plane-wave cutoff.**
- At 40 Ry the vc-relax collapses the cell (−12% volume in a 16-atom test).
- Fix: 60/480 Ry. `elastic_sqs_qe.py` stops if the relaxed volume leaves 10.8–12.2 Å³/atom.

**False "JOB DONE".**
- QE prints JOB DONE even when "convergence NOT achieved".
- Fix: the drivers check for this and stop instead of building on an unconverged step.
- Finished steps are reused on restart, so a requeued or resubmitted job resumes where it stopped.

**SCF convergence of difficult compositions** (Fe–80Cr, the shear strains at 65–80% Cr). All of these are environment options of `elastic_sqs_qe.py`; with none set, the input is unchanged.
- `QE_MIXING_BETA`, `QE_MIXING_NDIM`, `QE_MAXSTEP`: mixing and iteration count. Gentler mixing helps some cases and hurts others; at 80% Cr (experimental lattice) β = 0.05 stalled at 2.5×10⁻⁴ Ry while β = 0.1 reached 3×10⁻⁵ Ry.
- `QE_MAG_CR`: Cr starting moment. A negative value starts Cr antiparallel to Fe, as needed for Fe–80Cr.
- `QE_UPSCALE`: during relaxations QE tightens `conv_thr` by up to 100×, and Fe–80Cr then oscillated at 1–5×10⁻⁹ Ry for 600 iterations. `upscale = 1` keeps 10⁻⁸ Ry, which is ample for forces converged to 10⁻³ Ry/bohr.
- `QE_START_FROM_FILE`: restart from the saved potential and wavefunctions.
- Every failed attempt is kept, renamed `*_INVALID_<reason>` with a README of what was tried.

**Volume convention.**
- PBE relaxes every composition to about 11.4 Å³/atom, while the experimental volume is 11.8–12.0 Å³/atom.
- At the relaxed volume the bulk modulus is about 30% above experiment; at the experimental lattice parameter pure Fe gives K = 170 GPa against 168 measured.
- Fix: `QE_FIXED_A=exp` in both drivers. It fixes the cubic cell at the room-temperature experimental lattice parameter (Vegard between Fe 2.8665 Å and Cr 2.8846 Å; ions relaxed), the convention of Razumovskiy et al. (2011).
- Mixing energies change by at most 10 meV/atom between the two conventions.

**Diagonalisation.**
- Davidson failed in some band and NSCF steps (QE errors 227 and 1243, "S matrix not positive definite").
- Fix: those steps use CG.
- CG was then too slow for the 510-k-point unfolding band step: it timed out at 24 h on one node. Fix: `unfold2.sbatch` runs on two nodes with the 2-day limit. It first tries ParO and falls back to CG if ParO fails. PPCG was the first choice, but QE 7.5 rejects it ("PPCG diagonalization not supported anymore (Dec. 2024)").

**Projected DOS.**
- projwfc with Gaussian broadening on a tetrahedron NSCF gave all-zero PDOS. Fix: tetrahedra with `-nk 1` (`pdosT.*`).
- The energy window ended about 1.8 eV above E_F because of `nbnd`. Fix: `NBND=1.6` (wide window, `pdosW.*`).

**Unfolding.**
- BandUPpy returns `(bands, kpts)`, and the energy is the third band column; both were read wrongly at first.
- Unfolding was attempted after a failed band step. Fix: the job now checks first.
- The unfolded energies from BandUPpy are already relative to E_F. Subtracting E_F again in `bands_path.py sqs-plot` pushed every band out of the window (a blank plot). Fix: subtract only when the energies are absolute. `uf_fermi()` also maps the Engaging path in `scf_source.txt` to the local mirror.

**Phonons.**
- On Engaging, `q2r.x` and `matdyn.x` must run under `mpirun -np 1` (MPI_Init error otherwise).
- Non-magnetic Cr showed imaginary modes, −4.4 THz at k = 16³. Neither SCF mixing nor a 0.4% change of the lattice parameter changed them; k = 24³ or 0.04 Ry smearing removed them. They came from under-sampling the nested Fermi surface, not from an instability.
- **Fe DFPT is about 60% too stiff.** At the transverse mode q = (0, 0, ⅓):

  | Method | Frequency |
  |---|---|
  | Neutron scattering | 4.31 THz |
  | Frozen phonon, identical settings | 4.40 THz |
  | DFPT | 6.99 THz |

  The PAW dataset, cutoff and magnetism are therefore fine; ph.x is what fails for ferromagnetic Fe here. Cr DFPT, including AFM Cr, agrees with neutron data to 0.35–0.41 THz.
- Fix: Fe phonons come from finite displacements.
- The first finite-displacement reference used the 2×2×5 SQS-shaped cell. Its tetragonal image pattern made the symmetrised force constants inconsistent (on-site constants 10.0 vs 11.9 eV/Å² along x and z) and gave a spurious imaginary branch. Fix: a cubic 3×3×3 supercell (`phonopy_cubic.py`). Result at the experimental lattice parameter: RMS 0.23 THz against neutron scattering (DFPT: 1.1–1.4 THz), and a mean phonon frequency of 6.50 THz against 6.51 THz from NRIXS (DFPT: +19% to +27%).

**Mixing-energy references.** Pure Fe and Cr references for ΔH_mix use the same k-point density as the SQS: 8×8×10 for the 2-atom cell, against 4×4×2 for the 2×2×5 SQS. A 16³ reference would differ by a few meV/atom, comparable to the Fe-rich ΔH.

**Cluster operations.**
- *Preemption.* Jobs on `mit_preemptable` are requeued after preemption; `mit_normal` allows 12 h and 96 CPUs per user.
- *Queue priority.* After heavy use, priority, not the per-user CPU cap, is what holds jobs. Short jobs (the pure-element elastic constants take about 3 min) run at once on `mit_quicktest` (15 min, 48 CPUs).
- *Disk quota (2026-10-14).* The 1.1 TB pool filled with wavefunction scratch, and one DOS job died at 13.5 h with exit code 2.
  - Fix: scratch of finished runs was deleted, logged in `CLEANUP_LOG.txt` on Engaging. For elastic and phonon runs the whole `tmp/` went; for DOS runs only the wavefunctions, because the density maps read the charge densities. Usage went from 100% to 63%.
  - `wait_key.sh` now alarms above 90%.
  - *Second fill (2026-10-15).* Each unfolding band step writes about 290 GB of wavefunctions (510 k × 2 spins × 392 bands), and the unfolding reads them back. With the faster ParO solver, three ran at once and filled the pool within an hour. Fix: the unfolding jobs run one at a time (`--dependency=afterany` chain). `unfold2.sbatch` and `bandpath_sqs.sbatch` delete their wavefunctions once the unfolded spectra (`unfolded_*.npz`, `ebs_*.npz`) exist, keeping the charge density the band-structure step reuses. A manual cleanup must check that the band step did *not* finish (`JOB DONE` in `bands.out`) before deleting its wavefunctions. Two finished band steps were lost to this on 2026-10-15: an error check matched a harmless word, and the restarted jobs then skipped the band step and found no wavefunctions to unfold.
  - Jobs that ran while the disk was full can carry silently truncated wavefunction files. On this NFS pool a write may appear to succeed and fail later. Two jobs crashed only when reading them back (`davcio` error 10: "error writing" or "error reading file ./tmp/*.wfc*"). Fix: rename the failed output, delete that run's `*.wfc*`, keep the SCF charge density written before the disk filled, and rerun the step.
- *Connection.* Logins need Duo. The shared ControlMaster connection must be started detached (`ssh -f -N`), so closing Terminal does not end it. `watch_connection.sh` keeps it busy and the Mac awake.

## 3. Conventions
- **Never overwrite** a result or figure. New settings go in new directories (`_aexp`, `_crafm2`, `_v2`, …) and new figure versions (`_v2`, `fig5_embrittlement_v6`, …). Superseded files are kept in `v_before_*` or `v1_backup` folders.
- **Numbers in the manuscript** come from the files mirrored here (`engaging/runs/*/result.json`, `refdata/`, `validation_summary.csv`), read by the plotting scripts. None are typed in by hand.

## 4. Why the analysis is set up this way

- **Two DFT codes (VASP and QE).** VASP gave the original tensors. The VASP archive could not be recovered to extend or check it, so QE (open source, run on MIT Engaging) gave an independent, fully documented set at more compositions. Both are kept in the figures. A VASP cross-check of the QE tensors is planned once a licence is available.
- **40-atom SQS.** Large enough for pair correlations close to random over three neighbour shells and for composition steps of 2.5 at.%, small enough for stress–strain runs with ionic relaxation at many compositions.
- **Stress–strain at ±0.5% with ionic relaxation.** Central differences cancel the residual stress of the relaxed cell. Ionic relaxation is required in an SQS because internal strain lowers C44 and C′.
- **Two lattice-parameter conventions.** The main text keeps the PBE-relaxed series, the self-consistent choice for a DFT paper. The supplement adds the room-temperature experimental lattice parameter. That is the convention of the EMTO-CPA reference data, and it brings the bulk modulus within 1–7% of experiment against about 30% at the PBE volume (Fig. S17). Mixing free energies barely depend on the choice (≤ 10 meV/atom), so both appear in Fig. 1(a).
- **Validation before interpretation.** Every quantity in the supplement was replotted on published data it should reproduce, with the deviation stated in the caption (Section S3.1). Where a method failed that test (Fe DFPT), it was replaced, not explained away.
- **Finite displacements for Fe, DFPT for Cr.** DFPT reproduces neutron data for Cr (non-magnetic and AFM) but is about 60% too stiff for ferromagnetic Fe with these settings. Finite displacements in a cubic supercell reproduce the Fe neutron and NRIXS data. Each element therefore uses the method that passes validation.
- **Never overwrite.** Superseded results stay on disk (renamed `_INVALID_*`, `v1_backup`, `v_before_*`), so every number in an earlier draft can still be traced.

## 5. Changelog of script changes (with reasons)

| Date | Script | Change | Why |
|---|---|---|---|
| 2026-10-12 | `elastic_sqs_qe.py` | Mixing options; stop on "convergence NOT achieved" | Difficult SCFs; QE prints JOB DONE even when unconverged |
| 2026-10-12 | `phonons.py` | Gentle DFPT mixing; refuse imaginary modes in stable phases | NM Cr imaginary modes (later traced to k-point sampling) |
| 2026-10-12 | `dos_sqs.py` | Tetrahedra projwfc (`-nk 1`); `QE_NBND_FACTOR` | All-zero PDOS with Gaussian broadening; window cut at +1.8 eV |
| 2026-10-13 | `elastic_sqs_qe.py` | `QE_MAG_CR`; `QE_FIXED_A` (other session) | Fe–80Cr needs Cr antiparallel; experimental-lattice series |
| 2026-10-13 | `validation_overlays.py` | New | Overlay every variant on literature and experiment, with deviations |
| 2026-10-13 | `unfold2.sbatch` | Two nodes, 2-day limit, fast solver with CG fallback | Unfolding band step timed out at 24 h |
| 2026-10-14 | `unfold2.sbatch` | PPCG → ParO | PPCG was removed from QE 7.5 |
| 2026-10-14 | `elastic_sqs_qe.py` | `QE_UPSCALE` | Relaxation SCF stalled at a 10⁻⁹ Ry noise floor after QE tightened its threshold |
| 2026-10-14 | `elastic_qe.py`, `bcc2_variant.sbatch` | `QE_FIXED_A` for pure elements | End points and mixing references at the experimental lattice parameter |
| 2026-10-14 | `phonopy_cubic.py`, `phcubic.sbatch`, `phonopy_bcc_bands.py` | New | Fe phonons by finite displacements (DFPT too stiff; 2×2×5 cell inconsistent) |
| 2026-10-14 | `phonons.py`, `phonon_ref_plot.py --fe-fd` | Read finite-displacement Fe (`phonopy_derived/`) | Supplementary phonon figures use the validated Fe phonons (`S_phonons_v2`, `S_phonons_vs_exp_v2`) |
| 2026-10-14 | `validation_overlays.py` | `ph_dir()`; pure-element bulk moduli | Finite-displacement Fe; the measured end points test the volume convention |
| 2026-10-14 | `../scripts/make_figures.py` | `qe_mixing(aexp=True)`; experimental-lattice ΔG as diamonds in Fig. 1(a) | Complete the QE ΔG set; the two conventions agree to 10 meV/atom |
| 2026-10-14 | `../scripts/make_figures_kdense.py`, `../scripts/visualize_sqs.py` | Per-block widths in `add_strip`; tight unit-cell blocks | Larger DFT cells in Fig. 1(d) at one common scale |
| 2026-10-14 | `pull_engaging.sh` | Pull every script (subfolders included) and `phonopy_cubic/` | Some job scripts existed only on Engaging |
| 2026-10-14 | `watch_connection.sh`, `wait_key.sh` | New | Keep the Duo-authenticated connection alive; catch results, failures and a full disk |
| 2026-10-15 | `ett_analysis.py` | New: Lifshitz ETT analysis (SM S3.2, Fig. S27) | Locate the Fe-rich transition from our bands (critical majority level at H), test it against the EMTO-CPA moduli and measured bulk moduli with one shared onset, and split the lever-rule stiffening into α/α′ terms. A first version with free per-modulus onsets was discarded: poorly constrained, and it extrapolated the square-root term beyond its range of validity |
| 2026-10-15 | `bands_path.py`, `unfold_fs.py` | Panel titles state the lattice convention | Panels at different lattice parameters are shown side by side |
| 2026-10-15 | `thermal_check.py` | New: S1 thermal estimates for every tensor, validated (Fig. S_thermal, SM S8) | Debye temperature against calorimetry and phonon DOS; D_min against TGS in CF8. Found that the working sheet used h×1e11 in place of k_B and m = 2 (×7.6 on κ and D) |
| 2026-10-15 | `../scripts/reorganize_sm.py` | New: supplement rebuilt in topic order (v12) | Text first, then figures, pinned with \FloatBarrier; figures previously split unrelated text |
| 2026-10-15 | `../scripts/fix_thermal_panels.py` | New: corrected SM Fig. S7e,f (`DvsX_v5`, `KvsX_v3`) in the original Excel style | κ and D now follow Eq. (2) of S1 (k_B, m = 1), which divides the earlier values by 7.62. Inputs are the documented VASP tensors. The original 20% and 50% points (2–5% different) could not be traced to the current sheet |
| 2026-10-15 | `../lammps/decomposition/md_decomp_seeds.sbatch`, `../lammps/build_lammps_engaging.sbatch` | Seeds 2–5 at all twelve MD compositions (48 tasks) on Engaging; Figs 3–5, S_thermal and the text now use the five-alloy mean ± s.d. | One alloy per composition gave no error estimate. The cluster LAMMPS lacked the MC package (`fix atom/swap`), so the same LAMMPS version as the local runs was built with MANYBODY+MC; its seed-1 check at Fe–36Cr matches the local result within 0.14%. Seed-to-seed s.d. of the stiffening is ≤ 0.1 percentage points; Fe–36Cr changes from +5.1% (seed 1) to +5.0 ± 0.03% |
