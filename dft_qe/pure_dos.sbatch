#!/bin/bash
#SBATCH -p mit_normal
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=120G
#SBATCH -t 06:00:00
#SBATCH -o slurm-%x-%j.out
# Total DOS of the pure metals for validation against the literature: scf as in fermi/<CASE>/scf.in (same cell,
# cutoffs, smearing, magnetic start), then a tetrahedron (tetrahedra_opt) nscf on a dense mesh and dos.x.
cd /orcd/pool/004/sdajani/fecr_qe
source qe75_env.sh
BIN=$(dirname $QE_PW)
for CASE in Fe_fm Fe_fm_aexp Cr_nm Cr_nm_aexp Cr_afm Cr_afm_aexp; do
  d=pure_dos/$CASE; mkdir -p $d; cd $d
  [ -f dos.dat ] && { cd - >/dev/null; continue; }
  cp ../../fermi/$CASE/scf.in scf.in
  case $CASE in Cr_afm*) K=24; NB=40;; *) K=36; NB=20;; esac
  sed -e "s/calculation = .scf./calculation = 'nscf'/" -e "s/occupations = .smearing., smearing = .mv., degauss = 0.02,/occupations = 'tetrahedra_opt',/" \
      -e "s/nbnd = [0-9]*/nbnd = $NB/" -e "/K_POINTS/{n;s/.*/$K $K $K 0 0 0/}" \
      -e "s/conv_thr = 1e-9, mixing_beta = 0.3/conv_thr = 1e-9, mixing_beta = 0.3, diagonalization = 'cg'/" scf.in > nscf.in
  # CG: Davidson aborted in cdiaghg (MPI errors 74/76) in the Fe_fm and Cr_nm_aexp nscf runs (2026-10-15)
  grep -q "nbnd" nscf.in || sed -i "s/nspin = \([0-9]\)/nspin = \1, nbnd = $NB/" nscf.in
  printf "&dos\n  prefix = 'fs', outdir = './tmp', fildos = 'dos.dat', Emin = -12, Emax = 30, DeltaE = 0.01\n/\n" > dos.in
  $QE_MPIRUN -np 32 $QE_PW -nk 8 -in scf.in > scf.out && \
  $QE_MPIRUN -np 32 $QE_PW -nk 16 -in nscf.in > nscf.out && \
  $QE_MPIRUN -np 1 $BIN/dos.x -in dos.in > dos.out
  grep -q "convergence NOT" scf.out nscf.out && echo "$CASE: SCF not converged"
  rm -rf tmp/*.save/wfc*.dat
  cd - >/dev/null
done
echo TASK_DONE
