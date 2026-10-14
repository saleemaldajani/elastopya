#!/bin/bash
#SBATCH -p mit_normal
#SBATCH -A mit_general
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --mem=120G
#SBATCH -t 08:00:00
#SBATCH -o slurm-%x-%j.out
# Convergence of the Fe (FM) DFPT frequencies, which come out ~15% above the neutron data: same cell as
# phonon/Fe_fm_aexp, single q-points (0.5,0,0) [mid Gamma-H] and N (0.5,0.5,0), for ecut 60/480 vs 90/1080 Ry and
# k 16^3 vs 24^3. Run from /orcd/pool/004/sdajani/fecr_qe.
cd /orcd/pool/004/sdajani/fecr_qe
mkdir -p phonon_test/fe_conv
source qe75_env.sh
BIN=$(dirname $QE_PW)
SRC=/orcd/pool/004/sdajani/fecr_qe/phonon/Fe_fm_aexp/scf.in
PREFIX=$(grep -o "prefix *= *'[^']*'" $SRC | sed "s/.*= *//")
cd phonon_test/fe_conv
for S in e60_k16 e90_k16 e60_k24; do
  mkdir -p $S; cd $S
  case $S in e60_k16) EW=60; ER=480; K=16;; e90_k16) EW=90; ER=1080; K=16;; e60_k24) EW=60; ER=480; K=24;; esac
  sed -e "s/ecutwfc *= *[0-9.]*/ecutwfc = $EW/" -e "s/ecutrho *= *[0-9.]*/ecutrho = $ER/" \
      -e "/K_POINTS/{n;s/.*/$K $K $K 0 0 0/}" -e "s#outdir *= *'[^']*'#outdir = './tmp'#" $SRC > scf.in
  $QE_MPIRUN -np 32 $QE_PW -nk 8 -in scf.in > scf.out
  for Q in "0.5 0.0 0.0" "0.5 0.5 0.0"; do
    n=$(echo $Q | tr ' .' '__')
    cat > ph_$n.in <<EOF2
Fe q $Q
&inputph
  prefix = $PREFIX, outdir = './tmp', fildyn = 'dyn_$n', tr2_ph = 1e-15, amass(1) = 55.845,
  alpha_mix(1) = 0.2, nmix_ph = 12
/
$Q
EOF2
    $QE_MPIRUN -np 32 $BIN/ph.x -nk 8 -in ph_$n.in > ph_$n.out
    echo "$S q=($Q): $(grep -E 'freq \(' ph_$n.out | awk '{print $5}' | tr '\n' ' ') THz"
  done
  cd ..
done
echo TASK_DONE
