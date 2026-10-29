# Quantum ESPRESSO 7.5 built in qe75/ with the same options as the local build
source /etc/profile >/dev/null 2>&1
module load gcc/12.2.0 openmpi/4.1.4 openblas/0.3.26 fftw/3.3.10
export QE_PW=/orcd/pool/004/sdajani/fecr_qe/qe75/q-e-qe-7.5/build/bin/pw.x
# Fallback (2026-10-21): on some nodes (node2701) the module load left mpirun off PATH and four jobs died at once
SPK=/orcd/software/core/001/spack/pkg
for d in $SPK/openmpi/4.1.4/zuyo6jx $SPK/openblas/0.3.26/ro5tivv $SPK/fftw/3.3.10/2qziucy; do
  case ":$LD_LIBRARY_PATH:" in *":$d/lib:"*) ;; *) export LD_LIBRARY_PATH=$d/lib:$LD_LIBRARY_PATH;; esac
done
command -v mpirun >/dev/null 2>&1 || export PATH=$SPK/openmpi/4.1.4/zuyo6jx/bin:$PATH
export QE_MPIRUN=$(command -v mpirun)
export QE_PSEUDO=/orcd/pool/004/sdajani/fecr_qe/pseudo
export PY=/orcd/data/omarabu/001/sdajani/mamba/bin/python3
export OMP_NUM_THREADS=1
