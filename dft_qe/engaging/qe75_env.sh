# Quantum ESPRESSO 7.5 built in qe75/ with the same options as the local build
source /etc/profile >/dev/null 2>&1
module load gcc/12.2.0 openmpi/4.1.4 openblas/0.3.26 fftw/3.3.10
export QE_PW=/orcd/pool/004/sdajani/fecr_qe/qe75/q-e-qe-7.5/build/bin/pw.x
export QE_MPIRUN=$(which mpirun)
export QE_PSEUDO=/orcd/pool/004/sdajani/fecr_qe/pseudo
export PY=/orcd/data/omarabu/001/sdajani/mamba/bin/python3
export OMP_NUM_THREADS=1
