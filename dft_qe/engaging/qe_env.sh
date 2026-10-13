# Quantum ESPRESSO 7.3.1 on Engaging (ORCD community module), used by the job scripts
source /etc/profile >/dev/null 2>&1
module load community-modules
module load quantum-espresso/7.3.1
export QE_PW=$(which pw.x)
export QE_MPIRUN=/orcd/software/core/001/spack/pkg/openmpi/4.1.4/zahpnmk/bin/mpirun
export QE_PSEUDO=/orcd/pool/004/sdajani/fecr_qe/pseudo
export PY=/orcd/data/omarabu/001/sdajani/mamba/bin/python3
export OMP_NUM_THREADS=1
