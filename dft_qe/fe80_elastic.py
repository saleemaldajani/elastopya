#!/usr/bin/env python3
"""Elastic constants of the Fe-80Cr SQS (PBE-relaxed cell, Cr antiparallel to Fe) from converged steps only
(added 2026-10-17).

The four strain steps were run twice (elastic_sqs_qe.py, QE_STRAIN_ONLY, per-site starting moments):
  runs/sqs40_x0.80_k442_e60_crafm2        conv_thr 1e-8: s_xx_p, s_yz_p, s_yz_m converged; s_xx_m sat at a 2e-6 Ry
                                          floor and stopped at 600 iterations
  runs/sqs40_x0.80_k442_e60_crafm2_thr3   conv_thr 3e-6 for s_xx_m, s_yz_p, s_yz_m (s_xx_p copied from the 1e-8 set)
The shear pair, converged at both thresholds, measures the effect of the looser threshold: sigma_yz agrees to 0.05 kbar
(C44 107.3 vs 107.5 GPa), but every diagonal stress of the 3e-6 runs is ~3 kbar higher. The 3e-6 s_xx_m therefore
biases C11 and C12 of the central difference (364 and 111 GPa, result.json of the _thr3 directory).

Here C11 and C12 come from the forward difference between the converged s_xx_p and the converged relaxed cell
(last BFGS step of relax.out), and C44 from the converged central yz pair. The yz pair checks linearity: the relaxed
sigma_yz (-1.86 kbar) lies at the midpoint of the +-strain values to 0.01 kbar, and the forward and backward
responses agree to 0.2 %, so the one-sided difference carries no appreciable nonlinearity error.

  fe80_elastic.py   -> runs/sqs40_x0.80_k442_e60_crafm2/result.json (read by make_figures.qe_dft / qe_mixing)
"""
import json
import re
from pathlib import Path

import numpy as np

RUN = Path(__file__).resolve().parent / "engaging" / "runs"
A, B = RUN / "sqs40_x0.80_k442_e60_crafm2", RUN / "sqs40_x0.80_k442_e60_crafm2_thr3"
EPS = 0.005


def stress_kbar(text):
    if "convergence NOT achieved" in text.split("End of BFGS")[0] if "End of BFGS" in text else "convergence NOT achieved" in text:
        raise SystemExit("unconverged SCF in a step used here")
    blk = text.split("total   stress")[-1].splitlines()[1:4]
    return np.array([[float(v) for v in l.split()[3:6]] for l in blk])


def main():
    relax = (A / "relax.out").read_text()
    s0 = stress_kbar(relax.split("End of BFGS Geometry Optimization")[0])          # converged relaxed cell
    sxp = stress_kbar((A / "s_xx_p.out").read_text())
    syp, sym = stress_kbar((A / "s_yz_p.out").read_text()), stress_kbar((A / "s_yz_m.out").read_text())
    # QE prints compressive stress as positive; C = -d sigma / d eps (in GPa: kbar / 10)
    c11 = -(sxp[0, 0] - s0[0, 0]) / EPS / 10
    c12 = -0.5 * ((sxp[1, 1] - s0[1, 1]) + (sxp[2, 2] - s0[2, 2])) / EPS / 10
    c44 = -(syp[1, 2] - sym[1, 2]) / (2 * EPS) / 10
    lin = abs((syp[1, 2] - s0[1, 2]) / (s0[1, 2] - sym[1, 2]) - 1)
    mid = abs(0.5 * (syp[1, 2] + sym[1, 2]) - s0[1, 2])
    # bias of the 3e-6 runs: diagonal stresses of the yz steps, 3e-6 minus 1e-8
    off = np.mean([np.mean(np.diag(stress_kbar((B / f"s_yz_{s}.out").read_text()) - stress_kbar((A / f"s_yz_{s}.out").read_text())))
                   for s in "pm"])
    c44_b = -(stress_kbar((B / "s_yz_p.out").read_text())[1, 2] - stress_kbar((B / "s_yz_m.out").read_text())[1, 2]) / (2 * EPS) / 10
    mags = re.findall(r"total magnetization\s+=\s+([-\d.]+)", relax.split("End of BFGS")[0])
    vol = float(re.findall(r"new unit-cell volume =\s+[\d.]+ a\.u\.\^3 \(\s+([\d.]+) Ang\^3", relax)[-1])
    res = {"x_cr": 0.8, "k": [4, 4, 2], "ecut_Ry": 60.0, "volume_A3": vol,
           "moment_muB_per_atom": float(mags[-1]) / 40, "fixed_a": None,
           "C11": round(c11, 1), "C12": round(c12, 1), "C44": round(c44, 1),
           "method": "C11, C12: forward difference s_xx_p (1e-8) vs relaxed cell; C44: central yz pair (1e-8). See fe80_elastic.py",
           "checks": {"yz_linearity_rel": round(lin, 4), "yz_midpoint_minus_relaxed_kbar": round(mid, 3),
                      "thr3e-6_diag_offset_kbar": round(off, 2), "C44_thr3e-6": round(c44_b, 1),
                      "thr3e-6_central_C11_C12": json.loads((B / "result.json").read_text()) and
                      [round(json.loads((B / "result.json").read_text())[k], 1) for k in ("C11", "C12")]}}
    (A / "result.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
