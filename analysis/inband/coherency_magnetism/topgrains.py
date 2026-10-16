#!/usr/bin/env python3
"""Sub-grid evaluation of the highest per-grain (P)SAW frequency on the top-K grains of a 10k run.

For small tensor changes (coherency: < 1 %), the 10k-grain maximum moves by less than the solver's ~0.6 MHz grid
step, so the top-K grains of the unperturbed run are re-evaluated with density dithering (common.grain_fmax_dither)
for the base and the perturbed tensors (general C_ijkl allowed, e.g. a single tetragonal laminate variant).
"""
from functools import partial
from multiprocessing import Pool

import numpy as np

import common as cm


def _one(euler, C, rho):
    return cm.grain_fmax_dither(C, rho, euler)


def top_ids(base_csv, k=40):
    f = cm.highest(base_csv)
    return np.argsort(-np.nan_to_num(f, nan=-1))[:k]


def eval_top(ids, C_ijkl_pa, rho, workers=6):
    ang = np.loadtxt(cm.ANGLES, delimiter=",")[ids]
    with Pool(workers) as p:
        return np.array(p.map(partial(_one, C=C_ijkl_pa, rho=rho), list(ang)))
