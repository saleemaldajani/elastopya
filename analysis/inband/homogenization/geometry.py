#!/usr/bin/env python3
"""Synthetic periodic spinodal microstructures (voxel masks of alpha') for FFT homogenization.

Cahn's picture of spinodal decomposition: the composition field is a superposition of concentration waves with
|k| ~ k0; in an elastically anisotropic cubic crystal with soft <100> the waves are restricted to (or near) <100>.
  kind="axes"   : k exactly along +-[100],[010],[001]  ->  c(r) = f(x) + g(y) + h(z) (interconnected 'brick' structure)
  kind="cone"   : k within a Gaussian angular spread `spread_deg` around the nearest <100>
  kind="iso"    : k isotropically distributed (no elastic anisotropy; cf. the MD/MC cells)
The field is thresholded at its (1 - w) quantile so that the alpha' volume fraction is w (to one voxel).
"""
import numpy as np


def spinodal_mask(n, w, kind="axes", k0_cycles=6, dk=1.0, nwaves=300, spread_deg=15.0, seed=0):
    rng = np.random.default_rng(seed)
    g = np.arange(n) / n
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    c = np.zeros((n, n, n))
    for _ in range(nwaves):
        kmag = k0_cycles + dk * rng.standard_normal()
        if kind == "axes":
            d = np.eye(3)[rng.integers(3)]
        else:
            if kind == "iso":
                d = rng.standard_normal(3)
            else:
                ax = np.eye(3)[rng.integers(3)]
                d = ax + np.tan(np.radians(spread_deg)) * rng.standard_normal(3) * (1 - ax)
            d /= np.linalg.norm(d)
        kv = np.round(kmag * d)                   # integer wave vector -> exactly periodic in the cell
        if not kv.any():
            continue
        c += np.cos(2 * np.pi * (kv[0] * X + kv[1] * Y + kv[2] * Z) + rng.uniform(0, 2 * np.pi))
    thr = np.quantile(c, 1 - w)
    return c > thr
