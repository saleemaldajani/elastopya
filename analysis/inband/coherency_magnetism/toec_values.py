#!/usr/bin/env python3
"""Literature inputs for the coherency (acoustoelastic) estimate. GPa throughout. Sources in REPORT.md.

Pressure-derivative convention (Brugger/Thurston effective constants, P = 0, 3B = C11 + 2C12):
    -3B dC11/dP = 2C11 + 2C12 + (C111 + 2C112)
    -3B dC12/dP = -C11 - C12 + (C123 + 2C112)
    -3B dC44/dP = C11 + 2C12 + C44 + (C144 + 2C155)
Checked: Powell & Skove's Fe C111, C112, C123 (derived with Rotter & Smith's dCij/dP) reproduce Rotter & Smith's
dC11/dP and dC12/dP through these relations to 1-2 %.
"""
import numpy as np

# --- alpha-Fe, room temperature -------------------------------------------------------------------------------
# Ledbetter & Reed, J. Phys. Chem. Ref. Data 2, 531 (1973), doi:10.1063/1.3253127
#   Table 20: Powell & Skove, Phys. Rev. 174, 977 (1968): C111 -28.29, C112 -8.00, C123 -6.07 (10^12 dyn/cm^2)
#             Mathur & Sharma (1970), Morse central-force calc.: -16.44, -2.60, -3.00, C144 -3.00, C155 -2.60, C456 -3.00
#   Table 22: dC11/dP, dC12/dP, dC44/dP: Rotter & Smith (1966) 7.51, 5.19, 2.66; Guinan & Beshers (1968) 6.72, 4.58, 2.59
#   Table 5 : Guinan & Beshers (1968) RT C11, C12, C44 = 2.301, 1.346, 1.167 (10^12 dyn/cm^2)
FE_PS = dict(c111=-2829.0, c112=-800.0, c123=-607.0)
FE_MS = dict(c111=-1644.0, c112=-260.0, c123=-300.0, c144=-300.0, c155=-260.0, c456=-300.0)
FE_C_GB = np.array([230.1, 134.6, 116.7])
FE_DCDP = {"Rotter-Smith 1966": (7.51, 5.19, 2.66), "Guinan-Beshers 1968": (6.72, 4.58, 2.59)}

# --- Cr ---------------------------------------------------------------------------------------------------------
# No measured Cr TOEC set was found. Hydrostatic combinations from pressure/volume dependence:
#  IXS: Rodiere et al., arXiv:2607.18015 (2026), Table 1: (P GPa; c11, c44, c') = (0.1; 384, 127, 131), (3.8; 409, 134, 141)
#       (finite difference 0.1 -> 3.8 GPa; some entries there come from fitted force constants)
#  QE (this project): Cr AF and NM, PBE-relaxed vs experimental lattice parameter (runs Cr_*_k16x16x16_e60[_aexp])
CR_IXS = dict(P=(0.1, 3.8), c11=(384.0, 409.0), c44=(127.0, 134.0), cp=(131.0, 141.0))


def combos_from_dcdp(C, dcdp):
    """(K1, K2, K3) = (C111+2C112, C123+2C112, C144+2C155) from RT C and dCij/dP."""
    c11, c12, c44 = C
    b3 = c11 + 2 * c12
    d11, d12, d44 = dcdp
    return np.array([-b3 * d11 - 2 * c11 - 2 * c12, -b3 * d12 + c11 + c12, -b3 * d44 - c11 - 2 * c12 - c44])


def combos_from_volume(C_small, a_small, C_big, a_big):
    """Same combinations from two DFT tensors at lattice parameters a_small < a_big (tangent moduli, -3B d/dP = d/de)."""
    de = np.log(a_big / a_small)
    dC = (np.asarray(C_big) - np.asarray(C_small)) / de
    c11, c12, c44 = C_big
    return np.array([dC[0] - 2 * c11 - 2 * c12, dC[1] + c11 + c12, dC[2] - c11 - 2 * c12 - c44])


def cr_ixs_dcdp():
    dP = CR_IXS["P"][1] - CR_IXS["P"][0]
    d11 = (CR_IXS["c11"][1] - CR_IXS["c11"][0]) / dP
    d44 = (CR_IXS["c44"][1] - CR_IXS["c44"][0]) / dP
    dcp = (CR_IXS["cp"][1] - CR_IXS["cp"][0]) / dP
    return d11, d11 - 2 * dcp, d44


def toec_from_combos(K, d1, d3, c456=0.0):
    """Six TOECs with given hydrostatic combinations K and splits d1 = C111 - C112, d3 = C144 - C155.
    C123 is fixed by K2 once C112 is set; C456 does not enter for diagonal (laminate) strains."""
    c112 = (K[0] - d1) / 3
    c111 = c112 + d1
    c123 = K[1] - 2 * c112
    c155 = (K[2] - d3) / 3
    c144 = c155 + d3
    return (c111, c112, c123, c144, c155, c456)
