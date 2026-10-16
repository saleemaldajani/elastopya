#!/usr/bin/env python3
"""Per tensor: highest per-grain frequency and % of grains with a pseudo-SAW above 609 MHz for relative PSAW
intensity thresholds i_psaw/i_saw >= 1e-4, 1e-3, 2e-3 (and with no threshold). Writes thresholds.csv."""
import numpy as np

import common as cm

THR = (1e-4, 1e-3, 2e-3)
REF = {"SS304 (pristine CASS, 209/133/121, 7950)": cm.ROOT / "data" / "saw" / "SS304_10k_python.csv",
       "qeaexp_x0.35 (QE random alloy; paper)": cm.ROOT / "data" / "saw" / "qeaexp_x0.35.csv",
       "qeaexpdec_x0.35 (QE decomposed; paper)": cm.ROOT / "data" / "saw" / "qeaexpdec_x0.35.csv",
       "md_x0.36 (MD random; paper)": cm.ROOT / "data" / "saw" / "md_x0.36.csv",
       "mddec_x0.36 (MD decomposed; paper)": cm.ROOT / "data" / "saw" / "mddec_x0.36.csv",
       "lit_x0.35 (EMTO random; paper)": cm.ROOT / "data" / "saw" / "lit_x0.35.csv",
       "litdec_x0.36 (EMTO decomposed; data/saw)": cm.ROOT / "data" / "saw" / "litdec_x0.36.csv",
       "lit_x0.90 (EMTO Cr-10Fe; paper)": cm.ROOT / "data" / "saw" / "lit_x0.90.csv",
       "qeaexp_x0.95 (QE Cr-5Fe; paper)": cm.ROOT / "data" / "saw" / "qeaexp_x0.95.csv",
       "qeaexp_x1.00 (QE AF Cr a_exp; paper)": cm.ROOT / "data" / "saw" / "qeaexp_x1.00.csv"}


def stats(path):
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    fs, fp = d[:, 0] / cm.LAM / 1e6, d[:, 1] / cm.LAM / 1e6
    r = d[:, 3] / d[:, 2]
    out = {}
    for t in THR + (0.0,):
        keep = np.isfinite(fp) & (r >= t)
        hi = np.where(keep, np.fmax(fs, fp), fs)
        out[t] = (np.nanmax(hi), 100 * np.mean(keep & (fp > cm.BAND)))
    return out


def main():
    paths = dict(REF)
    for p in sorted(cm.OUT.glob("*.csv")):
        paths[p.stem] = p
    lines = ["tag," + ",".join(f"fmax_thr{t:g},pct_psaw_gt609_thr{t:g}" for t in THR + (0.0,))]
    print(f"{'tensor':44s}" + "".join(f"  thr={t:g}: fmax  %>609" for t in THR + (0.0,)))
    for k, p in paths.items():
        s = stats(p)
        lines.append(k.replace(",", ";") + "," + ",".join(f"{s[t][0]:.1f},{s[t][1]:.2f}" for t in THR + (0.0,)))
        print(f"{k:44s}" + "".join(f"  {s[t][0]:12.1f} {s[t][1]:6.2f}" for t in THR + (0.0,)))
    (cm.HERE / "thresholds.csv").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
