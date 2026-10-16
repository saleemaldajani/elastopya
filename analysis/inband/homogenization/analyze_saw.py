#!/usr/bin/env python3
"""Task 3 (SAW): per-grain highest mode and in-band fraction for every saw_runs/<model>_<scheme>.csv, plus a
grain-by-grain check of the Voigt upper bound. Writes fecr36_saw_summary.csv / .txt."""
import csv
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RUN = HERE / "saw_runs"
BAND = 609.0
ORDER = ["random", "voigt", "reuss", "hill", "laminate_tet", "laminate_V3", "sc_sphere", "fft_axes", "fft_iso"]


def load(tag):
    d = np.loadtxt(RUN / f"{tag}.csv", delimiter=",", skiprows=1)
    keep = d[:, 3] / d[:, 2] >= 2e-3
    saw = d[:, 0] / 6.4e-6 / 1e6
    psaw = np.where(keep, d[:, 1] / 6.4e-6 / 1e6, np.nan)
    top = np.where(keep, np.fmax(saw, psaw), saw)
    return saw, psaw, top


def main():
    tens = {(r["model"], r["scheme"]): r for r in csv.DictReader(open(HERE / "fecr36_tensors.csv"))}
    rows, lines = [], []
    lines.append(f"{'model':5s} {'scheme':14s} {'A_cub':>6s} {'f_iso':>6s} {'maxSAW':>7s} {'maxTop':>7s} {'p99':>6s} "
                 f"{'median':>7s} {'%>609':>6s} | vs Voigt: %grains SAW>V+1MHz, max(SAW-V), %grains top>V+1MHz, max(top-V)")
    for model in ("qe", "emto"):
        if not (RUN / f"{model}_voigt.csv").exists():
            continue
        sv, _, tv = load(f"{model}_voigt")
        for sch in ORDER:
            tag = f"{model}_{sch}"
            if not (RUN / f"{tag}.csv").exists():
                continue
            s, p, t = load(tag)
            ds, dt = s - sv, t - tv
            r = dict(model=model, scheme=sch, A_cub=float(tens[(model, sch)]["A_cub"]), f_iso=float(tens[(model, sch)]["f_iso"]),
                     max_saw=np.nanmax(s), max_top=np.nanmax(t), p99_top=np.nanpercentile(t, 99), median_top=np.nanmedian(t),
                     pct_inband=100 * np.mean(t > BAND), pct_saw_above_voigt=100 * np.mean(ds > 1.0), max_saw_minus_voigt=np.nanmax(ds),
                     pct_top_above_voigt=100 * np.mean(dt > 1.0), max_top_minus_voigt=np.nanmax(dt))
            rows.append(r)
            lines.append(f"{model:5s} {sch:14s} {r['A_cub']:6.3f} {r['f_iso']:6.1f} {r['max_saw']:7.1f} {r['max_top']:7.1f} "
                         f"{r['p99_top']:6.1f} {r['median_top']:7.1f} {r['pct_inband']:6.2f} | {r['pct_saw_above_voigt']:5.2f} "
                         f"{r['max_saw_minus_voigt']:+6.1f} {r['pct_top_above_voigt']:5.2f} {r['max_top_minus_voigt']:+6.1f}")
    with open(HERE / "fecr36_saw_summary.csv", "w", newline="") as fh:
        wr = csv.DictWriter(fh, list(rows[0])); wr.writeheader()
        wr.writerows([{k: (f"{v:.3f}" if isinstance(v, float) else v) for k, v in r.items()} for r in rows])
    txt = "\n".join(lines)
    print(txt)
    (HERE / "fecr36_saw_summary.txt").write_text(txt + "\n")


if __name__ == "__main__":
    main()
