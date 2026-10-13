# Working spreadsheet: original and corrected copies

- `FeCr_working_sheet_original.xlsx`: export of the working Google Sheet, 12 Oct 2026, unchanged.
- `FeCr_working_sheet_fixed.xlsx`: the same workbook with the formula errors below corrected,
  produced by `scripts/fix_spreadsheet.py`. Charts and formatting are untouched. The file
  recalculates when opened (Excel, Numbers, or File > Import in Google Sheets).

The manuscript figures do not read this workbook; `scripts/make_figures.py` uses the CSV
snapshots in `data/` and applies the same corrections in code.

| # | Sheet, cells | Original | Corrected | Effect |
|---|---|---|---|---|
| 1 | Analyzed data, B8:B21, K8:K21, B23, K23 (literature C11 and B) | `'Raw data'!H` read one row down from x ≈ 0.11; pure Cr read row 8 | same row as the composition | pure-Cr C11 388 → 439 GPa; all literature C11, B, f_SAW for x ≥ 0.11 |
| 2 | Analyzed data, W2:W21, W23 (literature E) | `9*K*T/(3*K+Q)`: Reuss G in the denominator | `9*K*T/(3*K+T)`: Hill G, as in the DFT columns | E, ν and f_SAW of the literature series |
| 3 | Analyzed data, AA24:AA33 (MD Poisson ratio) | divided by `$K$23` (pure-Cr literature B) | each row's own K | ν and f_SAW of the MD φ series |
| 4 | Analyzed data, AS28:AS33 (MD density) | φ in column A used as x_Cr | equiatomic Fe50Cr50 | density and f_SAW of the MD φ series |
| 5 | Analyzed data, DY (unallowed-range minimum) | 0.627 GHz (Acta Fig. 9) | 0.609 GHz (Acta text) | value adopted in the manuscript (12 Oct 2026) |
| 6 | Enthalpy of demixing, new AF:AG | — | three-point second derivative of the PWSCF–PAW demixing free energy on the non-uniform grid (duplicate x = 0.65 skipped) | replaces column P, whose two-point scheme used the wrong spacing; P is kept because Fig. 1b shows it for comparison |
