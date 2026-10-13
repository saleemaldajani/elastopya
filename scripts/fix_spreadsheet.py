#!/usr/bin/env python3
"""Corrected copy of the working Google Sheet (data/spreadsheet/).

Usage: python scripts/fix_spreadsheet.py
  in : data/spreadsheet/FeCr_working_sheet_original.xlsx   (Google Sheets export, 12 Oct 2026)
  out: data/spreadsheet/FeCr_working_sheet_fixed.xlsx

Only the faulty formulas are rewritten, directly in the sheet XML, so all charts and formatting
survive (openpyxl would drop the 82 charts). Changed cells lose their stale cached value and the
workbook is flagged to recalculate fully on open (Excel, Numbers, Google Sheets).

'Analyzed data'
  1. C11 (col B) and B (col K), rows 8-21 and 23: read 'Raw data'!H one row down from x~0.11,
     and pure Cr (row 23) read row 8. Now the same row.
  2. E (col W), literature rows 2-21 and 23: 9KT/(3K+Q) used G_Reuss (Q) in the denominator;
     now 9BG/(3B+G) with the Hill G (T), as in the DFT columns.
  3. Poisson ratio (col AA), MD rows 24-33: divided by $K$23 (pure-Cr literature B); now each
     row's own K.
  4. Density (col AS), MD rows 28-33: used phi (col A) as x_Cr; now equiatomic Fe50Cr50.
  5. Unallowed-range minimum (col DY): 0.627 GHz (Acta Fig. 9) -> 0.609 GHz (Acta text), the
     value used in the manuscript (decision of 12 Oct 2026).
'Enthalpy of demixing'
  6. New columns AF-AG: three-point second derivative of the PWSCF-PAW demixing free energy on
     the non-uniform grid (duplicate x = 0.65 skipped). The original columns O-Q (two-point,
     wrong spacing) are kept because Fig. 1b shows them for comparison.
"""
import re
import shutil
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "data" / "spreadsheet" / "FeCr_working_sheet_original.xlsx"
DST = ROOT / "data" / "spreadsheet" / "FeCr_working_sheet_fixed.xlsx"
ANALYZED, DEMIX = "xl/worksheets/sheet7.xml", "xl/worksheets/sheet9.xml"


def set_formula(xml, ref, formula):
    """Give cell `ref` an explicit (non-shared) formula and drop its cached value."""
    pat = re.compile(r'<c r="%s"([^>]*?)(/>|>(.*?)</c>)' % ref, re.S)
    m = pat.search(xml)
    assert m, ref
    attrs = re.sub(r'\s+t="[^"]*"', "", m.group(1))
    return xml[:m.start()] + f'<c r="{ref}"{attrs}><f>{formula}</f></c>' + xml[m.end():]


def set_shared_master(xml, ref, formula):
    """Rewrite the master formula of a shared-formula group (children follow it)."""
    pat = re.compile(r'(<c r="%s"[^>]*><f t="shared"[^>]*>)([^<]*)(</f>)(<v>[^<]*</v>)?' % ref)
    m = pat.search(xml)
    assert m, ref
    return xml[:m.start()] + m.group(1) + formula + m.group(3) + xml[m.end():]


def fix_analyzed(s):
    for r in list(range(8, 22)) + [23]:
        raw = r - 1 if r == 23 else r          # row 23 holds 'Raw data' row 22 (pure Cr)
        s = set_formula(s, f"B{r}", f"(3*'Raw data'!H{raw}+4*'Raw data'!B{raw})/3")
        s = set_formula(s, f"K{r}", f"'Raw data'!H{raw}")
    s = set_shared_master(s, "W2", "9*K2*T2/(3*K2+T2)")
    s = set_formula(s, "W23", "9*K23*T23/(3*K23+T23)")
    for r in range(24, 34):
        s = set_formula(s, f"AA{r}", f"(3*K{r}-W{r})/(6*K{r})")
    for r in range(28, 34):
        s = set_formula(s, f"AS{r}", "(0.5*B$39+0.5*B$38)*1000")
    s = set_shared_master(s, "DY2", "609/1000")
    return s


def fix_demix(s):
    rows = {}
    for r in range(2, 35):
        m = re.search(r'<c r="L%d"[^>]*><v>([^<]*)</v>' % r, s)
        if m:
            rows[r] = float(m.group(1))
    seen, grid = set(), []
    for r, x in sorted(rows.items(), key=lambda t: (t[1], t[0])):
        if round(x, 9) not in seen:              # skip the duplicated x = 0.65 row
            seen.add(round(x, 9)); grid.append(r)
    cells = {1: ('<c r="AF1" t="inlineStr"><is><t>X_Cr (unique grid)</t></is></c>'
                 '<c r="AG1" t="inlineStr"><is><t>Second derivative PWSCF-PAW, three-point non-uniform '
                 '(corrected; replaces col P)</t></is></c>')}
    for i, r in enumerate(grid):
        c = f'<c r="AF{r}"><f>L{r}</f></c>'
        if 0 < i < len(grid) - 1:
            p, n = grid[i - 1], grid[i + 1]
            c += (f'<c r="AG{r}"><f>2*((M{n}-M{r})/(L{n}-L{r})-(M{r}-M{p})/(L{r}-L{p}))/(L{n}-L{p})</f></c>')
        cells[r] = c
    for r, c in cells.items():
        s = re.sub(r'(<row r="%d"[^>]*>.*?)(</row>)' % r, lambda m: m.group(1) + c + m.group(2), s, count=1, flags=re.S)
    return s


def main():
    shutil.copy2(SRC, DST)
    with zipfile.ZipFile(SRC) as zin:
        items = {n: zin.read(n) for n in zin.namelist()}
    items[ANALYZED] = fix_analyzed(items[ANALYZED].decode()).encode()
    items[DEMIX] = fix_demix(items[DEMIX].decode()).encode()
    wb = items["xl/workbook.xml"].decode()
    wb = re.sub(r"<calcPr[^>]*/>", '<calcPr fullCalcOnLoad="1"/>', wb)
    items["xl/workbook.xml"] = wb.encode()
    with zipfile.ZipFile(DST, "w", zipfile.ZIP_DEFLATED) as zout:
        for n, b in items.items():
            zout.writestr(n, b)
    print("wrote", DST)


if __name__ == "__main__":
    main()
