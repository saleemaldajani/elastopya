#!/bin/bash
# Rebuild every derived dataset and figure set after new MD / SAW results, then sync the repo.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-python3}
$PY scripts/collect_md_random.py > /dev/null
$PY scripts/collect_md_decomposition.py > /dev/null
$PY scripts/phase_compositions.py > /dev/null
$PY scripts/visualize_decomposition.py
$PY scripts/make_figures.py > /dev/null
$PY scripts/figures_with_strips.py > /dev/null
$PY scripts/make_figures_kdense.py | sed -n '/figures ->/,$p'
echo REGEN_DONE
