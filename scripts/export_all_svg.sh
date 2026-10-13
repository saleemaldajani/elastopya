#!/bin/bash
# Regenerate every paper / SM / response figure with an SVG copy (editable text). Run from the repo root.
set -u
P=/private/tmp/claude-501/-Users-saleem-MIT-research-Scripta/37f8fe43-8e5f-4d11-883b-ebb44fa03560/scratchpad/pdfenv/bin/python
SM=manuscript/figures/kdense_scripta/supplementary
E=dft_qe/engaging
step() { echo "== $(date +%H:%M) $*"; }
step "main figures (KDense, all variants)"; (cd scripts && nice -n 10 $P make_figures_kdense.py > /dev/null 2>&1) && echo ok
step "Fig. 4 mechanism"; (cd scripts && nice -n 10 $P make_fig4_mechanism.py > /dev/null 2>&1) && echo ok
step "SM Fe-rich elastic"; (cd scripts && nice -n 10 $P make_si_figures.py > /dev/null 2>&1) && echo ok
step "response figures"; (cd scripts && nice -n 10 $P make_response_figures.py > /dev/null 2>&1) && echo ok
step "S_dos_sqs"; nice -n 10 $P dft_qe/dos_plot.py $SM/S_dos_sqs.png $(ls -d $E/dos/x* | grep -v _aexp | while read d; do ls $d/pdos?.pdos_tot >/dev/null 2>&1 && echo $d; done) > /dev/null 2>&1 && echo ok
step "S_phonons"; nice -n 10 $P dft_qe/phonons.py plot $SM/S_phonons.png $E/phonon/Fe_fm $E/phonon/Cr_nm_v5_k24 $E/phonon/Cr_afm > /dev/null 2>&1 && echo ok
step "S_phonons_Cr_nm_v* - svg only, PNG and PDF untouched"
T=$(mktemp -d)
for v in $SM/S_phonons_Cr_nm_v*.png; do
  b=$(basename $v .png); d=${b#S_phonons_}; d=$(ls -d $E/phonon/${d%%_*}_${d#*_} 2>/dev/null || true)
  case $b in *v1_default-mixing) d=$E/phonon/Cr_nm;; *v2_mix) d=$E/phonon/Cr_nm_v2_mix;; *v3_k20) d=$E/phonon/Cr_nm_v3_k20;;
             *v4_deg04) d=$E/phonon/Cr_nm_v4_deg04;; *v5_k24) d=$E/phonon/Cr_nm_v5_k24;; *v6_a2843) d=$E/phonon/Cr_nm_v6_a2843;; esac
  nice -n 10 $P dft_qe/phonons.py plot $T/$b.png $d --allow-imaginary > /dev/null 2>&1 && cp $T/$b.svg $SM/$b.svg && echo "  $b.svg"
done
step "S_fermi_Fe, S_fermi_pure"; nice -n 10 $P dft_qe/fermi_slices.py korzhavyi $SM/S_fermi_Fe.png $E/fermi/Fe_fm > /dev/null 2>&1 && \
  nice -n 10 $P dft_qe/fermi_slices.py plot $SM/S_fermi_pure.png $E/fermi/Fe_fm $E/fermi/Cr_nm $E/fermi/Cr_afm > /dev/null 2>&1 && echo ok
step "S_bands_pure"; nice -n 10 $P dft_qe/bands_path.py pure-plot $SM/S_bands_pure.png $E/bandpath/Fe_fm $E/bandpath/Cr_nm $E/bandpath/Cr_afm > /dev/null 2>&1 && echo ok
step "summary"; for f in $(grep -ho 'figures/kdense_scripta/[^}]*\.pdf' manuscript/finalpaper*.tex response_to_mike*.tex | sort -u); do
  s=manuscript/${f%.pdf}.svg; [ -f "$s" ] || s=${f%.pdf}.svg; [ -f "manuscript/${f%.pdf}.svg" ] || [ -f "${f%.pdf}.svg" ] || echo "  MISSING svg: $f"; done
echo DONE
