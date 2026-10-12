End-to-end pipeline test only (16-atom SQS, Fe50Cr50, k 2x2x2, 40 Ry). Not a result.
The cutoff is far below the Fe PAW minimum (64 Ry wfc / 782 Ry rho), so vc-relax collapsed the cell
to 10.2 A^3/atom (-12 %; the final fresh-basis SCF still shows P = -282 kbar) and C11 = 568 GPa is
an artefact of that compression. Production runs use 60 Ry, validated on bcc Fe (a = 2.827 A,
C11 within 1 GPa of 80 Ry), and elastic_sqs_qe.py now stops if the relaxed volume is unphysical.
