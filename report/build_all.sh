#!/usr/bin/env bash
# Rebuild the master-architecture report end to end.
#
#   figures_src/*.py  ->  figures/*.pdf        (schematics)
#   build_report.py   ->  master_architecture_report.tex
#                          (numbers harvested from the notebook + fixtures)
#   pdflatex          ->  master_architecture_report.pdf
#   copy              ->  main.tex             (the name Overleaf expects)
#
# Publishing to Dropbox/Overleaf is a separate step: copy main.tex and
# figures/*.pdf into  <Dropbox>/Apps/Overleaf/dbhe_master_architecture/ .
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p figures

for f in master water_circuit topologies qmap; do
    python3 "figures_src/$f.py" >/dev/null
done
echo "figures rebuilt"

python3 build_report.py

pdflatex -interaction=nonstopmode master_architecture_report.tex >/dev/null
pdflatex -interaction=nonstopmode master_architecture_report.tex >/dev/null
pdflatex -interaction=nonstopmode master_architecture_report.tex >/dev/null
n_err=$(grep -c '^!' master_architecture_report.log || true)
echo "latex errors: $n_err"
[ "$n_err" -eq 0 ] || { echo "REPORT DID NOT COMPILE CLEANLY"; exit 1; }

cp master_architecture_report.tex main.tex
echo "main.tex ready to publish"
