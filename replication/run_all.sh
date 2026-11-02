#!/usr/bin/env bash
# Regenerate the stored studies (strict-limit validation, tolerance grid,
# cross-language agreement, reference-band checks, calibration and power
# screens) with the 0.4.0 packages.
set -euo pipefail
cd "$(dirname "$0")/.."
PY="${PYTHON:-python3}"
export PYTHONPATH="$PWD/python/src"
"$PY" replication/scripts/run_study.py --B 1999 --seed 42 --calibration-B 999 --calibration-reps 300 \
  --power-B 499 --power-reps 150 --workers 4 --out data
"$PY" replication/scripts/run_cross_language_agreement.py --B 1999 --seed 2026 --n-grid 50 \
  --out data/cross_language_agreement.csv
"$PY" replication/scripts/run_diagnostic_study.py --seed 20260609 --illustration-M 999 --sample-M 199 \
  --dataset-reps 15 --probe-reps 100 --workers 4 --out data
"$PY" replication/scripts/write_comparison_log.py --out data/qnca_comparison.log --bottleneck-grid 101
for f in frontier_geometry pi_reduction outlier_fragility power_curve; do
  Rscript "replication/R/figures/$f.R"
done
echo "regeneration complete"
