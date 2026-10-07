#!/usr/bin/env bash
# repair_v9 sensitivity items S1 (IMF ceiling) and S2 (subgroup labels), as
# pre-registered in provenance/repair_v9_s1_s2_prereg.json (commit 3c571f6).
# Run detached from the repository root.
set -euo pipefail
export PYTHONPATH=scripts WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging
export CYGOB2_CHAIN=repair_v9 OMP_NUM_THREADS=1
PY="${PYTHON_BIN:-python3}"
LOG=data/processed/repair_v9_logs
( "$PY" scripts/repair_v9_s1_imf_ceiling.py > "$LOG/s1.log" 2>&1 ) & s1=$!
( "$PY" scripts/repair_v9_s2_label_uncertainty.py responsibilities && \
  "$PY" scripts/repair_v9_s2_label_uncertainty.py masses && \
  "$PY" scripts/repair_v9_s2_label_uncertainty.py fit && \
  "$PY" scripts/repair_v9_s2_label_uncertainty.py ledger ) > "$LOG/s2.log" 2>&1 & s2=$!
wait $s1; wait $s2
"$PY" scripts/repair_v9_s1_s2_score.py
echo "=== S1/S2 complete $(date '+%F %T')"
