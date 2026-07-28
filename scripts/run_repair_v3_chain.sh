#!/usr/bin/env bash
# repair_v3 chain: variogram-based anchor spatial prior.
#
# Re-runs WP3 -> WP4 -> WP5 with a per-star extinction prior width read off the
# fitted anchor variogram instead of one global width per R_V branch.  Fixes
# PROJECT_TRACE.md issue #1b (CygOB2-B mid-window WP5 residual).
#
# Every step writes its own *_repair_v3 artifacts and provenance record; no
# repair_v1 or repair_v2 file is touched.  Run from the repository root.
set -euo pipefail

export PYTHONPATH=scripts
export WP_REPAIR_VERSION=repair_v3
export WP3_ANCHOR_PRIOR_MODE=variogram

echo "=== [1/5] WP3 extinction (variogram prior) ==="
python3 scripts/wp3_extinction_repair.py

echo "=== [2/5] WP4 ages ==="
python3 scripts/wp4_fit_ages.py --repair-version repair_v3

echo "=== [3/5] WP4 mass posteriors ==="
python3 scripts/wp4_mass_posteriors_repair.py

echo "=== [4/5] WP5 injections ==="
python3 scripts/wp5_injections_repair.py --output-version repair_v3

echo "=== [5/5] WP5 IMF fit and gate ==="
python3 scripts/wp5_fit_imf.py --repair-version repair_v3 --wp5-version repair_v3

echo "=== repair_v3 chain complete ==="
