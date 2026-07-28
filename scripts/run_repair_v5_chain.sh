#!/usr/bin/env bash
# repair_v5 chain: kriged anchor prior mean + truth-side joint age-k fit.
#
# ADOPTED WORKING VERSION (2026-07-28).  See provenance/wp3_kriging_adoption.json
# for the decision record and reports/WP3_KRIGED_PRIOR_repair_v5.md for the
# evidence.
#
# What differs from repair_v3/v4:
#   WP3  ANCHOR_PRIOR_MODE=kriging -- the anchor spatial prior's MEAN is now a
#        simple-kriging estimate using the variogram already fitted in
#        repair_v3, instead of the plain median of the 8 nearest anchors.  The
#        prior WIDTH is unchanged.  No new parameter.
#   WP5  truth-side joint age-k marginalization (inherited from repair_v4).
#
# Every repair_v1..v4 artifact is preserved byte-identical; the older prior is
# still reachable with WP3_ANCHOR_PRIOR_MODE=variogram|global.
#
# Runtime: about 15 min WP3+WP4, about 30 min for 37 node injections, about
# 20 min for the 54-branch fit.  Run from the repository root.
set -euo pipefail

export PYTHONPATH=scripts
export WP_REPAIR_VERSION=repair_v5
export WP3_ANCHOR_PRIOR_MODE=kriging

echo "=== [1/5] WP3 extinction (kriged anchor prior) ==="
python3 scripts/wp3_extinction_repair.py

echo "=== [2/5] WP4 ages ==="
python3 scripts/wp4_fit_ages.py --repair-version repair_v5

echo "=== [3/5] WP4 mass posteriors ==="
python3 scripts/wp4_mass_posteriors_repair.py

echo "=== [4/5] WP5 truth-age node injections ==="
python3 scripts/wp5_injections_agenodes.py --output-version repair_v5

echo "=== [5/5] WP5 joint age-k fit and gate ==="
python3 scripts/wp5_fit_imf_joint.py \
    --upstream-version repair_v5 --wp5-version repair_v5 --compare-version repair_v4

echo "=== reporting ==="
python3 scripts/wp5_report.py --wp5-version repair_v5
python3 scripts/wp5_repair_v4_finalize.py --wp5-version repair_v5 --upstream-version repair_v5

echo "=== repair_v5 chain complete ==="
