#!/usr/bin/env bash
# repair_v4 chain: truth-side joint age--k fit (issue #1c, step 3c).
#
# repair_v4 is a WP5-only version.  WP3 extinction and WP4 ages and masses are
# consumed from repair_v3 unchanged -- the model change is entirely in the WP5
# injection truth age, so re-running the upstream would only add Monte-Carlo
# noise.  Every repair_v1..v3 artifact is preserved byte-identical.
#
# Evidence for this design: provenance/wp5_age_scan_execution.json (gate G2),
# provenance/wp5_joint_fit_baseline_check_execution.json (gate 3b).
#
# Runtime: roughly 3 hours for the node injections plus about an hour for the
# 54-branch fit.  Run from the repository root.
set -euo pipefail

export PYTHONPATH=scripts
export WP_REPAIR_VERSION=repair_v3
export WP3_ANCHOR_PRIOR_MODE=variogram

echo "=== [1/2] WP5 truth-age node injections ==="
python3 scripts/wp5_injections_agenodes.py --output-version repair_v4

echo "=== [2/2] WP5 joint age-k fit and gate G3 ==="
python3 scripts/wp5_fit_imf_joint.py \
    --upstream-version repair_v3 --wp5-version repair_v4

echo "=== repair_v4 chain complete ==="
