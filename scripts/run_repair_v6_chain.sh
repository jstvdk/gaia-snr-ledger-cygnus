#!/usr/bin/env bash
# repair_v6 chain: truth-age isochrone interpolation instead of grid snapping.
#
# Issue #13.  The repair_v4/v5 node rule snapped the nine WP4 age-posterior
# quantiles to native isochrone ages.  The native grid is coarse (0.05 dex), so
# the truth model was a step function of the WP4 age: a 0.0005 Myr shift of the
# posterior could move the node distribution by 0.0546 Myr in 1-Wasserstein
# distance, a 109x amplification, and add or delete a whole node carrying up to
# 1/9 of the weight.  Measured in
# provenance/wp5_node_rule_continuity_execution.json.
#
# repair_v6 keeps the nine quantiles where the posterior puts them, each with
# weight 1/9, and interpolates the truth isochrone between native ages using the
# same bracketing and linear blend the RECOVERY side already used
# (wp4_repair_common._interpolate_age_sequence).  The two sides were previously
# inconsistent; now they are not.  The rule is Lipschitz-1 in the posterior.
#
# repair_v6 is a WP5-ONLY version.  WP3 extinction (kriged prior) and WP4 ages
# and masses are consumed from repair_v5 unchanged -- no upstream step is re-run,
# exactly as repair_v4 consumed repair_v3.  Every repair_v1..v5 artifact stays
# byte-identical.
#
# Predictions were pre-registered BEFORE this chain was run:
#   provenance/wp5_node_interpolation_prereg.json
# Adoption rests on P3 (baseline still passes) and P4 (the change is local),
# both direction-free.  P2 -- the one regressing cell clearing -- is reported
# either way and does NOT gate adoption.
#
# Runtime: about 2.5 h for 162 node injections (37 -> 162; no reuse is possible
# because every node sits at a new age), then about 20 min for the 54-branch
# fit.  Run from the repository root.
set -euo pipefail

export PYTHONPATH=scripts
# Upstream is repair_v5: the kriged WP3 prior and the WP4 products built on it.
export WP_REPAIR_VERSION=repair_v5
export WP3_ANCHOR_PRIOR_MODE=kriging

echo "=== [1/4] WP5 truth-age node injections (interpolated isochrones) ==="
python3 scripts/wp5_injections_agenodes.py --output-version repair_v6

echo "=== [2/4] WP5 joint age-k fit and gate ==="
python3 scripts/wp5_fit_imf_joint.py \
    --upstream-version repair_v5 --wp5-version repair_v6 --compare-version repair_v5

# Verdict stability MUST precede finalization: the finalizer reports gate G3
# under both readings (CUTS section 14.7) and reads the refined one from here.
echo "=== [3/4] pre-registration scoring and verdict stability ==="
python3 scripts/wp5_node_interpolation_outcome.py
python3 scripts/wp5_verdict_stability.py \
    --versions repair_v3 repair_v4 repair_v5 repair_v6 --replicates 400

echo "=== [4/4] reporting, gate finalization, WP6 completeness ==="
python3 scripts/wp5_report.py --wp5-version repair_v6
python3 scripts/wp5_repair_v4_finalize.py --wp5-version repair_v6 --upstream-version repair_v5
python3 scripts/wp6_bright_completeness.py --version repair_v6

echo "=== repair_v6 chain complete ==="
