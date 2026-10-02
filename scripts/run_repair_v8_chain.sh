#!/usr/bin/env bash
# repair_v8 chain: issue #19 -- anchor masses at the repaired ages, per R_V.
#
# What changes vs repair_v7: the 150 spectroscopic-HRD anchor masses.  They were
# read off isochrones at the PRE-REPAIR WP4 ages (unversioned
# wp4_age_posteriors.parquet) and copied onto every R_V branch.  repair_v8 reads
# each branch's own repair_v5 age and de-reddened M_G0 (decision D1, F3(b)).
# Present-day Mass is kept; Mini (F3(a)) is a separate, later change.
#
# Consumed unchanged: WP3 extinction and WP4 ages (repair_v5); the WP5 node and
# WP6 mass-extension injection responses (repair_v7), which never read the
# anchor masses (provenance/issue19_repair_v8_prereg.json, "reuse").
#
# Pre-registered BEFORE this chain ran: provenance/issue19_repair_v8_prereg.json
# (predictions P1-P6, integrity checks I1-I4, adoption rule).
#
# Nothing at repair_v7 or earlier is overwritten.  Scripts that wrote repair_v7
# products without a suffix write _repair_v8 siblings (scripts/chain.py), and
# every versioned step below refuses to overwrite an existing repair_v8 output
# where the script supports it.
#
# Runtime: ~15 min WP4 masses, ~20 min WP5 fit (+20 min replay, in parallel),
# ~40 min WP6/WP7 incl. the 2e6-iteration ledger, ~15 min WP8-WP12.
# Run from the repository root.
set -euo pipefail

export PYTHONPATH=scripts
export WP_REPAIR_VERSION=repair_v5        # WP3/WP4 inputs
export WP3_ANCHOR_PRIOR_MODE=kriging
export CYGOB2_CHAIN=repair_v8             # downstream reads/writes (chain.py)
PY="${PYTHON_BIN:-python3}"

step() { printf '\n=== %s\n' "$1"; }
have() { [ -f "$1" ] && echo "    present: $1 (not regenerated)"; }

step "0  preregistration (refuses to overwrite)"
have provenance/issue19_repair_v8_prereg.json || "$PY" scripts/issue19_prereg.py

step "1  WP4 anchors at the repair_v5 ages, one mass per family x R_V"
have data/processed/wp4_anchor_hrd_repair_v8.parquet || \
  "$PY" scripts/wp4_anchors_hrd.py \
    --age-version repair_v5 --extinction-version repair_v5 --output-version repair_v8
"$PY" scripts/issue19_integrity.py I2

step "2  WP4 mass posteriors (photometric path unchanged)"
have data/processed/wp4_mass_posteriors_repair_v8.parquet || \
  "$PY" scripts/wp4_mass_posteriors_repair.py \
    --anchor-version repair_v8 --output-version repair_v8
"$PY" scripts/issue19_integrity.py I1

step "3  WP5 joint age-k fit on repair_v8 masses, repair_v7 responses; and the replay"
have data/processed/wp5_imf_normalization_repair_v7_replay_issue19.parquet || \
  "$PY" scripts/wp5_fit_imf_joint.py --upstream-version repair_v5 \
    --mass-version repair_v5 --age-version repair_v5 --response-version repair_v7 \
    --wp5-version repair_v7_replay_issue19 --compare-version repair_v6
"$PY" scripts/issue19_integrity.py I3
have data/processed/wp5_imf_normalization_repair_v8.parquet || \
  "$PY" scripts/wp5_fit_imf_joint.py --upstream-version repair_v5 \
    --mass-version repair_v8 --age-version repair_v5 --response-version repair_v7 \
    --wp5-version repair_v8 --compare-version repair_v7

step "4  WP6 census, closure (4.0 Msun floor), attribution, living ledger"
"$PY" scripts/wp6_massive_census.py
"$PY" scripts/wp6_closure_test.py --wp5-version repair_v8 --response-version repair_v7
"$PY" scripts/wp6_closure_attribution.py --wp5-version repair_v8
"$PY" scripts/wp6_ledger.py

step "5  WP7 ledger via the convergence scan (last pass = 2,000,000 iterations)"
"$PY" scripts/wp7_convergence_scan.py
"$PY" scripts/wp7_alpha_headline_adopt.py
"$PY" scripts/wp7_binary_bound.py

step "6  WP8 cross-checks, WP9 verdict, WP11 isotope forecast"
"$PY" scripts/wp8_crosschecks.py
"$PY" scripts/wp9_verdict.py
"$PY" scripts/wp11_isotope_forecast.py --iterations 500000   # the published count (repair_v7)

step "7  reconciliations and the alpha-plausibility record"
"$PY" scripts/wp4_wp5_age_reconciliation.py
"$PY" scripts/wp5_association_mass_reconciliation.py --compare-version repair_v7
"$PY" scripts/wp5_alpha_plausibility.py

step "8  WP12 on repair_v8: new hash record, then the analyses"
have provenance/wp12_revision_prereg_repair_v8.json || \
  "$PY" scripts/wp12_prereg.py --chain-record
"$PY" scripts/wp12_gate_landscape.py
"$PY" scripts/wp12_closure_slopes.py
"$PY" scripts/wp12_scenario_score.py
"$PY" scripts/wp12_neighbour_budget.py

step "9  integrity I4 (repair_v7 untouched) and scoring of P1-P6"
"$PY" scripts/issue19_integrity.py I4
"$PY" scripts/issue19_score.py

echo
echo "=== repair_v8 chain complete"
