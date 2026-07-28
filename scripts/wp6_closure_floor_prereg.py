#!/usr/bin/env python3
"""Issue #17: pre-register the closure-integral lower-limit fix.

Written BEFORE the corrected closure test is run over the grid.

THE DEFECT
----------
scripts/wp6_closure_test.py integrates the forward prediction from 8 Msun:

    predicted observed = k * integral[8, M_turnoff] dM M^-alpha R(obs > 8 | M)

while the observed side counts every member's P(M > 8) regardless of that
member's true mass.  The two sides are therefore not the same quantity.  Stars
whose TRUE mass is below 8 Msun but whose ESTIMATED mass lands above it are
counted on the observed side and are absent from the prediction entirely.

That up-scatter is not small.  R(estimated > 8 | true M) measured on the
accepted repair_v6 response is 0.23 at 7.0 Msun, 0.09 at 6.0 and 0.03 at 5.0.

CUTS section 16.2 SPECIFIES THE CORRECT FORM AND THE CODE CONTRADICTS IT
-----------------------------------------------------------------------
The binding specification writes the integral with no lower limit at all --

    expected observed count above threshold
        = k * integral dM M^-alpha R(observed above threshold | M)

-- and its own justification for forbidding a scalar completeness is that "the
response ALSO SCATTERS MASS ESTIMATES ACROSS THE THRESHOLD", going as far as to
report six cells where "net up-scatter across 8 Msun more than compensates the
recovery loss".  So the specification is right, was right when it was written,
and the implementation truncated the very effect the specification exists to
capture.  This is issue #3's error committed a second time in a different
place: a threshold applied to a quantity that the response smears across that
threshold.

WHY IT WAS NOT CAUGHT
---------------------
The upper limit received all the scrutiny -- it is physically meaningful (the
turnoff), it drove the WP6 mass extension, and it produced issue #14.  The
lower limit looked like a restatement of the 8 Msun supernova threshold rather
than an integration bound, and nothing downstream flagged it because the bias
has the same sign as the excess WP6 was already reporting.

Output: provenance/wp6_closure_floor_prereg.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp6_closure_floor_prereg.py
"""
from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w

# The floor the corrected estimator will adopt.  Chosen on the CONVERGENCE of
# the integral, not on what it does to the closure ratio: on the probe node the
# integral changes by 0.4% between a 4.0 and a 3.0 Msun floor, and R(est > 8 |
# true M) is 0.004 at 4 Msun.  4.0 is also comfortably inside the frozen
# MASS_GRID, so no new injections are required.
ADOPTED_FLOOR_MSUN = 4.0
CONVERGENCE_CHECK_FLOORS = [8.0, 7.0, 6.0, 5.0, 4.0, 3.0]


def main() -> None:
    closure = pd.read_csv(w.TABLES / "wp6_closure.csv")
    baseline = closure[
        closure.family.eq("PARSEC") & closure.R_V.eq(3.1) & closure.alpha.eq(2.3)
    ]

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp6_closure_floor_prereg.py",
        "status": "PREREGISTERED",
        "issue": "#17 — the closure integral truncates at the threshold it is meant to smear across",
        "defect": (
            "scripts/wp6_closure_test.py integrates the forward prediction from "
            "8 Msun, while the observed side counts every member's P(M > 8) "
            "regardless of true mass.  Stars with true mass below 8 Msun whose "
            "estimated mass lands above it are counted as observed and never "
            "predicted, so the two sides are not the same quantity."
        ),
        "specification_says_otherwise": (
            "CUTS section 16.2 writes the integral with NO lower limit, and its "
            "stated reason for forbidding a scalar completeness is that the "
            "response scatters mass estimates ACROSS the threshold — it even "
            "reports six cells where net up-scatter exceeds the recovery loss.  "
            "The specification is correct; the implementation contradicted it."
        ),
        "defect_class": (
            "identical to issue #3 — a threshold applied to a quantity the "
            "response smears across that threshold — committed a second time at "
            "the other end of the same integral"
        ),
        "why_not_caught": (
            "the upper limit absorbed the scrutiny: it is physical, it drove the "
            "mass extension, and it produced issue #14.  The lower limit read as "
            "a restatement of the 8 Msun supernova threshold rather than as an "
            "integration bound, and the bias has the SAME SIGN as the excess WP6 "
            "was already reporting, so nothing downstream contradicted it."
        ),
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "measured_before_registration": {
            "disclosure": (
                "the defect was found by probing ONE node (CygOB2-A PARSEC "
                "R_V=3.1, fifth truth-age node) while the issue #15 injections "
                "ran.  That probe is disclosed here rather than presented as a "
                "prediction: it is what motivated the fix, so it cannot also be "
                "evidence for it."
            ),
            "probe_node_up_scatter_R_est_above_8": {
                "7.00_Msun": 0.2306, "6.00_Msun": 0.0892,
                "5.00_Msun": 0.0294, "4.00_Msun": 0.0039,
            },
            "probe_node_integral_by_floor": {
                "8.0": 0.03163, "7.0": 0.03469, "6.0": 0.03681,
                "5.0": 0.03805, "4.0": 0.03848, "3.0": 0.03862,
            },
            "probe_node_effect_of_floor_8_to_4": "+21.7% predicted observable",
        },
        "adopted_floor_Msun": ADOPTED_FLOOR_MSUN,
        "floor_chosen_how": (
            "on convergence of the integral, NOT on what it does to the closure "
            "ratio: the probe integral changes by 0.4% between a 4.0 and a 3.0 "
            "Msun floor, and R(est > 8 | true 4 Msun) is 0.004.  4.0 also lies "
            "inside the frozen MASS_GRID, so the fix requires no new injections "
            "and touches no accepted artifact."
        ),
        "convergence_check_floors": CONVERGENCE_CHECK_FLOORS,
        "no_new_free_parameter": (
            "the floor is an integration bound driven to convergence, not a "
            "fitted quantity.  The convergence scan is reported in full so the "
            "choice is auditable."
        ),
        "predictions": [
            {
                "id": "F1",
                "statement": (
                    "The corrected closure ratio falls in ALL 54 cells, because "
                    "up-scatter can only add to the predicted observable count."
                ),
                "falsifies_if": "any cell's ratio rises or is unchanged",
                "if_falsified": (
                    "the response or the integration is doing something other "
                    "than what is claimed, and the fix must be withdrawn pending "
                    "diagnosis"
                ),
            },
            {
                "id": "F2",
                "statement": (
                    "The reduction is LARGER at steeper alpha, because a steeper "
                    "IMF puts more weight at the low-mass end where the "
                    "previously-omitted up-scatter lives.  Ordering: effect at "
                    "alpha 2.6 > 2.3 > 2.0."
                ),
                "falsifies_if": "the ordering of the mean reduction is not 2.6 > 2.3 > 2.0",
                "if_falsified": (
                    "the mechanism is not the IMF-weighted up-scatter it is "
                    "claimed to be"
                ),
            },
            {
                "id": "F3",
                "statement": (
                    "The effect is approximately subgroup-INDEPENDENT.  It is a "
                    "threshold effect at 8 Msun, so it should NOT follow the "
                    "turnoff ordering A < B < C that governs the upper limit.  "
                    "Quantitatively: the spread in fractional reduction across "
                    "the three subgroups at the baseline branch is under 5 "
                    "percentage points."
                ),
                "falsifies_if": (
                    "the fractional reduction spans more than 5 percentage "
                    "points across subgroups"
                ),
                "if_falsified": (
                    "the lower limit and the turnoff are entangled in a way not "
                    "understood, and the attribution in wp6_closure.md section 4 "
                    "must be redone"
                ),
                "note": (
                    "this is the deliberate contrast with issue #15's prediction "
                    "M2, which requires the multiplicity effect to FOLLOW the "
                    "turnoff ordering.  The two mechanisms are separable exactly "
                    "because they predict different subgroup patterns."
                ),
            },
            {
                "id": "F4",
                "statement": (
                    "Decisive for the interpretation: the grid-median closure "
                    "ratio at alpha = 2.3 falls from 1.444 to below 1.25, and "
                    "the grid-median closing alpha rises from 2.070 toward "
                    "Salpeter."
                ),
                "falsifies_if": "the grid-median ratio stays at or above 1.25",
                "if_falsified": (
                    "the single-node probe was unrepresentative and the lower "
                    "limit is a real but minor bug; the IMF reading then stands "
                    "on the strength it had before"
                ),
            },
        ],
        "decision_rule": (
            "The lower-limit fix is a BUG FIX, not a branch: if F1 holds the "
            "corrected estimator replaces the incumbent unconditionally, because "
            "the incumbent computes a quantity that does not match its own "
            "observed side.  It is adopted even if it makes the closure ratio "
            "worse.  F2/F3/F4 test whether the mechanism is understood, not "
            "whether to adopt it.  The published closure ratios (1.087/1.448/"
            "1.706, closing alpha 2.070) are WITHDRAWN on adoption and every "
            "downstream document must be updated."
        ),
        "interaction_with_issue_15": (
            "these are independent defects and must not be conflated.  Issue #15 "
            "changes the TRUTH model above 8 Msun; issue #17 changes the "
            "INTEGRATION BOUND below it.  The issue #15 diagnostic was scoped to "
            "M >= 8 deliberately, which — as issue #17 now makes clear — is the "
            "range where the multiplicity mechanism has the LEAST room to act, "
            "since a star already above 8 Msun cannot be scattered into the "
            "census by being made brighter.  That scope limit must be stated "
            "when the #15 result is reported."
        ),
        "consequence_for_the_15_prereg": (
            "wp6_multiplicity_prereg.json declared M3 as an ABSOLUTE threshold, "
            "1.222, arithmetically derived from the then-published 1.444 as "
            "'absorbs at least HALF of the baseline excess'.  Issue #17 "
            "invalidates that reference state.  The prereg is NOT amended — "
            "amending a pre-registration after seeing data is precisely what it "
            "exists to prevent.  Instead M3 is scored BOTH ways and both are "
            "reported: literally against 1.222, and in the relative form its own "
            "sentence states ('absorbs at least half of the excess'), which is "
            "invariant to the baseline correction.  The relative form is the "
            "faithful reading of what was written down."
        ),
        "reference_state_being_withdrawn": {
            "closure_ratio_alpha2.3": [
                {"subgroup": row.subgroup, "closure_ratio": round(row.closure_ratio, 3)}
                for row in baseline.itertuples()
            ],
            "grid_median_alpha2.3": round(
                float(closure[closure.alpha.eq(2.3)].closure_ratio.median()), 3
            ),
            "closing_alpha_grid_median": 2.070,
        },
        "scope": {
            "no_new_injections": (
                "the frozen MASS_GRID already spans 0.5-18 Msun and the accepted "
                "repair_v6 node responses already cover 4-8 Msun.  Only the "
                "integration bound changes."
            ),
            "wp5_untouched": (
                "the WP5 normalization k is fitted on 2-8 Msun COUNTS and does "
                "not involve this integral at all, so no accepted WP5 artifact "
                "moves.  Verified by re-running V1 after the fix."
            ),
        },
    }
    w.write_json(w.PROVENANCE / "wp6_closure_floor_prereg.json", record)

    print("issue #17 — closure-integral lower limit, pre-registered\n")
    print("  defect: predicted side integrates from 8 Msun; observed side")
    print("          counts P(M > 8) for members of ANY true mass")
    print(f"\n  adopted floor: {ADOPTED_FLOOR_MSUN} Msun (chosen on convergence)")
    print("  probe node, floor 8.0 -> 4.0: +21.7% predicted observable")
    print("\n  pre-declared:")
    for entry in record["predictions"]:
        print(f"    {entry['id']}: {entry['statement'][:66]}...")
    print("\nwrote provenance/wp6_closure_floor_prereg.json")


if __name__ == "__main__":
    main()
