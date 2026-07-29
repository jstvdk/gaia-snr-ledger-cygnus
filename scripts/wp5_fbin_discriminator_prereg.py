#!/usr/bin/env python3
"""Pre-registration: does multiplicity BELOW 8 Msun justify a repair_v7?

Written BEFORE any injection is generated.  Design fixed in
tasks/repair_v7_recommendation.md section 4 and not altered since.

THE QUESTION
------------
Issue #15 measured mass-dependent multiplicity ABOVE 8 Msun and found it absorbs
3.7% of the WP6 closure excess -- so it is not the explanation.  But that test
was scoped to the range where the mechanism has LEAST room to act: a star
already above 8 Msun cannot be scattered INTO the census by being made brighter.

Below 8 Msun the truth model still applies a constant f_bin = 0.40, and that is
where two things live:

  the WP5 calibration window (2-8 Msun), which sets the normalization k
  the up-scatter channel (4-8 Msun), worth ~24% of the predicted observable
  count above 8 Msun after issue #17

A full repair_v7 costs ~8 h plus G3 re-acceptance and would move an ACCEPTED
artifact.  This measures the size of the effect first, on three nodes, for a
couple of minutes of compute.  It is a go/no-go instrument, not a result.

THE MODEL UNDER TEST -- DELIBERATELY CONSERVATIVE
-------------------------------------------------
    f_bin = 0.40 at 2 Msun  ->  0.55 at 8 Msun  ->  0.70 at 16 Msun

log-linear between anchors, flat outside.  Duchene & Kraus 2013 give ~50%
multiplicity already at late-B/A masses (2-5 Msun), so holding 0.40 at 2 Msun
UNDERSTATES the likely change.  The measured effect is therefore a LOWER BOUND,
which is the right direction for a go/no-go test: if even this conservative
model clears the threshold, repair_v7 is justified a fortiori.

Output: provenance/wp5_fbin_discriminator_prereg.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp5_fbin_discriminator_prereg.py
"""
from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w

# The extended truth-side model.  Below-8 anchors are the new part; the 16 Msun
# anchor is inherited unchanged from issue #15 so the two tests are commensurate.
FBIN_ANCHORS_MSUN = np.array([2.0, 8.0, 16.0])
FBIN_ANCHORS = np.array([0.40, 0.55, 0.70])

# Windows the two questions live in.
CALIBRATION_WINDOW = (2.0, 8.0)   # sets k -- the WP5 question
UPSCATTER_WINDOW = (4.0, 8.0)     # feeds the >8 census -- the WP6 question
WEIGHT_ALPHA = 2.3                # IMF weighting for both means, baseline branch

# Pre-declared go/no-go threshold.  Chosen to match a systematic already carried
# in the project (the +2.4% Orellana distance offset), so it is calibrated
# against something real rather than invented for this test.
DECISION_THRESHOLD = 0.02

BASELINE_FAMILY, BASELINE_RV = "PARSEC", 3.1


def extended_binary_fraction(mass: np.ndarray) -> np.ndarray:
    """f_bin(M) with the below-8 Msun rise switched on."""
    mass = np.asarray(mass, dtype=float)
    return np.interp(
        np.log10(mass), np.log10(FBIN_ANCHORS_MSUN), FBIN_ANCHORS,
        left=FBIN_ANCHORS[0], right=FBIN_ANCHORS[-1],
    )


def main() -> None:
    grid = w.MASS_GRID
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_fbin_discriminator_prereg.py",
        "status": "PREREGISTERED",
        "question": (
            "Does mass-dependent multiplicity BELOW 8 Msun move the WP5 "
            "calibration or the 4-8 Msun up-scatter channel enough to justify a "
            "repair_v7 full-chain re-run?"
        ),
        "why_issue_15_does_not_answer_it": (
            "issue #15 was scoped to M >= 8 Msun, the range where the mechanism "
            "has least room to act, because a star already above the threshold "
            "cannot be scattered into the census by being brightened.  It "
            "measured 3.7% there.  Below 8 Msun the mechanism acts through two "
            "channels it could not touch: the normalization k, and up-scatter "
            "across the 8 Msun threshold."
        ),
        "instrument_not_result": (
            "this is a go/no-go discriminator on three nodes, not a measurement "
            "of the effect.  A positive outcome authorizes repair_v7; it does "
            "not itself change any number."
        ),
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "design": {
            "branch": f"{BASELINE_FAMILY} R_V={BASELINE_RV} (reporting branch)",
            "subgroups": list(w.SUBGROUPS),
            "node_per_subgroup": (
                "the single highest-prior-weight truth-age node, so the test "
                "sits where the age posterior actually concentrates"
            ),
            "mass_grid": "the frozen wp5_common.MASS_GRID, 0.5-18 Msun",
            "arms": {
                "ctl": "frozen constant F_BINARY = 0.40 at every mass",
                "trt": "extended f_bin(M) with the below-8 rise switched on",
            },
            "pairing": (
                "both arms run from a fresh default_rng(SEED) on an identical "
                "grid, and the per-star binary threshold consumes the same "
                "single rng.random(n) draw, so donor, extinction, photometric "
                "and QMC realizations are bit-identical.  The arms differ in "
                "exactly one thing: which stars got a companion.  Realization "
                "noise therefore cancels in the difference, which is why three "
                "nodes suffice."
            ),
            "recovery_side_unchanged": (
                "as in issue #15, the estimator keeps assuming 0.40.  The "
                "mismatch between the rate nature makes binaries at and the "
                "rate the estimator assumes IS the bias under test."
            ),
        },
        "model": {
            "anchors_Msun": FBIN_ANCHORS_MSUN.tolist(),
            "anchors_f_bin": FBIN_ANCHORS.tolist(),
            "interpolation": "linear in log mass, flat outside the anchors",
            "grid": {
                f"{m:g}": round(float(extended_binary_fraction(m)), 3)
                for m in [0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 12.0, 16.0, 18.0]
            },
            "deliberately_conservative": (
                "Duchene & Kraus 2013 give ~50% multiplicity already at "
                "late-B/A masses (2-5 Msun), so holding 0.40 at the 2 Msun "
                "anchor UNDERSTATES the likely change.  The measured effect is "
                "a LOWER BOUND, which is the correct direction for a go/no-go "
                "test: if even this model clears the threshold, repair_v7 is "
                "justified a fortiori."
            ),
            "relation_to_issue_15": (
                "the 16 Msun anchor of 0.70 is inherited unchanged from "
                "wp6_multiplicity_prereg.json, so the two tests are "
                "commensurate and differ only below 8 Msun"
            ),
            "relation_to_CUTS_row_20": (
                "row 20 carries f_bin as a Class E branch spanning 0.3-0.5.  "
                "Both 0.55 and 0.70 lie outside that range, which is the "
                "defect already recorded in CUTS section 18.0."
            ),
        },
        "measurements": {
            "D1": {
                "quantity": (
                    "IMF-weighted mean of R(recovered | M) over the 2-8 Msun "
                    f"calibration window, weight proportional to M^-{WEIGHT_ALPHA}"
                ),
                "window_Msun": list(CALIBRATION_WINDOW),
                "answers": (
                    "does the WP5 normalization k move?  k is fitted so that "
                    "k * integral M^-alpha R dM reproduces the observed counts "
                    "in this window, so a fractional shift in this response is "
                    "the fractional shift k must absorb."
                ),
            },
            "D2": {
                "quantity": (
                    "IMF-weighted mean of R(estimated > 8 Msun | M) over the "
                    "4-8 Msun window"
                ),
                "window_Msun": list(UPSCATTER_WINDOW),
                "answers": (
                    "does the up-scatter channel move?  Issue #17 showed this "
                    "range contributes ~24% of the predicted observable count "
                    "above 8 Msun, and it is where unresolved companions matter "
                    "most because these are the stars that can be pushed ACROSS "
                    "the threshold."
                ),
            },
        },
        "decision_rule": {
            "threshold": DECISION_THRESHOLD,
            "rule": (
                "if EITHER D1 or D2 shifts by more than 2% in absolute "
                "fractional terms, averaged over the three nodes, repair_v7 is "
                "JUSTIFIED and should be scheduled.  If both stay below 2%, the "
                "effect is recorded as a carried systematic and the accepted "
                "chain stands unchanged."
            ),
            "threshold_justification": (
                "2% is not invented for this test: it is the size of the "
                "Orellana distance systematic already carried in the project "
                "(+2.4% on the closure ratios).  An effect smaller than a "
                "systematic we already accept does not warrant re-running an "
                "accepted chain; an effect larger than it does."
            ),
            "either_not_both": (
                "D1 and D2 answer different questions -- k versus up-scatter -- "
                "and either alone would matter, so the rule is disjunctive.  "
                "Requiring both would let a large effect in one channel be "
                "excused by a small one in the other."
            ),
        },
        "predictions": [
            {
                "id": "G1",
                "statement": (
                    "Both D1 and D2 shift in the POSITIVE direction: more "
                    "companions means more flux at fixed true mass, so recovery "
                    "rises and up-scatter across 8 Msun rises."
                ),
                "falsifies_if": "either shifts negative by more than its node-to-node spread",
                "if_falsified": (
                    "the mechanism as reasoned is wrong and the discriminator's "
                    "sign convention must be re-derived before its magnitude is "
                    "trusted"
                ),
            },
            {
                "id": "G2",
                "statement": (
                    "D2 exceeds D1.  The up-scatter channel is a threshold "
                    "effect and should be more sensitive to a brightness shift "
                    "than the recovery fraction inside the calibration window, "
                    "where most stars are recovered either way."
                ),
                "falsifies_if": "D1 >= D2",
                "if_falsified": (
                    "the two channels are not separable the way this test "
                    "assumes, and a repair_v7 decision cannot be delegated to "
                    "the larger of the two"
                ),
            },
        ],
        "what_this_test_cannot_do": [
            "it cannot measure how much k would actually move -- only whether "
            "the response feeding k moves enough to be worth finding out",
            "it uses three nodes at one branch, so it says nothing about "
            "branch-to-branch or age-to-age variation",
            "a negative outcome does not show multiplicity is unimportant, only "
            "that it is smaller than a systematic already carried",
        ],
        "frozen_inputs": {
            "mass_grid_Msun": [float(m) for m in grid],
            "n_inject_per_mass": int(w.N_INJECT_PER_MASS),
            "seed": int(w.SEED),
            "F_BINARY_control": float(w.F_BINARY),
            "Q_MIN": float(w.Q_MIN),
        },
        "bit_preservation_obligation": (
            "the control arm must reproduce the accepted repair_v6 node "
            "response byte for byte, since it differs from it in nothing.  This "
            "is checked before the treatment arm is scored."
        ),
    }
    w.write_json(w.PROVENANCE / "wp5_fbin_discriminator_prereg.json", record)

    print("repair_v7 discriminator — pre-registered\n")
    print("  truth-side f_bin(M), below-8 rise switched on:")
    for mass, value in record["model"]["grid"].items():
        print(f"    {mass:>4s} Msun  {value:.3f}")
    print(f"\n  D1  R(recovered | M) over {CALIBRATION_WINDOW} Msun   -> does k move?")
    print(f"  D2  R(est > 8 | M)   over {UPSCATTER_WINDOW} Msun   -> does up-scatter move?")
    print(f"\n  decision: repair_v7 JUSTIFIED if either shifts > "
          f"{DECISION_THRESHOLD * 100:.0f}%")
    print("  (2% = the Orellana distance systematic already carried)")
    print("\nwrote provenance/wp5_fbin_discriminator_prereg.json")


if __name__ == "__main__":
    main()
