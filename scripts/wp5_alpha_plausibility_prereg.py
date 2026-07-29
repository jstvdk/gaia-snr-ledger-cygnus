#!/usr/bin/env python3
"""Pre-registration: how plausible is each IMF slope branch?

Written BEFORE the statistics are computed.  The statistics and the
prohibitions are fixed here.

THIS IS A DOCUMENTATION EXERCISE, NOT A DECISION
------------------------------------------------
WP7 reports N_SN = 1.93-28.74 across 54 all-explode branches, and a variance
decomposition attributes 91.9% of that spread to the IMF slope alpha alone
(R_V 4.0%, isochrone family 2.0%, star-formation duration 0.2%).  So the
question "can the range be narrowed" is really one question about alpha.

Two bodies of evidence already exist in the project and neither was designed as
an alpha measurement.  This script measures what they say.  It does NOT act on
the answer.

BINDING PROHIBITIONS
--------------------
  * NO branch is removed.  All 54 remain in tables/wp7_ledger.csv.
  * NO artifact is regenerated, refitted or reweighted.
  * NO branch weights are applied to any reported number anywhere.
  * The WP7, WP6 and WP5 results stand exactly as published.

The output is a plausibility statement for the principal investigator to act on
later, or not.  If it is ever acted on, that adoption is a separate decision
requiring its own pre-registration.

THE TWO LINES OF EVIDENCE
-------------------------
E1  the 2-8 Msun CALIBRATION WINDOW, internal to WP5.
    k is refitted at each alpha, so the Poisson chi-square of that fit measures
    whether the assumed slope matches the observed low-mass mass function.
    This costs NOTHING: it is a gate statistic WP5 already computes, and using
    it consumes no validation.

E2  the >8 Msun CENSUS CLOSURE, from WP6.
    k is fitted from 2-8 Msun counts alone and the >8 Msun census never enters
    the WP5 likelihood, so the closure ratio is an out-of-sample test of the
    extrapolation.  It is a far stronger discriminant than E1 -- and using it
    to select alpha has a REAL COST, stated below.

E1 and E2 probe DIFFERENT MASS RANGES, which is what makes their agreement (or
disagreement) informative.

Output: provenance/wp5_alpha_plausibility_prereg.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp5_alpha_plausibility_prereg.py
"""
from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w

CHI_SQUARE_DOF = 5


def main() -> None:
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_alpha_plausibility_prereg.py",
        "status": "PREREGISTERED",
        "work_package": "WP5/WP6/WP7 cross-cutting diagnostic",
        "question": (
            "how plausible is each of the three carried IMF slope branches, "
            "given evidence the project already has?"
        ),
        "motivation": {
            "variance_decomposition_of_log_N_SN": {
                "alpha": 0.919, "R_V": 0.040, "family": 0.020,
                "sf_duration": 0.002,
            },
            "consequence": (
                "the WP7 spread is one question about alpha, not four "
                "questions about four branch axes.  Narrowing R_V or the "
                "isochrone family would buy almost nothing."
            ),
        },
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "binding_prohibitions": [
            "NO branch is removed -- all 54 remain in tables/wp7_ledger.csv",
            "NO artifact is regenerated, refitted or reweighted",
            "NO branch weights are applied to any reported number anywhere",
            "the WP5, WP6 and WP7 results stand exactly as published",
            "this is a plausibility statement for a later decision, and any "
            "such decision requires its own separate pre-registration",
        ],
        "evidence_E1_calibration_window": {
            "range_Msun": [2.0, 8.0],
            "statistic": (
                "the Poisson chi-square of the WP5 fit, per cell, with "
                f"{CHI_SQUARE_DOF} degrees of freedom, already recorded in "
                "wp5_imf_normalization_repair_v7.parquet"
            ),
            "why_it_is_a_constraint_on_alpha": (
                "k is refitted at every alpha, so the NORMALIZATION is free and "
                "cannot absorb a slope error.  What the chi-square measures is "
                "whether the assumed SHAPE matches the observed mass function "
                "inside the window."
            ),
            "cost_of_using_it": (
                "none.  It is internal to WP5 and is already a published gate "
                "statistic (residual_gate_pass).  Reading it consumes no "
                "validation."
            ),
            "reported_as": [
                "per-cell chi-square for each alpha, all 18 cells",
                "which alpha wins in each cell, and the tally",
                "median delta chi-square relative to the best alpha",
                "per-subgroup breakdown, because the subgroups may disagree",
            ],
        },
        "evidence_E2_census_closure": {
            "range_Msun": [8.0, 120.0],
            "statistic": (
                "the closure ratio and its 68% interval per cell, from "
                "tables/wp6_closure_repair_v7.csv, expressed as the number of "
                "68% half-widths between the ratio and unity"
            ),
            "why_it_is_out_of_sample": (
                "k is fitted from 2-8 Msun counts alone; the >8 Msun census "
                "never enters the WP5 likelihood.  Testing the extrapolation "
                "against it is therefore a genuine out-of-sample test, which is "
                "what the branch grid was carried in order to receive."
            ),
            "THE_COST_OF_USING_IT": (
                "the census closure is currently this analysis's ONLY "
                "out-of-sample validation.  If it is used to SELECT alpha, it "
                "can no longer also be reported as an independent confirmation "
                "of the IMF -- the project would get one or the other, not "
                "both.  PROJECT_TRACE section 10c item 2 already records the "
                "standing decision that WP7 must not refit alpha for exactly "
                "this reason.  This diagnostic MEASURES what E2 says; it does "
                "not spend it."
            ),
        },
        "no_formal_combination": {
            "rule": (
                "E1 and E2 are reported SEPARATELY and are never multiplied "
                "into a joint likelihood or posterior over alpha."
            ),
            "why": (
                "the 18 cells are not independent -- they share the same stars "
                "across families and R_V values -- so multiplying per-cell "
                "likelihoods would manufacture significance.  E1 and E2 also "
                "share k, so they are not cleanly independent of one another "
                "either.  What is legitimate is to report each, note whether "
                "they agree, and let the agreement of two different mass ranges "
                "carry qualitative rather than numerical weight."
            ),
        },
        "predictions": [
            {
                "id": "A1",
                "statement": (
                    "E1 and E2 agree in DIRECTION: both disfavour alpha = 2.6 "
                    "relative to alpha = 2.3"
                ),
                "falsifies_if": (
                    "either line favours 2.6 over 2.3 on its own statistic"
                ),
                "if_falsified": (
                    "the two mass ranges disagree about the slope, which is "
                    "itself a reportable result and would mean the single "
                    "power-law assumption is the thing under strain"
                ),
            },
            {
                "id": "A2",
                "statement": (
                    "the subgroups DISAGREE about alpha, with CygOB2-C "
                    "preferring a shallower slope than CygOB2-A on both lines "
                    "of evidence"
                ),
                "falsifies_if": (
                    "C's preferred alpha is greater than or equal to A's on "
                    "either line"
                ),
                "if_falsified": (
                    "the 'C is different' signal that WP6's M2 and the closing "
                    "alpha both reported would not be reproducing here, and its "
                    "earlier appearances would need re-examination"
                ),
            },
            {
                "id": "A3",
                "statement": (
                    "E2 discriminates more sharply than E1, because the >8 Msun "
                    "census sits far from the calibration window and a slope "
                    "error compounds over that lever arm"
                ),
                "falsifies_if": (
                    "E1's separation between alpha branches exceeds E2's, "
                    "measured in each statistic's own units of uncertainty"
                ),
                "if_falsified": (
                    "the extrapolation is less sensitive to alpha than the "
                    "window fit, which would weaken the case for treating the "
                    "closure test as the primary IMF discriminant"
                ),
            },
        ],
        "what_this_diagnostic_cannot_do": [
            "it cannot measure alpha.  The project carries three values and has "
            "never run a branch between them; a plausibility ordering over three "
            "grid points is not a fitted slope with an interval",
            "it cannot justify narrowing N_SN unless the IMF is a SINGLE POWER "
            "LAW from 2 to 120 Msun.  E1 constrains the slope at 2-8 Msun and "
            "N_SN depends on the slope above 52 Msun; the link between them is "
            "the single-power-law assumption, which this diagnostic assumes and "
            "does not test",
            "it cannot resolve the subgroup disagreement, only measure it",
            "it cannot be used to reweight any published number without a "
            "separate pre-registered adoption decision",
        ],
        "deliverable": (
            "a plausibility table per alpha, the implied N_SN range under each "
            "candidate branch set, and an explicit statement of what adopting "
            "each set would cost.  The principal investigator decides later, or "
            "not at all."
        ),
    }
    w.write_json(w.PROVENANCE / "wp5_alpha_plausibility_prereg.json", record)

    print("alpha plausibility diagnostic — pre-registered\n")
    print("  91.9% of the WP7 N_SN spread is the alpha axis alone")
    print("\n  E1  2-8 Msun calibration-window chi-square   (free to use)")
    print("  E2  >8 Msun census closure, out-of-sample    (costs the validation)")
    print("\n  BINDING: no branch removed, no artifact regenerated,")
    print("           no weights applied.  Documentation only.")
    print("\n  predictions:")
    for entry in record["predictions"]:
        print(f"    {entry['id']}  {entry['statement'][:62]}")
    print("\nwrote provenance/wp5_alpha_plausibility_prereg.json")


if __name__ == "__main__":
    main()
