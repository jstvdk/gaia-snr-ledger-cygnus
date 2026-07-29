#!/usr/bin/env python3
"""How plausible is each IMF slope branch?  Measured, not acted on.

Pre-registered in provenance/wp5_alpha_plausibility_prereg.json.  Predictions
A1-A3 are scored, not amended.

NOTHING IS REMOVED, REGENERATED OR REWEIGHTED.  All 54 branches remain in
tables/wp7_ledger.csv and every published number stands.  This produces a
plausibility table and the implied N_SN under each CANDIDATE branch set, so the
principal investigator can decide later -- or not.

Outputs:
  tables/wp5_alpha_plausibility.csv       per-cell evidence, both lines
  tables/wp5_alpha_candidate_sets.csv     implied N_SN under each candidate set
  provenance/wp5_alpha_plausibility_execution.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp5_alpha_plausibility.py
"""
from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
from wp5_alpha_plausibility_prereg import CHI_SQUARE_DOF

WP5_VERSION = "repair_v7"


def main() -> None:
    norm = pd.read_parquet(w.PROC / f"wp5_imf_normalization_{WP5_VERSION}.parquet")
    closure = pd.read_csv(w.TABLES / f"wp6_closure_{WP5_VERSION}.csv")
    ledger = pd.read_csv(w.TABLES / "wp7_ledger.csv")

    # ---- E1: calibration-window chi-square, per cell -----------------------
    e1 = norm.pivot_table(
        index=["subgroup", "family", "R_V"], columns="alpha",
        values="poisson_chi_square",
    )
    e1_best = e1.idxmin(axis=1)
    e1_tally = e1_best.value_counts().sort_index()

    # ---- E2: closure ratio distance from unity, in 68% half-widths ---------
    closure = closure.copy()
    half = 0.5 * (closure.closure_ratio_hi68 - closure.closure_ratio_lo68)
    closure["unity_distance_sigma"] = (closure.closure_ratio - 1.0).abs() / half
    e2 = closure.pivot_table(
        index=["subgroup", "family", "R_V"], columns="alpha",
        values="unity_distance_sigma",
    )
    e2_best = e2.idxmin(axis=1)
    e2_tally = e2_best.value_counts().sort_index()

    rows = []
    for key in e1.index:
        subgroup, family, rv = key
        for alpha in sorted(e1.columns):
            rows.append(
                {
                    "subgroup": subgroup, "family": family, "R_V": rv,
                    "alpha": alpha,
                    "E1_chi_square": round(float(e1.loc[key, alpha]), 3),
                    "E1_delta_vs_best": round(
                        float(e1.loc[key, alpha] - e1.loc[key].min()), 3
                    ),
                    "E1_is_best": bool(e1_best.loc[key] == alpha),
                    "E2_closure_ratio": round(
                        float(
                            closure[
                                closure.subgroup.eq(subgroup)
                                & closure.family.eq(family)
                                & closure.R_V.eq(rv) & closure.alpha.eq(alpha)
                            ].closure_ratio.iloc[0]
                        ), 3
                    ),
                    "E2_unity_distance_sigma": round(
                        float(e2.loc[key, alpha]), 2
                    ),
                    "E2_is_best": bool(e2_best.loc[key] == alpha),
                }
            )
    table = pd.DataFrame(rows)
    table.to_csv(w.TABLES / "wp5_alpha_plausibility.csv", index=False)

    # ---- per-subgroup preference, for prediction A2 ------------------------
    per_subgroup = {}
    for subgroup in w.SUBGROUPS:
        block = table[table.subgroup.eq(subgroup)]
        e1_pref = float(
            block.groupby("alpha").E1_chi_square.median().idxmin()
        )
        e2_pref = float(
            block.groupby("alpha").E2_unity_distance_sigma.median().idxmin()
        )
        per_subgroup[subgroup] = {
            "E1_preferred_alpha": e1_pref,
            "E2_preferred_alpha": e2_pref,
            "E1_median_chi_square_by_alpha": {
                f"{a:g}": round(float(v), 2)
                for a, v in block.groupby("alpha").E1_chi_square.median().items()
            },
            "E2_median_sigma_by_alpha": {
                f"{a:g}": round(float(v), 2)
                for a, v in
                block.groupby("alpha").E2_unity_distance_sigma.median().items()
            },
        }

    # ---- candidate branch sets: what each WOULD imply, if ever adopted ----
    assoc = ledger[
        ledger.scope.eq("association") & ledger.explodability.eq("all_explode")
    ]
    candidates = []
    for name, mask, cost in (
        (
            "as published (all 54 branches)",
            assoc.alpha.notna(),
            "none -- this is the current reported range",
        ),
        (
            "drop alpha=2.6",
            assoc.alpha.ne(2.6),
            "E1 only, which is free; E2 agrees but need not be spent",
        ),
        (
            "alpha=2.3 only",
            assoc.alpha.eq(2.3),
            "spends E2 -- the census could no longer be reported as an "
            "independent confirmation of the IMF",
        ),
        (
            "alpha in {2.0, 2.3}, drop R_V=3.5",
            assoc.alpha.ne(2.6) & assoc.R_V.ne(3.5),
            "requires independent extinction evidence for R_V, which this "
            "diagnostic does NOT provide",
        ),
    ):
        block = assoc[mask]
        candidates.append(
            {
                "candidate_set": name,
                "n_branches": int(len(block)),
                "N_SN_min": round(float(block.N_SN_mean.min()), 2),
                "N_SN_max": round(float(block.N_SN_mean.max()), 2),
                "N_SN_median": round(float(block.N_SN_mean.median()), 2),
                "factor": round(
                    float(block.N_SN_mean.max() / block.N_SN_mean.min()), 2
                ),
                "evidential_cost": cost,
            }
        )
    pd.DataFrame(candidates).to_csv(
        w.TABLES / "wp5_alpha_candidate_sets.csv", index=False
    )

    # ---- score the pre-registered predictions ------------------------------
    e1_med = table.groupby("alpha").E1_chi_square.median()
    e2_med = table.groupby("alpha").E2_unity_distance_sigma.median()
    a1 = bool(e1_med[2.6] > e1_med[2.3] and e2_med[2.6] > e2_med[2.3])
    a2 = bool(
        per_subgroup["CygOB2-C"]["E1_preferred_alpha"]
        < per_subgroup["CygOB2-A"]["E1_preferred_alpha"]
        and per_subgroup["CygOB2-C"]["E2_preferred_alpha"]
        < per_subgroup["CygOB2-A"]["E2_preferred_alpha"]
    )
    # A3: separation between the best and worst alpha, in each statistic's own
    # units of uncertainty.  chi-square differences are in units of its own
    # sqrt(2*dof) sampling spread; the closure statistic is already in sigma.
    e1_separation = float(
        (e1_med.max() - e1_med.min()) / np.sqrt(2.0 * CHI_SQUARE_DOF)
    )
    e2_separation = float(e2_med.max() - e2_med.min())
    a3 = bool(e2_separation > e1_separation)

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_alpha_plausibility.py",
        "status": "SUCCESS",
        "work_package": "WP5/WP6/WP7 cross-cutting diagnostic",
        "prereg": "provenance/wp5_alpha_plausibility_prereg.json",
        "nothing_was_changed": (
            "no branch removed, no artifact regenerated, no weights applied.  "
            "All 54 branches remain in tables/wp7_ledger.csv and every "
            "published number stands unaltered.  This record is documentation "
            "for a later decision that has NOT been taken."
        ),
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "E1_calibration_window": {
            "cells": int(len(e1)),
            "median_chi_square_by_alpha": {
                f"{a:g}": round(float(v), 2) for a, v in e1_med.items()
            },
            "dof": CHI_SQUARE_DOF,
            "wins_by_alpha": {
                f"{a:g}": int(n) for a, n in e1_tally.items()
            },
            "reading": (
                "alpha = 2.6 is the best fit in only "
                f"{int(e1_tally.get(2.6, 0))} of {len(e1)} cells and is the "
                "worst on the median.  alpha = 2.0 and 2.3 are genuinely "
                "competitive with each other."
            ),
        },
        "E2_census_closure": {
            "cells": int(len(e2)),
            "median_unity_distance_sigma_by_alpha": {
                f"{a:g}": round(float(v), 2) for a, v in e2_med.items()
            },
            "wins_by_alpha": {
                f"{a:g}": int(n) for a, n in e2_tally.items()
            },
            "reading": (
                "the closure ratio is furthest from unity at alpha = 2.6 and "
                "closest at alpha = 2.3.  This is the out-of-sample line and it "
                "has NOT been spent -- it is reported here, not used."
            ),
        },
        "per_subgroup": per_subgroup,
        "predictions": [
            {
                "id": "A1",
                "statement": "E1 and E2 agree in direction, both disfavouring 2.6",
                "outcome": "PASS" if a1 else "FAIL",
                "E1_median_2p3_vs_2p6": [
                    round(float(e1_med[2.3]), 2), round(float(e1_med[2.6]), 2)
                ],
                "E2_median_2p3_vs_2p6": [
                    round(float(e2_med[2.3]), 2), round(float(e2_med[2.6]), 2)
                ],
            },
            {
                "id": "A2",
                "statement": "the subgroups disagree, C preferring a shallower slope than A",
                "outcome": "PASS" if a2 else "FAIL",
                "C_preferred": [
                    per_subgroup["CygOB2-C"]["E1_preferred_alpha"],
                    per_subgroup["CygOB2-C"]["E2_preferred_alpha"],
                ],
                "A_preferred": [
                    per_subgroup["CygOB2-A"]["E1_preferred_alpha"],
                    per_subgroup["CygOB2-A"]["E2_preferred_alpha"],
                ],
            },
            {
                "id": "A3",
                "statement": "E2 discriminates more sharply than E1",
                "outcome": "PASS" if a3 else "FAIL",
                "E1_separation_in_own_units": round(e1_separation, 2),
                "E2_separation_in_sigma": round(e2_separation, 2),
            },
        ],
        "candidate_branch_sets": candidates,
        "what_the_PI_is_being_asked_to_decide": {
            "the_question": (
                "whether to report N_SN over all three alpha branches, or over "
                "a subset justified by the evidence measured here"
            ),
            "the_cheapest_defensible_move": (
                "dropping alpha = 2.6.  E1 alone supports it, E1 is internal to "
                "WP5 and free, and E2 independently agrees without needing to "
                "be spent.  This cuts the reported spread from a factor of "
                f"{candidates[0]['factor']} to {candidates[1]['factor']}."
            ),
            "the_move_with_a_price": (
                "restricting to alpha = 2.3.  It cuts the spread to a factor of "
                f"{candidates[2]['factor']}, but E1 alone does not support it -- "
                "alpha = 2.0 wins in a substantial minority of cells -- so it "
                "requires E2, and spending E2 costs the project its only "
                "out-of-sample validation of the IMF."
            ),
            "the_move_not_supported_here": (
                "dropping R_V = 3.5.  This diagnostic provides no extinction "
                "evidence and the R_V axis carries only 4% of the variance.  It "
                "is listed for completeness, not recommended."
            ),
            "what_no_choice_can_fix": (
                "CygOB2-C contributes 0 or about 7 supernovae depending on "
                "isochrone family, because its fitted age straddles the age at "
                "which the turnoff crosses the 120 Msun IMF ceiling.  That is "
                "independent of alpha and survives every candidate set."
            ),
        },
        "outputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in [
                w.TABLES / "wp5_alpha_plausibility.csv",
                w.TABLES / "wp5_alpha_candidate_sets.csv",
            ]
        },
    }
    w.write_json(
        w.PROVENANCE / "wp5_alpha_plausibility_execution.json", record
    )

    print("alpha plausibility — measured, nothing acted on\n")
    print("  E1  2-8 Msun calibration window (free)")
    print(f"     median chi2 (dof {CHI_SQUARE_DOF}): " + "  ".join(
        f"a={a:g}: {v:5.2f}" for a, v in e1_med.items()))
    print("     wins: " + "  ".join(
        f"a={a:g}: {n}/18" for a, n in e1_tally.items()))
    print("\n  E2  >8 Msun census closure (out-of-sample, NOT spent)")
    print("     median |ratio-1| in 68% half-widths: " + "  ".join(
        f"a={a:g}: {v:5.2f}" for a, v in e2_med.items()))
    print("     wins: " + "  ".join(
        f"a={a:g}: {n}/18" for a, n in e2_tally.items()))
    print("\n  per subgroup, preferred alpha (E1 / E2):")
    for subgroup, entry in per_subgroup.items():
        print(f"     {subgroup:10s} {entry['E1_preferred_alpha']:.1f} / "
              f"{entry['E2_preferred_alpha']:.1f}")
    print("\n  candidate branch sets (NOT adopted):")
    for entry in candidates:
        print(f"     {entry['candidate_set']:34s} "
              f"{entry['N_SN_min']:6.2f}-{entry['N_SN_max']:6.2f}  "
              f"factor {entry['factor']:5.2f}")
    print("\n  predictions:")
    for entry in record["predictions"]:
        print(f"     {entry['id']}  {entry['outcome']}")
    print("\nwrote provenance/wp5_alpha_plausibility_execution.json")


if __name__ == "__main__":
    main()
