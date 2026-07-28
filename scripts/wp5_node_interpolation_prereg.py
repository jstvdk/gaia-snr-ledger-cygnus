#!/usr/bin/env python3
"""Issue #13, part c: pre-register the repair_v6 predictions before running.

Written and committed BEFORE any repair_v6 injection is generated, so the
outcome cannot be read back into the design.  Same discipline as
``provenance/wp3_kriging_prior_prereg.json``.

Output: provenance/wp5_node_interpolation_prereg.json
"""
from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
import wp5_joint_age_fit as J

BASELINE = ("PARSEC", 3.1, 2.3)


def main() -> None:
    age_posterior = pd.read_parquet(w.PROC / "wp4_age_posteriors_repair_v5.parquet")
    native = {family: J.native_isochrone_ages(family) for family in w.FAMILIES}

    node_plan = []
    for family in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            for subgroup in w.SUBGROUPS:
                snapped = J.truth_age_nodes(
                    age_posterior, subgroup, family, rv, native[family]
                )
                unsnapped = J.truth_age_nodes(
                    age_posterior, subgroup, family, rv, native[family], snap=False
                )
                prior_mean_snapped = sum(a * v for a, v in snapped.items())
                prior_mean_unsnapped = sum(a * v for a, v in unsnapped.items())
                node_plan.append(
                    {
                        "subgroup": subgroup,
                        "family": family,
                        "R_V": rv,
                        "n_nodes_repair_v5": len(snapped),
                        "n_nodes_repair_v6": len(unsnapped),
                        "prior_mean_age_repair_v5_Myr": round(prior_mean_snapped, 4),
                        "prior_mean_age_repair_v6_Myr": round(prior_mean_unsnapped, 4),
                        "prior_mean_age_shift_Myr": round(
                            prior_mean_unsnapped - prior_mean_snapped, 4
                        ),
                    }
                )
    plan = pd.DataFrame(node_plan)

    v5_bins = pd.read_parquet(w.PROC / "wp5_mass_function_bins_repair_v5.parquet")
    v5_mass = pd.read_parquet(w.PROC / "wp5_association_mass_repair_v5.parquet")
    regressing = v5_bins[
        v5_bins.subgroup.eq("CygOB2-C")
        & v5_bins.family.eq("MIST")
        & v5_bins.R_V.eq(3.5)
        & v5_bins.alpha.eq(2.0)
    ].sort_values("bin_index")

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_node_interpolation_prereg.py",
        "status": "PREREGISTERED",
        "issue": "#13 — truth-age node snapping replaced by isochrone interpolation",
        "version": "repair_v6",
        "upstream_unchanged": (
            "repair_v6 is a WP5-only version: WP3 extinction (kriged prior) and "
            "WP4 ages and masses are consumed from repair_v5 unchanged"
        ),
        "change": (
            "The nine WP4 age-posterior quantiles are no longer snapped to native "
            "isochrone ages.  Each keeps weight 1/9 at the age the posterior puts "
            "it, and wp5_common.load_isochrone_between_ages builds the truth "
            "isochrone there by the same bracketing and linear age blend the "
            "recovery side (wp4_repair_common._interpolate_age_sequence) already "
            "uses.  No new parameter: the node count is still N_AGE_NODES = 9."
        ),
        "motivation": (
            "Measured in provenance/wp5_node_rule_continuity_execution.json: the "
            "snapping rule amplifies a 0.0005 Myr shift of the WP4 posterior into "
            "a 0.055 Myr move of the truth-age node distribution (109x), and the "
            "node count varies between 1 and 4 within a single branch under a "
            "+-0.1 Myr translation.  The unsnapped rule is Lipschitz-1 — it moves "
            "in 1-Wasserstein distance by exactly the shift that caused it, which "
            "is the tightest continuity any rule can have."
        ),
        "not_tuned": (
            "The defect is identifiable without reference to any gate outcome: it "
            "is a discontinuity in a map from a continuous input, and it is an "
            "internal inconsistency between the truth and recovery sides of the "
            "same pipeline.  No gate statistic, threshold or window is changed. "
            "The adoption decision rests on P3 and P4 below, which are direction-"
            "free, and NOT on P2, which is the specific cell that motivated the "
            "investigation."
        ),
        "cost": {
            "nodes_repair_v5": int(plan["n_nodes_repair_v5"].sum()),
            "nodes_repair_v6": int(plan["n_nodes_repair_v6"].sum()),
            "reusable": 0,
            "note": (
                "every repair_v6 node is at a new age, so no repair_v4/v5 node "
                "response can be reused; repair_v1-v5 artifacts are untouched"
            ),
        },
        "predictions": [
            {
                "id": "P1",
                "statement": (
                    "The repair_v6 node rule is continuous in the WP4 posterior: "
                    "W1 jump per 0.0005 Myr offset step is 0.0005 Myr, ratio 1.0."
                ),
                "status": "ALREADY MEASURED — not a forecast",
                "evidence": "provenance/wp5_node_rule_continuity_execution.json",
            },
            {
                "id": "P2",
                "statement": (
                    "The one genuine repair_v4 -> repair_v5 regression, CygOB2-C "
                    "MIST R_V=3.5 alpha=2.0, clears: max|Pearson residual| < 3.0. "
                    "Its repair_v5 value of "
                    f"{float(np.max(np.abs(regressing.pearson_residual))):.2f} came "
                    "from deleting the 2.248 Myr node that held 41.6% of the "
                    "posterior weight; the unsnapped mixture spans 1.782-2.099 Myr "
                    "continuously."
                ),
                "falsifies_if": "max|r| stays >= 3.0",
                "if_falsified": (
                    "that cell's failure is NOT a node-set artifact.  It is then "
                    "a genuine branch failure, reported under the CUTS section 13 "
                    "retention policy, and issue #13's fix is judged on P3 and P4 "
                    "alone.  This prediction failing does not by itself reject the "
                    "fix — continuity is established independently by P1."
                ),
            },
            {
                "id": "P3",
                "statement": (
                    "The baseline (PARSEC, R_V=3.1, alpha=2.3) still passes the "
                    "residual gate for all three subgroups under both the "
                    "incumbent and the section 14 replacement trend statistic, and "
                    "CygOB2-B's |T| stays below 2."
                ),
                "falsifies_if": "any baseline subgroup fails either statistic",
                "if_falsified": (
                    "the finer node set is telling us repair_v5's baseline passed "
                    "partly by accident of the coarse mixture.  Issue #1c reopens "
                    "and the CygOB2-B result must be re-examined before anything "
                    "downstream is claimed."
                ),
            },
            {
                "id": "P4",
                "statement": (
                    "The fix is a local correction, not a wholesale re-scoring: "
                    "every branch's truth-age posterior mean moves by less than "
                    "0.15 Myr against repair_v5, the association stellar mass "
                    "stays within +-5%, and the 54-cell grid moves by at most 6 "
                    "cells in either direction."
                ),
                "falsifies_if": (
                    "any branch's posterior mean age moves >= 0.15 Myr, or the "
                    "association mass moves >= 5%, or |grid change| > 6"
                ),
                "if_falsified": (
                    "snapping was BIASING the truth model, not merely making it "
                    "jumpy.  That is a larger finding than issue #13 as written: "
                    "repair_v4's conclusions would inherit the same bias and would "
                    "have to be re-derived.  Report it as such rather than "
                    "adopting repair_v6 quietly."
                ),
                "direction_free": (
                    "P4 bounds the SIZE of the change, not its sign, so it cannot "
                    "be satisfied by moving the result toward passing"
                ),
            },
        ],
        "adoption_rule": (
            "repair_v6 is adopted if P1 (established), P3 and P4 hold.  P2 is "
            "reported either way and does not gate adoption.  If P3 or P4 fails "
            "the fix is not adopted, repair_v5 stands, and the failure is written "
            "up as a finding."
        ),
        "prior_mean_age_shift_forecast": {
            "max_abs_Myr": round(float(plan["prior_mean_age_shift_Myr"].abs().max()), 4),
            "median_abs_Myr": round(
                float(plan["prior_mean_age_shift_Myr"].abs().median()), 4
            ),
            "note": (
                "computed from the node rules alone, before any injection; this "
                "is the PRIOR mean shift, whereas P4 bounds the POSTERIOR mean "
                "shift, which also carries the likelihood reweighting"
            ),
        },
        "node_plan": node_plan,
        "repair_v5_reference": {
            "regressing_cell": "CygOB2-C MIST R_V=3.5 alpha=2.0",
            "repair_v5_residuals": [
                round(float(x), 3) for x in regressing.pearson_residual
            ],
            "repair_v5_max_abs_residual": round(
                float(np.max(np.abs(regressing.pearson_residual))), 3
            ),
            "repair_v5_baseline_association_mass_Msun": round(
                float(
                    v5_mass[
                        v5_mass.family.eq(BASELINE[0])
                        & v5_mass.R_V.eq(BASELINE[1])
                        & v5_mass.alpha.eq(BASELINE[2])
                    ]["multiplicity_adjusted_mass_median_Msun"].iloc[0]
                ),
                1,
            ),
        },
        "inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in [
                w.PROC / "wp4_age_posteriors_repair_v5.parquet",
                w.PROC / "wp5_mass_function_bins_repair_v5.parquet",
                w.PROC / "wp5_association_mass_repair_v5.parquet",
                w.ROOT / "scripts" / "wp5_common.py",
                w.ROOT / "scripts" / "wp5_joint_age_fit.py",
                w.ROOT / "scripts" / "wp5_injections_repair.py",
            ]
        },
    }
    w.write_json(w.PROVENANCE / "wp5_node_interpolation_prereg.json", record)
    print(
        f"nodes {record['cost']['nodes_repair_v5']} -> "
        f"{record['cost']['nodes_repair_v6']}, all requiring fresh injection"
    )
    print(
        "prior-mean truth age shift: max "
        f"{record['prior_mean_age_shift_forecast']['max_abs_Myr']:.4f} Myr, median "
        f"{record['prior_mean_age_shift_forecast']['median_abs_Myr']:.4f} Myr"
    )
    print("wrote provenance/wp5_node_interpolation_prereg.json")


if __name__ == "__main__":
    main()
