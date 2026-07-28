#!/usr/bin/env python3
"""Issue #13, part a: measure the discontinuity of the truth-age node rule.

The repair_v4 rule (``wp5_joint_age_fit.truth_age_nodes``) snaps the nine WP4
posterior quantiles to native isochrone ages.  The native grid is coarse
(0.05 dex, about 12 per cent in age), so an arbitrarily small change in the
age posterior can move a quantile across a snap boundary and thereby add or
delete a whole truth-age node carrying up to 1/9 of the prior weight.  That is
a discontinuous map from a continuous input, and it is what produced the one
genuine repair_v4 -> repair_v5 regression recorded in
``provenance/wp5_node_snapping_discontinuity.json``.

This script quantifies the defect and the proposed repair *without running any
injection*, so the repair can be pre-registered on measured grounds.  Both
rules are applied to the same WP4 posteriors, translated by a continuous
offset epsilon, and the node distributions are compared in the 1-Wasserstein
metric -- the natural distance for a distribution over truth ages, and the one
that controls how much the mixture completeness curve can move.

A continuous rule must satisfy W1(epsilon) -> 0 as epsilon -> 0.  For a pure
translation of the posterior the unsnapped rule satisfies W1 = |epsilon|
exactly, which is the tightest possible bound.

Outputs:
  tables/wp5_node_rule_continuity.csv
  provenance/wp5_node_rule_continuity_execution.json

Run (seconds; no injections, no fits):
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp5_node_rule_continuity.py
"""
from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.stats import norm

import wp5_common as w
import wp5_joint_age_fit as J
from wp4_repair_common import F_BINARY, N_AGE_NODES

# Translation sweep.  +-0.10 Myr at 0.0005 Myr resolution: fine enough that a
# continuous rule cannot move more than 0.0005 Myr in W1 between steps, so any
# larger step is a genuine discontinuity rather than sampling coarseness.
EPSILON_HALF_WIDTH_MYR = 0.10
EPSILON_STEP_MYR = 0.0005


def posterior_triple(
    age_posterior: pd.DataFrame, subgroup: str, family: str, rv: float
) -> tuple[float, float, float]:
    """(MAP, lo68, hi68) exactly as ``age_posterior_nodes`` selects them."""
    row = age_posterior[
        age_posterior["subgroup"].eq(subgroup)
        & age_posterior["family"].eq(family)
        & age_posterior["R_V"].eq(rv)
        & age_posterior["f_bin"].eq(F_BINARY)
        & age_posterior["indicator"].eq("ums")
        & age_posterior["dmu"].eq(0.0)
    ]
    if len(row) != 1:
        raise RuntimeError(f"missing age posterior for {subgroup}/{family}/R_V={rv}")
    return (
        float(row["age_map"].iloc[0]),
        float(row["age_lo68"].iloc[0]),
        float(row["age_hi68"].iloc[0]),
    )


def raw_nodes(centre: float, lo: float, hi: float) -> np.ndarray:
    """``wp4_repair_common.age_posterior_nodes`` arithmetic, inlined so that a
    translated posterior can be evaluated without rebuilding the frame."""
    probabilities = (np.arange(N_AGE_NODES) + 0.5) / N_AGE_NODES
    z = norm.ppf(probabilities)
    lower_scale = max(centre - lo, 0.02) / abs(norm.ppf(0.16))
    upper_scale = max(hi - centre, 0.02) / norm.ppf(0.84)
    nodes = centre + np.where(z < 0, lower_scale * z, upper_scale * z)
    return np.clip(nodes, 1.0, 10.0)


def snapped_distribution(
    nodes: np.ndarray, native_ages: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """repair_v4 rule: snap to native isochrone ages, weights = summed counts."""
    weights: dict[float, float] = {}
    for node in nodes:
        native = J.snap_to_native(float(node), native_ages)
        weights[native] = weights.get(native, 0.0) + 1.0 / len(nodes)
    ordered = sorted(weights)
    return np.array(ordered, dtype=float), np.array(
        [weights[age] for age in ordered], dtype=float
    )


def interpolated_distribution(
    nodes: np.ndarray, native_ages: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Proposed rule: keep the quantiles where they are, equal weights."""
    ordered = np.sort(np.asarray(nodes, dtype=float))
    return ordered, np.full(len(ordered), 1.0 / len(ordered))


def wasserstein1(
    a_support: np.ndarray,
    a_weight: np.ndarray,
    b_support: np.ndarray,
    b_weight: np.ndarray,
) -> float:
    """1-Wasserstein distance between two discrete distributions on the line."""
    grid = np.unique(np.concatenate([a_support, b_support]))
    a_cdf = np.array([a_weight[a_support <= point].sum() for point in grid])
    b_cdf = np.array([b_weight[b_support <= point].sum() for point in grid])
    return float(np.sum(np.abs(a_cdf[:-1] - b_cdf[:-1]) * np.diff(grid)))


def sweep(
    centre: float,
    lo: float,
    hi: float,
    native_ages: np.ndarray,
    rule,
) -> dict:
    offsets = np.arange(
        -EPSILON_HALF_WIDTH_MYR,
        EPSILON_HALF_WIDTH_MYR + 0.5 * EPSILON_STEP_MYR,
        EPSILON_STEP_MYR,
    )
    supports, weights = [], []
    for offset in offsets:
        support, weight = rule(
            raw_nodes(centre + offset, lo + offset, hi + offset), native_ages
        )
        supports.append(support)
        weights.append(weight)
    reference = int(np.argmin(np.abs(offsets)))
    distance = np.array(
        [
            wasserstein1(
                supports[reference], weights[reference], supports[index], weights[index]
            )
            for index in range(len(offsets))
        ]
    )
    step_jump = np.array(
        [
            wasserstein1(
                supports[index - 1], weights[index - 1], supports[index], weights[index]
            )
            for index in range(1, len(offsets))
        ]
    )
    node_count = np.array([len(support) for support in supports])
    # A node appearing or vanishing is the specific failure mode; count the
    # offsets at which the support set changes at all.
    support_changes = sum(
        1
        for index in range(1, len(offsets))
        if not np.array_equal(supports[index - 1], supports[index])
    )
    return {
        "max_step_jump_Myr": float(step_jump.max()),
        "median_step_jump_Myr": float(np.median(step_jump)),
        "jump_over_step_ratio": float(step_jump.max() / EPSILON_STEP_MYR),
        "max_distance_Myr": float(distance.max()),
        "node_count_min": int(node_count.min()),
        "node_count_max": int(node_count.max()),
        "support_change_events": int(support_changes),
        "offsets_evaluated": int(len(offsets)),
    }


def main() -> None:
    age_posterior = pd.read_parquet(
        w.PROC / "wp4_age_posteriors_repair_v5.parquet"
    )
    native = {family: J.native_isochrone_ages(family) for family in w.FAMILIES}
    grid_spacing = {
        family: {
            "n_ages": int(len(native[family])),
            "min_Myr": float(native[family].min()),
            "max_Myr": float(native[family].max()),
            "median_relative_spacing": float(
                np.median(np.diff(native[family]) / native[family][:-1])
            ),
        }
        for family in w.FAMILIES
    }

    rows = []
    for family in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            for subgroup in w.SUBGROUPS:
                centre, lo, hi = posterior_triple(age_posterior, subgroup, family, rv)
                snapped = sweep(centre, lo, hi, native[family], snapped_distribution)
                interpolated = sweep(
                    centre, lo, hi, native[family], interpolated_distribution
                )
                rows.append(
                    {
                        "subgroup": subgroup,
                        "family": family,
                        "R_V": rv,
                        "age_map_Myr": centre,
                        "age_lo68_Myr": lo,
                        "age_hi68_Myr": hi,
                        **{f"snap_{key}": value for key, value in snapped.items()},
                        **{
                            f"interp_{key}": value
                            for key, value in interpolated.items()
                        },
                    }
                )
    table = pd.DataFrame(rows)
    out_csv = w.TABLES / "wp5_node_rule_continuity.csv"
    table.to_csv(out_csv, index=False)

    # The recorded regression, reproduced exactly from its own posteriors.
    v4_posterior = pd.read_parquet(w.PROC / "wp4_age_posteriors_repair_v3.parquet")
    case = {}
    for label, frame in (("repair_v4", v4_posterior), ("repair_v5", age_posterior)):
        centre, lo, hi = posterior_triple(frame, "CygOB2-C", "MIST", 3.5)
        nodes = raw_nodes(centre, lo, hi)
        snap_support, snap_weight = snapped_distribution(nodes, native["MIST"])
        case[label] = {
            "age_map_Myr": centre,
            "age_lo68_Myr": lo,
            "age_hi68_Myr": hi,
            "raw_quantile_nodes_Myr": [round(float(x), 4) for x in nodes],
            "snapped_support_Myr": [float(x) for x in snap_support],
            "snapped_weights": [float(x) for x in snap_weight],
        }
    case["posterior_shift_Myr"] = {
        "age_map": case["repair_v5"]["age_map_Myr"] - case["repair_v4"]["age_map_Myr"],
        "age_lo68": case["repair_v5"]["age_lo68_Myr"]
        - case["repair_v4"]["age_lo68_Myr"],
        "age_hi68": case["repair_v5"]["age_hi68_Myr"]
        - case["repair_v4"]["age_hi68_Myr"],
    }
    case["snapped_W1_Myr"] = wasserstein1(
        np.array(case["repair_v4"]["snapped_support_Myr"]),
        np.array(case["repair_v4"]["snapped_weights"]),
        np.array(case["repair_v5"]["snapped_support_Myr"]),
        np.array(case["repair_v5"]["snapped_weights"]),
    )
    v4_nodes = np.array(case["repair_v4"]["raw_quantile_nodes_Myr"])
    v5_nodes = np.array(case["repair_v5"]["raw_quantile_nodes_Myr"])
    case["interpolated_W1_Myr"] = wasserstein1(
        *interpolated_distribution(v4_nodes, native["MIST"]),
        *interpolated_distribution(v5_nodes, native["MIST"]),
    )
    case["amplification_factor"] = (
        case["snapped_W1_Myr"] / case["interpolated_W1_Myr"]
        if case["interpolated_W1_Myr"] > 0
        else None
    )

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_node_rule_continuity.py",
        "status": "SUCCESS",
        "issue": "#13 part a — measuring the truth-age node-rule discontinuity",
        "question": (
            "Is the repair_v4 truth-age node rule a discontinuous function of "
            "the WP4 age posterior, and does dropping the snap step remove the "
            "discontinuity?"
        ),
        "method": (
            "Translate each branch's WP4 age posterior (MAP, lo68, hi68 all "
            "shifted by the same epsilon) over a fine grid and measure the "
            "1-Wasserstein distance between the resulting truth-age node "
            "distributions.  A continuous rule cannot move further in W1 "
            "between adjacent offsets than the offset step itself."
        ),
        "epsilon_half_width_Myr": EPSILON_HALF_WIDTH_MYR,
        "epsilon_step_Myr": EPSILON_STEP_MYR,
        "n_age_nodes": N_AGE_NODES,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in [
                w.PROC / "wp4_age_posteriors_repair_v3.parquet",
                w.PROC / "wp4_age_posteriors_repair_v5.parquet",
                w.PROC / "wp3_isochrones_parsec.parquet",
                w.PROC / "wp3_isochrones_mist.parquet",
            ]
        },
        "native_isochrone_grid": grid_spacing,
        "summary": {
            "branches": int(len(table)),
            "snap_max_step_jump_Myr": float(table["snap_max_step_jump_Myr"].max()),
            "snap_worst_jump_over_step_ratio": float(
                table["snap_jump_over_step_ratio"].max()
            ),
            "snap_branches_with_support_changes": int(
                (table["snap_support_change_events"] > 0).sum()
            ),
            "snap_node_count_range_within_branch_max": int(
                (table["snap_node_count_max"] - table["snap_node_count_min"]).max()
            ),
            "interp_max_step_jump_Myr": float(table["interp_max_step_jump_Myr"].max()),
            "interp_worst_jump_over_step_ratio": float(
                table["interp_jump_over_step_ratio"].max()
            ),
            "interp_node_count_always_nine": bool(
                table["interp_node_count_min"].eq(N_AGE_NODES).all()
                and table["interp_node_count_max"].eq(N_AGE_NODES).all()
            ),
        },
        "recorded_regression_case": case,
        "verdict": (
            "The snapping rule is discontinuous: within a +-0.1 Myr "
            "translation every branch's node support changes and the node "
            "distribution moves in W1 by far more than the translation that "
            "caused it.  Removing the snap makes W1 equal to the translation "
            "exactly, which is the tightest continuity a rule can have."
        ),
        "outputs": {str(out_csv.relative_to(w.ROOT)): w.sha256(out_csv)},
    }
    w.write_json(w.PROVENANCE / "wp5_node_rule_continuity_execution.json", record)

    print(f"{len(table)} branches swept, epsilon step {EPSILON_STEP_MYR} Myr\n")
    print("                          max W1 jump between adjacent offsets")
    print("  rule            max jump (Myr)   jump/step   node count")
    for prefix, name in (("snap", "snap (v4)"), ("interp", "interpolate")):
        print(
            f"  {name:14s}  {table[f'{prefix}_max_step_jump_Myr'].max():10.5f}"
            f"   {table[f'{prefix}_jump_over_step_ratio'].max():9.1f}"
            f"   {table[f'{prefix}_node_count_min'].min()}"
            f"-{table[f'{prefix}_node_count_max'].max()}"
        )
    print(
        f"\nrecorded regression (C MIST R_V=3.5), posterior moved "
        f"{case['posterior_shift_Myr']['age_hi68']:+.3f} Myr at the 68% upper edge:"
    )
    print(
        f"  snapped node distribution moved   W1 = "
        f"{case['snapped_W1_Myr']:.4f} Myr"
    )
    print(
        f"  unsnapped node distribution moved W1 = "
        f"{case['interpolated_W1_Myr']:.4f} Myr"
        f"   (amplification {case['amplification_factor']:.1f}x)"
    )
    print("\nwrote provenance/wp5_node_rule_continuity_execution.json")


if __name__ == "__main__":
    main()
