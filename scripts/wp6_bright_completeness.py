#!/usr/bin/env python3
"""Issue #3: the bright-mass completeness WP6's closure test must divide by.

WP6 step 2(a) compares the observed number of massive stars against the number
the WP5 IMF normalization predicts, and reads the shortfall as the population
that has already exploded.  The step as designed assumes every massive star in
the association is in the catalogue.  It is not: the injection experiment
recovers only about four in five, and the loss is *flat in mass* all the way to
18 Msun, so it is the WP2 quality filter rather than any magnitude limit.  Left
uncorrected, that missing fraction is read as a real deficit and inflates N_SN
by the same factor.

This script measures the correction from the injection response that WP5
already produced, and states the estimator WP6 must use instead.  It does not
run WP6 -- WP6 has not started -- but it removes the defect from its design and
quantifies what it would have cost.

The correct comparison is a FORWARD one.  Do not divide an observed count by a
scalar completeness: the response also scatters mass estimates across the
threshold, and above 8 Msun the IMF is steep enough that up-scatter from below
is not negligible.  With k the WP5 normalization, alpha the IMF slope, and
R(observed above threshold | true mass M) the injection response,

    expected observed count above threshold
        = k * integral dM M^-alpha R(above threshold | M)

which is the same construction WP5 already uses inside 2-8 Msun, extended to
one open-ended bin.  The deficit is then observed minus expected, and only that
difference is attributable to stars that have left the main sequence.

Outputs:
  tables/wp6_bright_completeness.csv
  provenance/wp6_bright_completeness_execution.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp6_bright_completeness.py
"""
from __future__ import annotations

import argparse
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
import wp5_joint_age_fit as J

# The census-closure threshold.  8 Msun is the conventional core-collapse floor
# and is also WP5's calibration ceiling, so the closure test starts exactly
# where the normalization stops being constrained by counts.
CLOSURE_THRESHOLD_MSUN = 8.0
PLATEAU_LO_MSUN = 8.0


def branch_response(
    version: str, subgroup: str, family: str, rv: float,
    age_posterior: pd.DataFrame, native: np.ndarray, upstream: str,
) -> tuple[dict[float, pd.DataFrame], dict[float, float]]:
    prior = J.truth_age_nodes(
        age_posterior, subgroup, family, rv, native,
        snap=not J.uses_age_interpolation(version),
    )
    responses = {}
    for age in prior:
        path = J.node_response_path(subgroup, family, rv, age, version)
        if not path.exists():
            tag = f"{age:.3f}".replace(".", "p")
            fallback = w.PROC / f"wp5_age_scan_B_response_age{tag}_repair_v3.parquet"
            if upstream == "repair_v3" and fallback.exists():
                path = fallback
            else:
                raise RuntimeError(f"missing node response {path.name}")
        responses[age] = pd.read_parquet(path)
    return responses, prior


def above_threshold_response(
    response: pd.DataFrame, threshold: float, draw_columns: list[str]
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """R(recovered above threshold | true mass) and the raw recovery fraction."""
    frame = response[
        response["true_primary_mass"].between(
            float(w.MASS_GRID.min()), float(w.MASS_GRID.max())
        )
    ]
    masses = np.sort(frame["true_primary_mass"].unique())
    above = np.zeros(len(masses))
    recovered_any = np.zeros(len(masses))
    for index, value in enumerate(masses):
        rows = frame[frame["true_primary_mass"].eq(value)]
        active = [c for c in draw_columns if rows[c].notna().any()]
        if not active or not len(rows):
            # No injected star at this true mass was recovered at all, so the
            # response to "observed above threshold" is exactly zero.
            above[index] = 0.0
            recovered_any[index] = 0.0
            continue
        draws = rows[active].to_numpy(float)
        finite = np.isfinite(draws)
        above[index] = np.sum(finite & (draws >= threshold)) / (
            len(active) * len(rows)
        )
        recovered_any[index] = float(np.mean(finite.any(axis=1)))
    return masses, above, recovered_any


def imf_trapezoid(masses: np.ndarray, alpha: float) -> np.ndarray:
    weight = np.empty(len(masses))
    weight[0] = 0.5 * (masses[1] - masses[0])
    weight[-1] = 0.5 * (masses[-1] - masses[-2])
    weight[1:-1] = 0.5 * (masses[2:] - masses[:-2])
    return weight * masses ** (-alpha)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", default="repair_v5")
    parser.add_argument("--threshold", type=float, default=CLOSURE_THRESHOLD_MSUN)
    args = parser.parse_args()
    version = args.version
    upstream = "repair_v5" if version in {"repair_v5", "repair_v6"} else "repair_v3"

    normalization = pd.read_parquet(
        w.PROC / f"wp5_imf_normalization_{version}.parquet"
    )
    curves = pd.read_parquet(w.PROC / f"wp5_completeness_curves_{version}.parquet")
    age_posterior = pd.read_parquet(
        w.PROC / f"wp4_age_posteriors_{upstream}.parquet"
    )
    native = {family: J.native_isochrone_ages(family) for family in w.FAMILIES}

    rows = []
    for family in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            for subgroup in w.SUBGROUPS:
                responses, prior = branch_response(
                    version, subgroup, family, rv, age_posterior,
                    native[family], upstream,
                )
                ages = list(prior)
                draw_columns = sorted(
                    c for c in responses[ages[0]].columns
                    if c.startswith("recovered_mass_draw_")
                )
                masses = None
                above = np.zeros(1)
                recovered = np.zeros(1)
                for age in ages:
                    node_masses, node_above, node_recovered = above_threshold_response(
                        responses[age], args.threshold, draw_columns
                    )
                    if masses is None:
                        masses = node_masses
                        above = np.zeros(len(masses))
                        recovered = np.zeros(len(masses))
                    above += prior[age] * node_above
                    recovered += prior[age] * node_recovered

                plateau = masses >= PLATEAU_LO_MSUN
                curve = curves[
                    curves.subgroup.eq(subgroup) & curves.family.eq(family)
                    & curves.R_V.eq(rv)
                ].sort_values("primary_mass")
                isotonic_plateau = float(
                    curve.loc[
                        curve.primary_mass.ge(PLATEAU_LO_MSUN), "recovery_isotonic"
                    ].median()
                )

                for alpha in w.IMF_SLOPES:
                    cell = normalization[
                        normalization.subgroup.eq(subgroup)
                        & normalization.family.eq(family)
                        & normalization.R_V.eq(rv)
                        & normalization.alpha.eq(alpha)
                    ]
                    if len(cell) != 1:
                        continue
                    k = float(cell["k_median"].iloc[0])
                    weight = imf_trapezoid(masses, float(alpha))
                    window = masses >= args.threshold
                    naive = float(k * np.sum(weight[window]))
                    corrected = float(k * np.sum(weight * above))
                    rows.append(
                        {
                            "version": version, "subgroup": subgroup,
                            "family": family, "R_V": float(rv), "alpha": float(alpha),
                            "k_median": k,
                            "threshold_Msun": args.threshold,
                            "expected_true_above_threshold": naive,
                            "expected_observed_above_threshold": corrected,
                            "effective_completeness": corrected / naive if naive else np.nan,
                            "bright_recovery_plateau_response": float(
                                np.mean(recovered[plateau])
                            ),
                            "bright_recovery_plateau_isotonic": isotonic_plateau,
                            "spurious_deficit_fraction": (
                                1.0 - corrected / naive if naive else np.nan
                            ),
                        }
                    )
    table = pd.DataFrame(rows)
    out_csv = w.TABLES / "wp6_bright_completeness.csv"
    table.to_csv(out_csv, index=False)

    baseline = table[
        table.family.eq("PARSEC") & table.R_V.eq(3.1) & table.alpha.eq(2.3)
    ]
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp6_bright_completeness.py",
        "status": "SUCCESS",
        "issue": "#3 — WP6 step 2(a) assumes bright-mass completeness ~1.0",
        "wp5_version": version,
        "defect": (
            "WP6's census-closure step compares the observed massive-star count "
            "against k * integral M^-alpha over M >= threshold, i.e. against the "
            "TRUE number, and reads the shortfall as stars that have already "
            "exploded.  The catalogue does not contain every massive star: the "
            "injection response recovers a fraction that is flat in mass out to "
            "18 Msun, so the loss is the WP2 quality filter rather than a "
            "magnitude limit and does not vanish at the bright end."
        ),
        "consequence_if_uncorrected": (
            "the un-recovered fraction is counted as a real deficit, so the "
            "inferred number of exploded stars — and therefore N_SN — is "
            "inflated by that fraction"
        ),
        "required_estimator": (
            "expected OBSERVED count = k * integral dM M^-alpha * "
            "R(recovered above threshold | M), with R the WP5 injection "
            "response.  This is the same forward construction WP5 already uses "
            "inside 2-8 Msun, extended to one open-ended bin.  Do NOT divide an "
            "observed count by a scalar completeness: the response also "
            "scatters mass estimates across the threshold, and the difference "
            "between the two treatments is reported below as "
            "effective_completeness vs bright_recovery_plateau_response."
        ),
        "threshold_Msun": args.threshold,
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in [
                w.PROC / f"wp5_imf_normalization_{version}.parquet",
                w.PROC / f"wp5_completeness_curves_{version}.parquet",
                w.PROC / f"wp4_age_posteriors_{upstream}.parquet",
            ]
        },
        "grid_summary": {
            "cells": int(len(table)),
            "effective_completeness_min": round(
                float(table.effective_completeness.min()), 4
            ),
            "effective_completeness_median": round(
                float(table.effective_completeness.median()), 4
            ),
            "effective_completeness_max": round(
                float(table.effective_completeness.max()), 4
            ),
            "spurious_deficit_median": round(
                float(table.spurious_deficit_fraction.median()), 4
            ),
            "bright_plateau_recovery_median": round(
                float(table.bright_recovery_plateau_response.median()), 4
            ),
        },
        "baseline_PARSEC_rv3.1_alpha2.3": [
            {
                "subgroup": row.subgroup,
                "k_median": round(row.k_median, 2),
                "expected_true_above_8Msun": round(
                    row.expected_true_above_threshold, 2
                ),
                "expected_observed_above_8Msun": round(
                    row.expected_observed_above_threshold, 2
                ),
                "effective_completeness": round(row.effective_completeness, 4),
                "bright_recovery_plateau": round(
                    row.bright_recovery_plateau_response, 4
                ),
            }
            for row in baseline.itertuples()
        ],
        "outputs": {str(out_csv.relative_to(w.ROOT)): w.sha256(out_csv)},
    }
    w.write_json(w.PROVENANCE / "wp6_bright_completeness_execution.json", record)

    print(f"{version}: effective completeness above {args.threshold} Msun")
    print("  subgroup      expected true   expected observed   effective   plateau")
    for entry in record["baseline_PARSEC_rv3.1_alpha2.3"]:
        print(
            f"  {entry['subgroup']:12s} {entry['expected_true_above_8Msun']:12.2f}"
            f"   {entry['expected_observed_above_8Msun']:16.2f}"
            f"   {entry['effective_completeness']:9.3f}"
            f"   {entry['bright_recovery_plateau']:7.3f}"
        )
    summary = record["grid_summary"]
    print(
        f"\nover all {summary['cells']} cells: effective completeness "
        f"{summary['effective_completeness_min']:.3f}-"
        f"{summary['effective_completeness_max']:.3f}, median "
        f"{summary['effective_completeness_median']:.3f}"
    )
    print(
        f"a WP6 closure test assuming completeness 1.0 would report a spurious "
        f"deficit of {summary['spurious_deficit_median']:.1%} and inflate N_SN "
        f"by the same factor"
    )
    print("wrote provenance/wp6_bright_completeness_execution.json")


if __name__ == "__main__":
    main()
