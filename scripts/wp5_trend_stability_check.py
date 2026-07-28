#!/usr/bin/env python3
"""Is the repair_v4 A/C gate flip a model regression or Monte-Carlo noise?

The repair_v4 G3 check flagged two CygOB2-C cells (MIST, R_V=3.1, alpha=2.0
and 2.3) as pass -> fail.  Both are **single-truth-age-node** cells, where
``fit_joint`` provably reduces to the unmodified ``wp5_fit_imf.fit_one``
bit-for-bit, so the estimator there is identical to repair_v3's and the flip
cannot have been caused by the age marginalization.  Both also *improved* in
chi-square and in max absolute residual; only the rank-based trend statistic
moved.

This script measures the alternative explanation directly.  For every
single-node cell it fits the identical estimator to two independent Monte-Carlo
realizations of the same model -- the stored repair_v3 injection response and
the freshly generated repair_v4 node response, which for these cells were
generated at the *same* truth age -- and reports how far each gate statistic
moves.  It also scans the fit RNG seed to separate injection noise from
posterior-draw noise.

The 6-bin two-sided Spearman p-value is quantized (rho = 1.000 -> p = 0.0028,
0.943 -> 0.0048, 0.886 -> 0.0333, 0.829 -> 0.0583, 0.771 -> 0.1028, ...), so a
small perturbation of a well-fitting residual vector can cross p = 0.05 while
chi-square and max residual barely move.  That is the hypothesis under test.

Output: provenance/wp5_trend_stability_check_execution.json

Run:
  WP_REPAIR_VERSION=repair_v3 WP3_ANCHOR_PRIOR_MODE=variogram \
  PYTHONPATH=scripts python3 scripts/wp5_trend_stability_check.py
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy import stats

import wp5_common as w
import wp5_joint_age_fit as J

UPSTREAM = "repair_v3"
VERSION = "repair_v4"
SEED_SCAN = [w.SEED, w.SEED + 1, w.SEED + 2, w.SEED + 3, w.SEED + 4]


def main() -> None:
    masses = pd.read_parquet(w.PROC / f"wp4_mass_posteriors_{UPSTREAM}.parquet")
    store = np.load(w.PROC / f"wp4_mass_posterior_samples_{UPSTREAM}.npz")
    age_posterior = pd.read_parquet(w.PROC / f"wp4_age_posteriors_{UPSTREAM}.parquet")
    old_curves = pd.read_parquet(w.PROC / f"wp5_completeness_curves_{UPSTREAM}.parquet")
    old_responses = pd.read_parquet(w.PROC / f"wp5_injection_response_{UPSTREAM}.parquet")
    old_norm = pd.read_parquet(w.PROC / f"wp5_imf_normalization_{UPSTREAM}.parquet")
    new_norm = pd.read_parquet(w.PROC / f"wp5_imf_normalization_{VERSION}.parquet")
    draw_columns = sorted(
        c for c in old_responses.columns if c.startswith("recovered_mass_draw_")
    )
    native = {family: J.native_isochrone_ages(family) for family in w.FAMILIES}

    # Quantization of the 6-bin two-sided Spearman p-value.
    centers = np.log10(np.sqrt(
        np.geomspace(2.0, 8.0, 7)[:-1] * np.geomspace(2.0, 8.0, 7)[1:]
    ))
    lattice = {}
    for permutation in range(200):
        rng = np.random.default_rng(permutation)
        values = rng.permutation(6).astype(float)
        result = stats.spearmanr(centers, values)
        lattice[round(float(result.statistic), 3)] = round(float(result.pvalue), 4)
    quantization = dict(sorted(lattice.items()))

    flips = []
    for row in new_norm.itertuples():
        before = old_norm[
            old_norm.subgroup.eq(row.subgroup) & old_norm.family.eq(row.family)
            & old_norm.R_V.eq(row.R_V) & old_norm.alpha.eq(row.alpha)
        ].iloc[0]
        if bool(before.residual_gate_pass) == bool(row.residual_gate_pass):
            continue
        nodes = J.truth_age_nodes(
            age_posterior, row.subgroup, row.family, row.R_V, native[row.family]
        )
        flips.append(
            {
                "subgroup": row.subgroup,
                "family": row.family,
                "R_V": float(row.R_V),
                "alpha": float(row.alpha),
                "direction": (
                    "pass_to_fail" if before.residual_gate_pass else "fail_to_pass"
                ),
                "n_truth_age_nodes": len(nodes),
                "model_identical_to_repair_v3": len(nodes) == 1,
                "before": {
                    "chi2_p": float(before.poisson_chi_square_p),
                    "trend_p": float(before.residual_trend_p),
                    "rho": float(before.residual_spearman_rho),
                    "max_abs_residual": float(before.max_abs_pearson_residual),
                },
                "after": {
                    "chi2_p": float(row.poisson_chi_square_p),
                    "trend_p": float(row.residual_trend_p),
                    "rho": float(row.residual_spearman_rho),
                    "max_abs_residual": float(row.max_abs_pearson_residual),
                },
            }
        )

    # Paired refit of every single-node flipped cell on both MC realizations.
    paired = []
    for flip in flips:
        if not flip["model_identical_to_repair_v3"]:
            continue
        subgroup, family, rv, alpha = (
            flip["subgroup"], flip["family"], flip["R_V"], flip["alpha"]
        )
        nodes = J.truth_age_nodes(age_posterior, subgroup, family, rv, native[family])
        age = next(iter(nodes))
        branch_samples = store["samples"][
            :,
            w.FAMILIES.index(family) * len(w.R_V_BRANCHES)
            + list(w.R_V_BRANCHES).index(rv),
            :,
        ]
        old_curve = old_curves[
            old_curves.family.eq(family) & old_curves.R_V.eq(rv)
            & old_curves.subgroup.eq(subgroup)
        ]
        old_response = old_responses[
            old_responses.family.eq(family) & old_responses.R_V.eq(rv)
            & old_responses.subgroup.eq(subgroup)
        ]
        new_curve = pd.read_parquet(J.node_curve_path(subgroup, family, rv, age, VERSION))
        new_response = pd.read_parquet(
            J.node_response_path(subgroup, family, rv, age, VERSION)
        )
        realizations = {}
        for label, curve, response in [
            ("repair_v3_injection_realization", old_curve, old_response),
            ("repair_v4_injection_realization", new_curve, new_response),
        ]:
            per_seed = []
            for seed in SEED_SCAN:
                summary, _, _, _ = J.fit_joint(
                    masses, {age: curve}, {age: response}, {age: 1.0},
                    subgroup, family, rv, alpha,
                    np.random.default_rng(seed), branch_samples, draw_columns,
                )
                per_seed.append(
                    {
                        "seed": int(seed),
                        "chi2_p": float(summary["poisson_chi_square_p"]),
                        "trend_p": float(summary["residual_trend_p"]),
                        "rho": float(summary["residual_spearman_rho"]),
                        "max_abs_residual": float(summary["max_abs_pearson_residual"]),
                        "gate_pass": bool(summary["residual_gate_pass"]),
                    }
                )
            realizations[label] = {
                "truth_age_Myr": float(age),
                "per_seed": per_seed,
                "trend_p_range": [
                    min(entry["trend_p"] for entry in per_seed),
                    max(entry["trend_p"] for entry in per_seed),
                ],
                "gate_pass_fraction": float(
                    np.mean([entry["gate_pass"] for entry in per_seed])
                ),
            }
        paired.append(
            {
                "subgroup": subgroup, "family": family, "R_V": rv, "alpha": alpha,
                "same_truth_age_both_versions": True,
                "realizations": realizations,
                "verdict": (
                    "gate outcome differs between two independent Monte-Carlo "
                    "realizations of an IDENTICAL model"
                    if realizations["repair_v3_injection_realization"][
                        "gate_pass_fraction"
                    ]
                    != realizations["repair_v4_injection_realization"][
                        "gate_pass_fraction"
                    ]
                    else "gate outcome stable across realizations"
                ),
            }
        )

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_trend_stability_check.py",
        "status": "SUCCESS",
        "issue": "#1c step 3c — are the flagged A/C gate flips model or Monte Carlo?",
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "spearman_p_quantization_n6": quantization,
        "gate_flips_v3_to_v4": flips,
        "paired_realization_test": paired,
        "inputs": {
            str(p.relative_to(w.ROOT)): w.sha256(p)
            for p in [
                w.PROC / f"wp5_imf_normalization_{UPSTREAM}.parquet",
                w.PROC / f"wp5_imf_normalization_{VERSION}.parquet",
                w.PROC / f"wp5_injection_response_{UPSTREAM}.parquet",
            ]
        },
    }
    w.write_json(w.PROVENANCE / "wp5_trend_stability_check_execution.json", record)
    print("Spearman p quantization (n=6):", quantization)
    print(f"\n{len(flips)} gate flips v3 -> v4:")
    for flip in flips:
        print(
            f"  {flip['subgroup']} {flip['family']} rv{flip['R_V']} a{flip['alpha']}: "
            f"{flip['direction']}, nodes={flip['n_truth_age_nodes']}, "
            f"model_identical={flip['model_identical_to_repair_v3']}"
        )
        print(
            f"     chi2p {flip['before']['chi2_p']:.4f}->{flip['after']['chi2_p']:.4f}  "
            f"trendp {flip['before']['trend_p']:.4f}->{flip['after']['trend_p']:.4f}  "
            f"max|r| {flip['before']['max_abs_residual']:.3f}->{flip['after']['max_abs_residual']:.3f}"
        )
    print("\npaired realization test:")
    for entry in paired:
        print(f"  {entry['subgroup']} {entry['family']} rv{entry['R_V']} a{entry['alpha']}")
        for label, data in entry["realizations"].items():
            print(
                f"     {label}: trend_p range {data['trend_p_range'][0]:.4f}"
                f"-{data['trend_p_range'][1]:.4f}, gate pass fraction "
                f"{data['gate_pass_fraction']:.2f}"
            )
        print(f"     -> {entry['verdict']}")
    print("\nwrote provenance/wp5_trend_stability_check_execution.json")


if __name__ == "__main__":
    main()
