#!/usr/bin/env python3
"""Validate the replacement trend statistic against criteria R1-R3.

Pre-declared in CUTS_AND_THRESHOLDS.md section 14.4.  The replacement is
adopted only if all four criteria hold; R4 (equal-footing re-scoring of every
version) is run by scripts/wp5_rescore_all_versions.py.

R1 Calibration -- false-positive rate on simulated null data drawn from the
   real fitted expectations of every repair_v4 cell must lie in [0.04, 0.06].
   Reported next to the incumbent's rate on the *same* data, so the comparison
   is on realistic low-count Poisson bins rather than an idealized permutation
   null.

R2 Power -- against residual drifts of injected known slope, the replacement
   must detect at least as well as the incumbent **at matched false-positive
   rate**.  Matching the false-positive rate is essential: a test that merely
   passes more cells is a weakened gate, and section 14.4 rejects it.

R3 Stability -- on the four identical-model cell pairs of issue #11, both
   Monte-Carlo realizations must return the same verdict.

Output: provenance/wp5_trend_replacement_validation_execution.json

Run:
  WP_REPAIR_VERSION=repair_v3 WP3_ANCHOR_PRIOR_MODE=variogram \
  PYTHONPATH=scripts python3 scripts/wp5_trend_replacement_validation.py
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
import wp5_residual_trend as T

N_CALIBRATION = 4_000
N_POWER = 2_000
SLOPE_GRID = [0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8]
VERSION = "repair_v4"


def cell_inputs(bins: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    ordered = bins.sort_values("bin_index")
    log_mass = np.log10(ordered["mass_geometric_center"].to_numpy(float))
    expected = ordered["expected_count_at_k_median"].to_numpy(float)
    rate = ordered["completeness_weighted_imf_integral_per_k"].to_numpy(float)
    return log_mass, expected, rate


def simulate_datasets(
    expected: np.ndarray,
    rate: np.ndarray,
    log_mass: np.ndarray,
    n: int,
    rng: np.random.Generator,
    slope: float = 0.0,
) -> np.ndarray:
    """Poisson datasets; ``slope`` tilts the truth to test power."""
    centred = log_mass - log_mass.mean()
    tilt = np.exp(slope * centred)
    lam_true = np.maximum(expected * tilt, 1e-12)
    simulated = rng.poisson(lam_true, size=(n, len(expected)))
    total_rate = float(np.sum(rate))
    totals = simulated.sum(axis=1).astype(float)
    k = np.array([T.jeffreys_k(value, total_rate) for value in totals])
    lam_fit = np.maximum(k[:, None] * rate[None, :], 1e-12)
    return (simulated - lam_fit) / np.sqrt(lam_fit)


def main() -> None:
    bins_all = pd.read_parquet(w.PROC / f"wp5_mass_function_bins_{VERSION}.parquet")
    cells = list(bins_all.groupby(["subgroup", "family", "R_V", "alpha"]))
    rng = np.random.default_rng(T.BOOTSTRAP_SEED)

    # ---------------- R1 calibration ----------------
    incumbent_reject = []
    replacement_reject = []
    per_cell = []
    for (subgroup, family, rv, alpha), cell in cells:
        log_mass, expected, rate = cell_inputs(cell)
        null = T.bootstrap_null(expected, rate, log_mass)
        residuals = simulate_datasets(
            expected, rate, log_mass, N_CALIBRATION, rng, slope=0.0
        )
        inc = np.array(
            [T.incumbent_trend_p(log_mass, row) for row in residuals]
        )
        obs = residuals @ (log_mass - log_mass.mean()) / np.sqrt(
            np.sum((log_mass - log_mass.mean()) ** 2)
        )
        rep = np.array([np.mean(np.abs(null) >= abs(value)) for value in obs])
        incumbent_reject.append(float(np.mean(inc < T.TREND_THRESHOLD)))
        replacement_reject.append(float(np.mean(rep < T.TREND_THRESHOLD)))
        per_cell.append(
            {
                "subgroup": subgroup, "family": family, "R_V": float(rv),
                "alpha": float(alpha),
                "incumbent_false_positive_rate": incumbent_reject[-1],
                "replacement_false_positive_rate": replacement_reject[-1],
                "min_expected_count": float(expected.min()),
            }
        )
    r1 = {
        "n_cells": len(per_cell),
        "n_simulated_null_datasets_per_cell": N_CALIBRATION,
        "incumbent_false_positive_rate_mean": float(np.mean(incumbent_reject)),
        "incumbent_false_positive_rate_range": [
            float(np.min(incumbent_reject)), float(np.max(incumbent_reject))
        ],
        "replacement_false_positive_rate_mean": float(np.mean(replacement_reject)),
        "replacement_false_positive_rate_range": [
            float(np.min(replacement_reject)), float(np.max(replacement_reject))
        ],
        "required_interval": [0.04, 0.06],
        "pass": bool(0.04 <= float(np.mean(replacement_reject)) <= 0.06),
    }

    # ---------------- R2 power, at MATCHED false-positive rate ----------------
    # The incumbent's achievable rate is lattice-constrained, so both tests are
    # compared at their own empirical 5% critical values rather than at a
    # nominal threshold neither hits exactly.
    power_cells = [
        entry for entry in cells
        if entry[0][1] == "PARSEC" and entry[0][2] == 3.1 and entry[0][3] == 2.3
    ]
    power_rows = []
    for (subgroup, family, rv, alpha), cell in power_cells:
        log_mass, expected, rate = cell_inputs(cell)
        null_residuals = simulate_datasets(
            expected, rate, log_mass, 20_000, rng, slope=0.0
        )
        centred = log_mass - log_mass.mean()
        null_T = np.abs(
            null_residuals @ centred / np.sqrt(np.sum(centred**2))
        )
        null_rho = np.abs(
            np.array([stats.spearmanr(log_mass, row).statistic for row in null_residuals])
        )
        crit_T = float(np.quantile(null_T, 0.95))
        crit_rho = float(np.quantile(null_rho, 0.95))
        for slope in SLOPE_GRID:
            residuals = simulate_datasets(
                expected, rate, log_mass, N_POWER, rng, slope=slope
            )
            stat_T = np.abs(residuals @ centred / np.sqrt(np.sum(centred**2)))
            stat_rho = np.abs(
                np.array([stats.spearmanr(log_mass, row).statistic for row in residuals])
            )
            power_rows.append(
                {
                    "subgroup": subgroup,
                    "injected_slope": slope,
                    "incumbent_power": float(np.mean(stat_rho >= crit_rho)),
                    "replacement_power": float(np.mean(stat_T >= crit_T)),
                }
            )
    non_null = [row for row in power_rows if row["injected_slope"] > 0]
    deficits = [
        row["incumbent_power"] - row["replacement_power"] for row in non_null
    ]
    r2 = {
        "matched_false_positive_rate": 0.05,
        "critical_values": "empirical 95th percentile of each statistic under its own null",
        "per_slope": power_rows,
        "max_power_deficit_vs_incumbent": float(max(deficits)),
        "mean_power_advantage_of_replacement": float(-np.mean(deficits)),
        "pass": bool(max(deficits) <= 0.01),
    }

    # ---------------- R3 stability on the issue #11 cell pairs ----------------
    stability_record = json.loads(
        (w.PROVENANCE / "wp5_trend_stability_check_execution.json").read_text(
            encoding="utf-8"
        )
    )
    identical_cells = [
        flip for flip in stability_record["gate_flips_v3_to_v4"]
        if flip["model_identical_to_repair_v3"]
    ]
    old_bins = pd.read_parquet(w.PROC / "wp5_mass_function_bins_repair_v3.parquet")
    new_bins = pd.read_parquet(w.PROC / f"wp5_mass_function_bins_{VERSION}.parquet")
    r3_rows = []
    for flip in identical_cells:
        verdicts = {}
        for label, frame in [("repair_v3", old_bins), ("repair_v4", new_bins)]:
            cell = frame[
                frame.subgroup.eq(flip["subgroup"]) & frame.family.eq(flip["family"])
                & frame.R_V.eq(flip["R_V"]) & frame.alpha.eq(flip["alpha"])
            ]
            log_mass, expected, rate = cell_inputs(cell)
            residual = cell.sort_values("bin_index")["pearson_residual"].to_numpy(float)
            p_new, t_obs = T.replacement_trend_p(residual, expected, rate, log_mass)
            verdicts[label] = {
                "incumbent_p": T.incumbent_trend_p(log_mass, residual),
                "replacement_p": p_new,
                "T": t_obs,
                "incumbent_pass": T.incumbent_trend_p(log_mass, residual) >= T.TREND_THRESHOLD,
                "replacement_pass": p_new >= T.TREND_THRESHOLD,
            }
        r3_rows.append(
            {
                "subgroup": flip["subgroup"], "family": flip["family"],
                "R_V": flip["R_V"], "alpha": flip["alpha"],
                "realizations": verdicts,
                "incumbent_consistent": verdicts["repair_v3"]["incumbent_pass"]
                == verdicts["repair_v4"]["incumbent_pass"],
                "replacement_consistent": verdicts["repair_v3"]["replacement_pass"]
                == verdicts["repair_v4"]["replacement_pass"],
            }
        )
    r3 = {
        "n_identical_model_cell_pairs": len(r3_rows),
        "incumbent_consistent_pairs": sum(
            row["incumbent_consistent"] for row in r3_rows
        ),
        "replacement_consistent_pairs": sum(
            row["replacement_consistent"] for row in r3_rows
        ),
        "pairs": r3_rows,
        "pass": bool(all(row["replacement_consistent"] for row in r3_rows)),
    }

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_trend_replacement_validation.py",
        "status": "SUCCESS",
        "issue": "#11 — validation of the pre-declared replacement trend statistic",
        "predeclaration": "CUTS_AND_THRESHOLDS.md section 14",
        "seed": T.BOOTSTRAP_SEED,
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "inputs": {
            str(p.relative_to(w.ROOT)): w.sha256(p) for p in [
                w.PROC / f"wp5_mass_function_bins_{VERSION}.parquet",
                w.PROC / "wp5_mass_function_bins_repair_v3.parquet",
                w.PROVENANCE / "wp5_trend_stability_check_execution.json",
            ]
        },
        "R1_calibration": r1,
        "R2_power": r2,
        "R3_stability": r3,
        "all_criteria_pass": bool(r1["pass"] and r2["pass"] and r3["pass"]),
        "per_cell_calibration": per_cell,
    }
    w.write_json(
        w.PROVENANCE / "wp5_trend_replacement_validation_execution.json", record
    )
    print("R1 calibration:")
    print(f"   incumbent  FPR mean {r1['incumbent_false_positive_rate_mean']:.4f} "
          f"range {r1['incumbent_false_positive_rate_range']}")
    print(f"   replacement FPR mean {r1['replacement_false_positive_rate_mean']:.4f} "
          f"range {r1['replacement_false_positive_rate_range']}  -> "
          f"{'PASS' if r1['pass'] else 'FAIL'}")
    print("\nR2 power at matched 5% false-positive rate:")
    for row in power_rows:
        print(f"   {row['subgroup']} slope {row['injected_slope']:.1f}: "
              f"incumbent {row['incumbent_power']:.3f}  "
              f"replacement {row['replacement_power']:.3f}")
    print(f"   max deficit {r2['max_power_deficit_vs_incumbent']:+.4f} -> "
          f"{'PASS' if r2['pass'] else 'FAIL'}")
    print("\nR3 stability on identical-model cell pairs:")
    for row in r3_rows:
        v3, v4 = row["realizations"]["repair_v3"], row["realizations"]["repair_v4"]
        print(f"   {row['subgroup']} {row['family']} rv{row['R_V']} a{row['alpha']}: "
              f"incumbent p {v3['incumbent_p']:.4f}/{v4['incumbent_p']:.4f} "
              f"({'consistent' if row['incumbent_consistent'] else 'FLIPS'})  "
              f"replacement p {v3['replacement_p']:.4f}/{v4['replacement_p']:.4f} "
              f"({'consistent' if row['replacement_consistent'] else 'FLIPS'})")
    print(f"   -> {'PASS' if r3['pass'] else 'FAIL'}")
    print(f"\nR1-R3 all pass: {record['all_criteria_pass']}")
    print("wrote provenance/wp5_trend_replacement_validation_execution.json")


if __name__ == "__main__":
    main()
