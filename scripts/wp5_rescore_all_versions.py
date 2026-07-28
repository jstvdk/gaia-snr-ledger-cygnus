#!/usr/bin/env python3
"""R4: re-score every WP5 version under both trend statistics, equal footing.

Pre-declared in CUTS_AND_THRESHOLDS.md section 14.4.  A replacement diagnostic
that improves only the newest version is tuning in disguise, so every version
that exists on disk is re-scored with identical code and the grid change is
reported for all of them.

Only the trend statistic differs between the two scorings; chi-square, the
max-residual criterion, the 0.05 threshold and the three-way conjunction are
untouched, and every number is recomputed from each version's own stored
residuals so no version gets a different code path.

Output: provenance/wp5_rescore_all_versions_execution.json

Run:
  WP_REPAIR_VERSION=repair_v3 WP3_ANCHOR_PRIOR_MODE=variogram \
  PYTHONPATH=scripts python3 scripts/wp5_rescore_all_versions.py
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

VERSIONS = ["frozen", "repair_v1", "repair_v2", "repair_v3", "repair_v4", "repair_v5"]
BASELINE = ("PARSEC", 3.1, 2.3)


def bins_path(version: str):
    if version == "frozen":
        return w.PROC / "wp5_mass_function_bins.parquet"
    return w.PROC / f"wp5_mass_function_bins_{version}.parquet"


def score(version: str) -> tuple[pd.DataFrame, dict]:
    frame = pd.read_parquet(bins_path(version))
    rows = []
    for (subgroup, family, rv, alpha), cell in frame.groupby(
        ["subgroup", "family", "R_V", "alpha"]
    ):
        ordered = cell.sort_values("bin_index")
        residual = ordered["pearson_residual"].to_numpy(float)
        log_mass = np.log10(ordered["mass_geometric_center"].to_numpy(float))
        expected = ordered["expected_count_at_k_median"].to_numpy(float)
        rate = ordered["completeness_weighted_imf_integral_per_k"].to_numpy(float)
        chi_square = float(np.sum(residual**2))
        chi_p = float(stats.chi2.sf(chi_square, max(w.N_IMF_BINS - 1, 1)))
        max_abs = float(np.max(np.abs(residual)))
        incumbent_p = T.incumbent_trend_p(log_mass, residual)
        replacement_p, statistic = T.replacement_trend_p(
            residual, expected, rate, log_mass
        )
        rows.append(
            {
                "version": version, "subgroup": subgroup, "family": family,
                "R_V": float(rv), "alpha": float(alpha),
                "chi2_p": chi_p, "max_abs_residual": max_abs,
                "incumbent_trend_p": incumbent_p,
                "replacement_trend_p": replacement_p,
                "replacement_T": statistic,
                "gate_incumbent": T.gate_pass(chi_p, incumbent_p, max_abs),
                "gate_replacement": T.gate_pass(chi_p, replacement_p, max_abs),
            }
        )
    table = pd.DataFrame(rows)
    base = table[
        table.family.eq(BASELINE[0]) & table.R_V.eq(BASELINE[1])
        & table.alpha.eq(BASELINE[2])
    ]
    summary = {
        "version": version,
        "cells": int(len(table)),
        "grid_incumbent": int(table.gate_incumbent.sum()),
        "grid_replacement": int(table.gate_replacement.sum()),
        "grid_change": int(table.gate_replacement.sum() - table.gate_incumbent.sum()),
        "baseline_all_pass_incumbent": bool(base.gate_incumbent.all()),
        "baseline_all_pass_replacement": bool(base.gate_replacement.all()),
        "baseline_detail": [
            {
                "subgroup": row.subgroup,
                "chi2_p": row.chi2_p,
                "max_abs_residual": row.max_abs_residual,
                "incumbent_trend_p": row.incumbent_trend_p,
                "replacement_trend_p": row.replacement_trend_p,
                "gate_incumbent": bool(row.gate_incumbent),
                "gate_replacement": bool(row.gate_replacement),
            }
            for row in base.itertuples()
        ],
    }
    return table, summary


def main() -> None:
    tables, summaries = [], []
    for version in VERSIONS:
        if not bins_path(version).exists():
            print(f"  (skipping {version}: {bins_path(version).name} not on disk)")
            continue
        table, summary = score(version)
        tables.append(table)
        summaries.append(summary)
        print(
            f"{version:10s} grid {summary['grid_incumbent']:2d} -> "
            f"{summary['grid_replacement']:2d} of {summary['cells']} "
            f"({summary['grid_change']:+d})   baseline all-pass "
            f"{summary['baseline_all_pass_incumbent']} -> "
            f"{summary['baseline_all_pass_replacement']}"
        )
    combined = pd.concat(tables, ignore_index=True)
    out_csv = w.TABLES / "wp5_trend_statistic_rescore.csv"
    combined.to_csv(out_csv, index=False)

    changes = combined[combined.gate_incumbent != combined.gate_replacement]
    only_newest = bool(
        len(changes) > 0
        and set(changes["version"].unique()) == {VERSIONS[-1]}
    )
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_rescore_all_versions.py",
        "status": "SUCCESS",
        "issue": "#11 R4 — equal-footing re-scoring under the replacement statistic",
        "predeclaration": "CUTS_AND_THRESHOLDS.md section 14",
        "seed": T.BOOTSTRAP_SEED,
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "inputs": {
            str(bins_path(v).relative_to(w.ROOT)): w.sha256(bins_path(v))
            for v in VERSIONS if bins_path(v).exists()
        },
        "unchanged": (
            "chi-square criterion, max-residual criterion, the 0.05 trend "
            "threshold and the three-way conjunction are identical in both "
            "scorings; only the trend statistic differs"
        ),
        "per_version": summaries,
        "cells_changing_verdict": [
            {
                "version": row.version, "subgroup": row.subgroup,
                "family": row.family, "R_V": row.R_V, "alpha": row.alpha,
                "incumbent_trend_p": row.incumbent_trend_p,
                "replacement_trend_p": row.replacement_trend_p,
                "chi2_p": row.chi2_p,
                "max_abs_residual": row.max_abs_residual,
                "direction": (
                    "fail_to_pass" if row.gate_replacement else "pass_to_fail"
                ),
            }
            for row in changes.itertuples()
        ],
        "improves_only_newest_version": only_newest,
        "R4_pass": not only_newest,
        "outputs": {str(out_csv.relative_to(w.ROOT)): w.sha256(out_csv)},
    }
    w.write_json(w.PROVENANCE / "wp5_rescore_all_versions_execution.json", record)
    print(f"\n{len(changes)} of {len(combined)} cells change verdict, across versions "
          f"{sorted(changes['version'].unique())}")
    print(f"improves only the newest version: {only_newest}  -> R4 "
          f"{'PASS' if record['R4_pass'] else 'FAIL'}")
    print("wrote provenance/wp5_rescore_all_versions_execution.json")


if __name__ == "__main__":
    main()
