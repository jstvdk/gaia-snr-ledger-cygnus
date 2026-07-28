#!/usr/bin/env python3
"""Issue #13: score the repair_v6 pre-registration against the outcome.

Predictions were written to provenance/wp5_node_interpolation_prereg.json
before any repair_v6 injection existed.  This script reads them back and
reports each one, with no discretion left at scoring time.

Adoption rests on P3 and P4, which bound the change without prescribing its
direction.  P2 -- the one regressing cell clearing -- is reported either way
and does NOT gate adoption, precisely because it is the cell that motivated
the investigation.

Outputs: provenance/wp5_node_interpolation_outcome.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp5_node_interpolation_outcome.py
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy import stats

import wp5_common as w
import wp5_residual_trend as T

BASELINE = ("PARSEC", 3.1, 2.3)
REGRESSING_CELL = ("CygOB2-C", "MIST", 3.5, 2.0)


def score_grid(version: str) -> pd.DataFrame:
    frame = pd.read_parquet(w.PROC / f"wp5_mass_function_bins_{version}.parquet")
    rows = []
    for (subgroup, family, rv, alpha), cell in frame.groupby(
        ["subgroup", "family", "R_V", "alpha"]
    ):
        ordered = cell.sort_values("bin_index")
        residual = ordered["pearson_residual"].to_numpy(float)
        log_mass = np.log10(ordered["mass_geometric_center"].to_numpy(float))
        expected = ordered["expected_count_at_k_median"].to_numpy(float)
        rate = ordered["completeness_weighted_imf_integral_per_k"].to_numpy(float)
        chi_p = float(stats.chi2.sf(float(np.sum(residual**2)), max(w.N_IMF_BINS - 1, 1)))
        max_abs = float(np.max(np.abs(residual)))
        incumbent = T.incumbent_trend_p(log_mass, residual)
        replacement, statistic = T.replacement_trend_p(residual, expected, rate, log_mass)
        rows.append(
            {
                "version": version, "subgroup": subgroup, "family": family,
                "R_V": float(rv), "alpha": float(alpha),
                "chi2_p": chi_p, "max_abs_residual": max_abs,
                "incumbent_trend_p": incumbent, "replacement_trend_p": replacement,
                "replacement_T": statistic,
                "residuals": [round(float(x), 3) for x in residual],
                "gate_incumbent": T.gate_pass(chi_p, incumbent, max_abs),
                "gate_replacement": T.gate_pass(chi_p, replacement, max_abs),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--version", default="repair_v6")
    parser.add_argument("--reference", default="repair_v5")
    args = parser.parse_args()
    version, reference = args.version, args.reference

    prereg = json.loads(
        (w.PROVENANCE / "wp5_node_interpolation_prereg.json").read_text(encoding="utf-8")
    )
    new = score_grid(version)
    old = score_grid(reference)

    # ---- P2: the regressing cell ------------------------------------------
    def cell(frame, key):
        subgroup, family, rv, alpha = key
        row = frame[
            frame.subgroup.eq(subgroup) & frame.family.eq(family)
            & frame.R_V.eq(rv) & frame.alpha.eq(alpha)
        ]
        return row.iloc[0] if len(row) == 1 else None

    before, after = cell(old, REGRESSING_CELL), cell(new, REGRESSING_CELL)
    p2 = {
        "id": "P2",
        "gates_adoption": False,
        "statement": prereg["predictions"][1]["statement"],
        f"{reference}_max_abs_residual": round(float(before.max_abs_residual), 3),
        f"{version}_max_abs_residual": round(float(after.max_abs_residual), 3),
        f"{reference}_residuals": before.residuals,
        f"{version}_residuals": after.residuals,
        f"{version}_trend_p": round(float(after.replacement_trend_p), 4),
        f"{version}_chi2_p": round(float(after.chi2_p), 4),
        f"{version}_gate_replacement": bool(after.gate_replacement),
        "confirmed": bool(after.max_abs_residual < 3.0),
    }

    # ---- P3: the baseline --------------------------------------------------
    baseline = new[
        new.family.eq(BASELINE[0]) & new.R_V.eq(BASELINE[1]) & new.alpha.eq(BASELINE[2])
    ]
    b_row = baseline[baseline.subgroup.eq("CygOB2-B")].iloc[0]
    p3 = {
        "id": "P3",
        "gates_adoption": True,
        "statement": prereg["predictions"][2]["statement"],
        "per_subgroup": [
            {
                "subgroup": row.subgroup,
                "residuals": row.residuals,
                "chi2_p": round(float(row.chi2_p), 4),
                "max_abs_residual": round(float(row.max_abs_residual), 3),
                "incumbent_trend_p": round(float(row.incumbent_trend_p), 4),
                "replacement_trend_p": round(float(row.replacement_trend_p), 4),
                "replacement_T": round(float(row.replacement_T), 3),
                "gate_incumbent": bool(row.gate_incumbent),
                "gate_replacement": bool(row.gate_replacement),
            }
            for row in baseline.itertuples()
        ],
        "all_pass_incumbent": bool(baseline.gate_incumbent.all()),
        "all_pass_replacement": bool(baseline.gate_replacement.all()),
        "B_abs_T": round(abs(float(b_row.replacement_T)), 3),
        "confirmed": bool(
            baseline.gate_incumbent.all()
            and baseline.gate_replacement.all()
            and abs(float(b_row.replacement_T)) < 2.0
        ),
    }

    # ---- P4: the change is local ------------------------------------------
    new_norm = pd.read_parquet(w.PROC / f"wp5_imf_normalization_{version}.parquet")
    old_norm = pd.read_parquet(w.PROC / f"wp5_imf_normalization_{reference}.parquet")
    keys = ["subgroup", "family", "R_V", "alpha"]
    merged = new_norm.merge(old_norm, on=keys, suffixes=("_new", "_old"))
    age_shift = (
        merged["truth_age_posterior_mean_Myr_new"]
        - merged["truth_age_posterior_mean_Myr_old"]
    ).abs()
    worst = merged.loc[age_shift.idxmax()]

    new_mass = pd.read_parquet(w.PROC / f"wp5_association_mass_{version}.parquet")
    old_mass = pd.read_parquet(w.PROC / f"wp5_association_mass_{reference}.parquet")

    def baseline_mass(frame):
        row = frame[
            frame.family.eq(BASELINE[0]) & frame.R_V.eq(BASELINE[1])
            & frame.alpha.eq(BASELINE[2])
        ]
        return float(row["multiplicity_adjusted_mass_median_Msun"].iloc[0])

    mass_new, mass_old = baseline_mass(new_mass), baseline_mass(old_mass)
    mass_change = (mass_new - mass_old) / mass_old
    grid_change = int(new.gate_replacement.sum()) - int(old.gate_replacement.sum())
    p4 = {
        "id": "P4",
        "gates_adoption": True,
        "statement": prereg["predictions"][3]["statement"],
        "max_abs_truth_age_posterior_mean_shift_Myr": round(float(age_shift.max()), 4),
        "max_shift_cell": (
            f"{worst.subgroup} {worst.family} R_V={worst.R_V} alpha={worst.alpha}"
        ),
        "median_abs_truth_age_posterior_mean_shift_Myr": round(
            float(age_shift.median()), 4
        ),
        "limit_truth_age_shift_Myr": 0.15,
        f"{reference}_baseline_association_mass_Msun": round(mass_old, 1),
        f"{version}_baseline_association_mass_Msun": round(mass_new, 1),
        "association_mass_fractional_change": round(float(mass_change), 4),
        "limit_association_mass_fractional_change": 0.05,
        "grid_replacement_before": int(old.gate_replacement.sum()),
        "grid_replacement_after": int(new.gate_replacement.sum()),
        "grid_incumbent_before": int(old.gate_incumbent.sum()),
        "grid_incumbent_after": int(new.gate_incumbent.sum()),
        "grid_change": grid_change,
        "limit_abs_grid_change": 6,
        "confirmed": bool(
            age_shift.max() < 0.15
            and abs(mass_change) < 0.05
            and abs(grid_change) <= 6
        ),
    }

    adopt = bool(p3["confirmed"] and p4["confirmed"])
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_node_interpolation_outcome.py",
        "status": "SUCCESS",
        "issue": "#13 — outcome of the pre-registered repair_v6 node-interpolation fix",
        "version": version,
        "reference": reference,
        "prereg": "provenance/wp5_node_interpolation_prereg.json",
        "prereg_sha256": w.sha256(w.PROVENANCE / "wp5_node_interpolation_prereg.json"),
        "adoption_rule": prereg["adoption_rule"],
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in [
                w.PROC / f"wp5_mass_function_bins_{version}.parquet",
                w.PROC / f"wp5_imf_normalization_{version}.parquet",
                w.PROC / f"wp5_association_mass_{version}.parquet",
            ]
        },
        "P1": {
            "id": "P1",
            "gates_adoption": True,
            "statement": prereg["predictions"][0]["statement"],
            "status": "established before the run",
            "evidence": "provenance/wp5_node_rule_continuity_execution.json",
            "confirmed": True,
        },
        "P2": p2,
        "P3": p3,
        "P4": p4,
        "adopted": adopt,
        "verdict": (
            "ADOPTED — P3 and P4 both hold" if adopt
            else "NOT ADOPTED — repair_v5 stands; the failure is a finding, "
                 "see the falsification clauses in the pre-registration"
        ),
    }
    w.write_json(w.PROVENANCE / "wp5_node_interpolation_outcome.json", record)

    print(f"pre-registration outcome for {version} against {reference}\n")
    print(f"  P1 continuity        {'CONFIRMED' if record['P1']['confirmed'] else 'FAILED'}"
          "   (established before the run)")
    print(f"  P2 regressing cell   {'CONFIRMED' if p2['confirmed'] else 'NOT CONFIRMED'}"
          f"   max|r| {p2[f'{reference}_max_abs_residual']:.2f} -> "
          f"{p2[f'{version}_max_abs_residual']:.2f}   (does not gate adoption)")
    print(f"  P3 baseline          {'CONFIRMED' if p3['confirmed'] else 'FAILED'}"
          f"   all-pass incumbent {p3['all_pass_incumbent']}, "
          f"replacement {p3['all_pass_replacement']}, B |T| {p3['B_abs_T']:.2f}")
    print(f"  P4 change is local   {'CONFIRMED' if p4['confirmed'] else 'FAILED'}"
          f"   age shift max {p4['max_abs_truth_age_posterior_mean_shift_Myr']:.4f} Myr, "
          f"mass {p4['association_mass_fractional_change']:+.2%}, "
          f"grid {p4['grid_replacement_before']} -> {p4['grid_replacement_after']}")
    print(f"\n{record['verdict']}")
    print("wrote provenance/wp5_node_interpolation_outcome.json")


if __name__ == "__main__":
    main()
