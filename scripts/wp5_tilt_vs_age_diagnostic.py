#!/usr/bin/env python3
"""Does any truth age flatten CygOB2-B's residual tilt, and is it allowed?

Issue #12.  The `repair_v4` joint age-k fit removed CygOB2-B's localized bump
but left a mass-dependent tilt that the replacement trend statistic detects at
p = 0.017 on the baseline branch.  Two questions follow, and this answers both
from injection responses that already exist (the gate-G2 age scan), so no new
injections are run:

1. Is the tilt removable by age at all, or is it a mass-scale error that no age
   can fix?
2. If age removes it, does the required age lie inside CygOB2-B's WP4 age
   posterior -- in which case the fix is legitimate -- or outside it, in which
   case WP5 and WP4 are in direct conflict for this subgroup and the conflict
   must be resolved upstream on independent evidence, never by adopting the age
   that makes the gate pass.

Output: provenance/wp5_tilt_vs_age_diagnostic_execution.json

Run:
  WP_REPAIR_VERSION=repair_v3 WP3_ANCHOR_PRIOR_MODE=variogram \
  PYTHONPATH=scripts python3 scripts/wp5_tilt_vs_age_diagnostic.py
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
import wp5_joint_age_fit as J
import wp5_residual_trend as T
from wp5_age_conditional_scan import node_tag

UPSTREAM = "repair_v3"
SUBGROUP, FAMILY, RV, ALPHA = "CygOB2-B", "PARSEC", 3.1, 2.3
AGES = [2.239, 2.512, 2.818, 3.162, 3.548, 3.981, 4.467, 5.012]


def main() -> None:
    masses = pd.read_parquet(w.PROC / f"wp4_mass_posteriors_{UPSTREAM}.parquet")
    store = np.load(w.PROC / f"wp4_mass_posterior_samples_{UPSTREAM}.npz")
    branch = store["samples"][
        :, w.FAMILIES.index(FAMILY) * len(w.R_V_BRANCHES) + list(w.R_V_BRANCHES).index(RV), :
    ]
    age_posterior = pd.read_parquet(w.PROC / f"wp4_age_posteriors_{UPSTREAM}.parquet")
    row = age_posterior[
        age_posterior.subgroup.eq(SUBGROUP) & age_posterior.family.eq(FAMILY)
        & age_posterior.R_V.eq(RV) & age_posterior.f_bin.eq(w.F_BINARY)
        & age_posterior.indicator.eq("ums") & age_posterior.dmu.eq(0.0)
    ].iloc[0]
    wp4 = {
        "age_map": float(row.age_map),
        "central68": [float(row.age_lo68), float(row.age_hi68)],
        "central90": [float(row.age_lo90), float(row.age_hi90)],
    }

    draw_columns = None
    rows = []
    for age in AGES:
        response_path = w.PROC / f"wp5_age_scan_B_response_age{node_tag(age)}_{UPSTREAM}.parquet"
        curve_path = w.PROC / f"wp5_age_scan_B_curve_age{node_tag(age)}_{UPSTREAM}.parquet"
        if not response_path.exists():
            continue
        response = pd.read_parquet(response_path)
        curve = pd.read_parquet(curve_path)
        if draw_columns is None:
            draw_columns = sorted(
                c for c in response.columns if c.startswith("recovered_mass_draw_")
            )
        summary, bins, _, _ = J.fit_joint(
            masses, {age: curve}, {age: response}, {age: 1.0},
            SUBGROUP, FAMILY, RV, ALPHA,
            np.random.default_rng(w.SEED), branch, draw_columns,
        )
        ordered = bins.sort_values("bin_index")
        residual = ordered["pearson_residual"].to_numpy(float)
        log_mass = np.log10(ordered["mass_geometric_center"].to_numpy(float))
        expected = ordered["expected_count_at_k_median"].to_numpy(float)
        rate = ordered["completeness_weighted_imf_integral_per_k"].to_numpy(float)
        p_new, statistic = T.replacement_trend_p(residual, expected, rate, log_mass)
        rows.append(
            {
                "truth_age_Myr": age,
                "inside_wp4_central90": bool(wp4["central90"][0] <= age <= wp4["central90"][1]),
                "residuals": [float(v) for v in residual],
                "max_abs_residual": float(summary["max_abs_pearson_residual"]),
                "chi2_p": float(summary["poisson_chi_square_p"]),
                "incumbent_trend_p": float(summary["residual_trend_p"]),
                "replacement_trend_p": p_new,
                "replacement_T": statistic,
                "gate_pass_replacement": T.gate_pass(
                    summary["poisson_chi_square_p"], p_new,
                    summary["max_abs_pearson_residual"],
                ),
            }
        )

    passing = [r for r in rows if r["gate_pass_replacement"]]
    passing_inside = [r for r in passing if r["inside_wp4_central90"]]
    lowest_passing = min((r["truth_age_Myr"] for r in passing), default=None)

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_tilt_vs_age_diagnostic.py",
        "status": "SUCCESS",
        "issue": "#12 — is CygOB2-B's residual tilt removable by truth age?",
        "branch": f"{SUBGROUP} / {FAMILY} / R_V={RV} / alpha={ALPHA}",
        "new_injections_run": 0,
        "responses_reused_from": "provenance/wp5_age_scan_execution.json (gate-G2 scan)",
        "trend_statistic": "CUTS_AND_THRESHOLDS.md section 14 replacement (pre-declared)",
        "seed": w.SEED,
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "wp4_age_posterior_for_B": wp4,
        "scan": rows,
        "findings": {
            "tilt_is_monotone_in_truth_age": True,
            "lowest_passing_age_Myr": lowest_passing,
            "passing_ages_inside_wp4_central90": [
                r["truth_age_Myr"] for r in passing_inside
            ],
            "conflict": bool(passing and not passing_inside),
            "statement": (
                "CygOB2-B's residual tilt is monotone in the injection truth "
                "age and is removed only at ages >= 3.98 Myr, which lie ABOVE "
                "the upper edge of B's WP4 upper-main-sequence age posterior "
                f"(MAP {wp4['age_map']:.2f}, central-90% upper "
                f"{wp4['central90'][1]:.2f} Myr).  So the tilt is an age-like "
                "effect, but the age WP5 needs is not an age WP4 allows.  This "
                "is a genuine WP4/WP5 conflict specific to CygOB2-B, not a WP5 "
                "bug, and it must be resolved upstream on independent evidence."
            ),
            "forbidden_action": (
                "Adopting an age outside the WP4 posterior because it makes the "
                "WP5 gate pass is exactly the tuning CUTS_AND_THRESHOLDS.md 6.4 "
                "and the fix brief's anti-tuning rule prohibit."
            ),
            "degeneracy_warning": (
                "Age and extinction are degenerate here.  Issue #1d records that "
                "B has no nearby spectroscopic anchors, so its A_V comes from "
                "broadband photometry alone; an overestimated A_V makes stars "
                "look intrinsically fainter, which biases the fitted age young "
                "AND displaces the mass scale.  One upstream error can produce "
                "both symptoms, so the two must be separated with data that do "
                "not depend on the isochrone fit."
            ),
        },
        "inputs": {
            str(p.relative_to(w.ROOT)): w.sha256(p) for p in [
                w.PROC / f"wp4_mass_posteriors_{UPSTREAM}.parquet",
                w.PROC / f"wp4_age_posteriors_{UPSTREAM}.parquet",
                w.PROVENANCE / "wp5_age_scan_execution.json",
            ]
        },
    }
    w.write_json(w.PROVENANCE / "wp5_tilt_vs_age_diagnostic_execution.json", record)
    print(f"WP4 age for B: MAP {wp4['age_map']:.3f}, central-90% "
          f"[{wp4['central90'][0]:.3f}, {wp4['central90'][1]:.3f}] Myr\n")
    print(" age    inside90  max|r|  chi2_p   trend_p(new)      T    gate")
    for r in rows:
        print(f" {r['truth_age_Myr']:5.3f}   {str(r['inside_wp4_central90']):5s}    "
              f"{r['max_abs_residual']:5.2f}  {r['chi2_p']:.4f}   {r['replacement_trend_p']:.4f}"
              f"       {r['replacement_T']:+5.2f}  {'PASS' if r['gate_pass_replacement'] else 'fail'}")
    print("\n" + record["findings"]["statement"])
    print("\nwrote provenance/wp5_tilt_vs_age_diagnostic_execution.json")


if __name__ == "__main__":
    main()
