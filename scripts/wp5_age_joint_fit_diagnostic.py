#!/usr/bin/env python3
"""Joint age--k fit forecast for CygOB2-B (step 2 addendum, informs step 3).

The gate-G2 scan (provenance/wp5_age_scan_execution.json) confirmed the age
mechanism: B's bump residual is strictly monotone in the injection truth age
and vanishes toward older ages, with two posterior-supported nodes (3.162 and
3.548 Myr) passing the full gate.  But the plain truth-side marginalization
over the WP4 posterior fails, because B's UMS posterior is extremely
bottom-heavy (sigma_lo ~ 0.02 Myr vs sigma_hi ~ 0.43 Myr) and puts 6/9 of its
weight on the MAP node.

The fix brief's anti-tuning rule (section 4) authorizes exactly one other
adoption: "a joint age--k fit with that [WP4] posterior as prior, applied
identically to all subgroups".  This diagnostic forecasts that fit for B on
the baseline branch, read-only over the saved node response snapshots:

* per age node j, the Poisson marginal likelihood of B's official observed
  soft-bin counts with k integrated out analytically under the Jeffreys
  k^(-1/2) prior:  log ML_j = sum_i n_i log r_ij - (N + 1/2) log R_j + const
  (n_i fixed across nodes, so the constant cancels in weight ratios);
* posterior node weights  w_j  proportional to  prior_j * ML_j  with the same
  9-node prior used by the recovery side (and the 1 Myr SF-spread prior as a
  variant);
* the posterior-predictive composite response (largest-remainder allocation of
  the 400 injections per true mass), fit with the unmodified official
  wp5_fit_imf.fit_one.

This is a diagnostic forecast of the repair_v4 joint fit, not an adoption run;
no stored artifact is overwritten.

Run:
  WP_REPAIR_VERSION=repair_v3 WP3_ANCHOR_PRIOR_MODE=variogram \
  PYTHONPATH=scripts python3 scripts/wp5_age_joint_fit_diagnostic.py
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import scipy

import wp5_common as w
import wp5_fit_imf as F
from wp5_age_conditional_scan import (
    SUBGROUP,
    FAMILY,
    RV,
    BUMP_BIN,
    compose_mixture,
    curve_from_response,
    fit_branch,
    node_tag,
)

REPAIR_VERSION = "repair_v3"
ALPHA_BASELINE = 2.3
BIN_EDGES = np.geomspace(w.CALIBRATION_NOMINAL_LO, w.CALIBRATION_HI, w.N_IMF_BINS + 1)


def rate_bins_for_response(response: pd.DataFrame, alpha: float) -> np.ndarray:
    """Point-estimate forward rates per k, the same math as fit_one."""
    draw_cols = sorted(
        c for c in response.columns if c.startswith("recovered_mass_draw_")
    )
    true_mass = np.sort(response["true_primary_mass"].unique())
    probability = np.zeros((len(true_mass), w.N_IMF_BINS))
    for index, value in enumerate(true_mass):
        rows = response[response["true_primary_mass"].eq(value)]
        active = [c for c in draw_cols if rows[c].notna().any()]
        if not active:
            continue
        draws = rows[active].to_numpy(float)
        for bin_index, (lo, hi) in enumerate(zip(BIN_EDGES[:-1], BIN_EDGES[1:], strict=True)):
            in_bin = (draws >= lo) & (
                draws < (hi if bin_index < w.N_IMF_BINS - 1 else hi + 1e-10)
            )
            probability[index, bin_index] = np.sum(in_bin) / (
                len(active) * len(rows)
            )
    weight = np.empty(len(true_mass))
    weight[0] = 0.5 * (true_mass[1] - true_mass[0])
    weight[-1] = 0.5 * (true_mass[-1] - true_mass[-2])
    weight[1:-1] = 0.5 * (true_mass[2:] - true_mass[:-2])
    return probability.T @ (weight * true_mass ** (-alpha))


def main() -> None:
    scan = json.loads(
        (w.PROVENANCE / "wp5_age_scan_execution.json").read_text(encoding="utf-8")
    )
    stored_bins = pd.read_parquet(
        w.PROC / f"wp5_mass_function_bins_{REPAIR_VERSION}.parquet"
    )
    observed = (
        stored_bins[
            stored_bins.subgroup.eq(SUBGROUP)
            & stored_bins.family.eq(FAMILY)
            & stored_bins.R_V.eq(RV)
            & stored_bins.alpha.eq(ALPHA_BASELINE)
        ]
        .sort_values("bin_index")["membership_weighted_count"]
        .to_numpy(float)
    )
    total = float(observed.sum())

    masses = pd.read_parquet(w.PROC / f"wp4_mass_posteriors_{REPAIR_VERSION}.parquet")
    sample_store = np.load(w.PROC / f"wp4_mass_posterior_samples_{REPAIR_VERSION}.npz")
    branch_samples = sample_store["samples"][
        :, w.FAMILIES.index(FAMILY) * len(w.R_V_BRANCHES) + list(w.R_V_BRANCHES).index(RV), :
    ]

    node_responses: dict[float, pd.DataFrame] = {}
    node_paths = {}
    for key, record in scan["age_nodes"].items():
        age = float(record["truth_age_Myr"])
        path = (
            w.PROC
            / f"wp5_age_scan_B_response_age{node_tag(age)}_{REPAIR_VERSION}.parquet"
        )
        node_responses[age] = pd.read_parquet(path)
        node_paths[key] = path

    log_ml = {}
    for age, response in node_responses.items():
        rates = rate_bins_for_response(response, ALPHA_BASELINE)
        log_ml[age] = float(
            np.sum(observed * np.log(rates)) - (total + 0.5) * np.log(rates.sum())
        )

    draw_cols = sorted(
        c
        for c in next(iter(node_responses.values())).columns
        if c.startswith("recovered_mass_draw_")
    )

    results = {}
    for prior_name, weight_key in [
        ("posterior_9node", "posterior_weight_9node"),
        ("posterior_9node_sf1myr", "posterior_weight_sf1myr"),
    ]:
        prior = {
            float(rec["truth_age_Myr"]): float(rec[weight_key])
            for rec in scan["age_nodes"].values()
            if float(rec[weight_key]) > 0.0
        }
        log_post = {
            age: np.log(pw) + log_ml[age] for age, pw in prior.items()
        }
        peak = max(log_post.values())
        posterior = {age: float(np.exp(v - peak)) for age, v in log_post.items()}
        norm = sum(posterior.values())
        posterior = {age: v / norm for age, v in posterior.items()}
        usable = {
            age: v for age, v in posterior.items() if v >= 1.0 / w.N_INJECT_PER_MASS
        }
        renorm = sum(usable.values())
        usable = {age: v / renorm for age, v in sorted(usable.items())}

        response = compose_mixture(node_responses, usable)
        curve = curve_from_response(response)
        fits = fit_branch(masses, curve, response, branch_samples, draw_cols)
        results[prior_name] = {
            "prior_weights": {f"{a:.3f}": v for a, v in sorted(prior.items())},
            "posterior_weights": {f"{a:.3f}": v for a, v in sorted(posterior.items())},
            "composited_weights_after_1_per_400_floor": {
                f"{a:.3f}": v for a, v in usable.items()
            },
            "fits": fits,
            "baseline_gate_pass": bool(
                fits[f"alpha_{ALPHA_BASELINE:.1f}"]["residual_gate_pass"]
            ),
        }
        print(f"=== joint fit, prior {prior_name} ===")
        print("posterior weights:", results[prior_name]["posterior_weights"])
        print(json.dumps(fits[f"alpha_{ALPHA_BASELINE:.1f}"], indent=2))

    forecast = {
        name: r["baseline_gate_pass"] for name, r in results.items()
    }
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_age_joint_fit_diagnostic.py",
        "status": "SUCCESS",
        "issue": "#1c step 2 addendum — joint age–k fit forecast for step 3",
        "repair_version": REPAIR_VERSION,
        "stored_artifacts_overwritten": False,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
        "method": {
            "marginal_likelihood": (
                "Poisson bins with Jeffreys k prior integrated analytically: "
                "log ML_j = sum_i n_i log r_ij - (N+1/2) log R_j + const; "
                "n_i are the official repair_v3 soft-bin membership-weighted "
                "counts (alpha=2.3 baseline), r_ij the point-estimate forward "
                "rates of node j's saved response snapshot"
            ),
            "posterior": "prior (recovery-side age discretization) times ML, normalized",
            "posterior_predictive": (
                "largest-remainder composite of the node responses at the "
                "posterior weights, fit with unmodified wp5_fit_imf.fit_one"
            ),
            "anti_tuning_status": (
                "the age posterior update is driven by the WP5 likelihood "
                "under the WP4 prior — the joint-fit adoption the brief "
                "section 4 explicitly authorizes; no gate statistic enters "
                "the weight computation"
            ),
        },
        "observed_counts_baseline": [float(v) for v in observed],
        "log_marginal_likelihood_by_node": {
            f"{a:.3f}": v for a, v in sorted(log_ml.items())
        },
        "inputs": {
            **{
                str(p.relative_to(w.ROOT)): w.sha256(p)
                for p in [
                    w.PROVENANCE / "wp5_age_scan_execution.json",
                    w.PROC / f"wp5_mass_function_bins_{REPAIR_VERSION}.parquet",
                    w.PROC / f"wp4_mass_posteriors_{REPAIR_VERSION}.parquet",
                    w.PROC / f"wp4_mass_posterior_samples_{REPAIR_VERSION}.npz",
                ]
            },
            **{
                str(p.relative_to(w.ROOT)): w.sha256(p)
                for p in node_paths.values()
            },
        },
        "results": results,
        "step3_forecast": {
            "joint_fit_baseline_gate_pass": forecast,
            "implication": (
                "if any prior variant passes, repair_v4 should implement the "
                "joint age–k fit (WP4 posterior as prior) uniformly for all "
                "three subgroups; if none passes, step 4 contingencies apply"
            ),
        },
    }
    w.write_json(
        w.PROVENANCE / "wp5_age_joint_fit_diagnostic_execution.json", record
    )
    print(json.dumps(record["step3_forecast"], indent=2))
    print("wrote provenance/wp5_age_joint_fit_diagnostic_execution.json")


if __name__ == "__main__":
    main()
