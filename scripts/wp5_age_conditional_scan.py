#!/usr/bin/env python3
"""Step 2 of the #1c gated plan: age-conditional refit of CygOB2-B (gate G2).

The 2026-07-27 shape diagnostic established that the WP5 injection truth model
generates synthetic photometry at a single age (the upper-MS MAP), making the
end-to-end design structurally blind to an age error (brief F3), and that B's
bump bin sits on the PMS/Henyey isochrone fold at B's MAP age (brief F2).

This scan re-runs B-only baseline injections (PARSEC, R_V=3.1) with the TRUTH
age monkeypatched to native-isochrone-age nodes spanning B's WP4 age
posterior, plus older nodes probing the direction of B's (unmeasurable,
grid-railed) PMS indicator.  The recovery side is untouched: synthetic sources
are still fit through the repaired WP3 posterior estimator and the WP4 mass
sampler that marginalizes the production 9-node age posterior.  Each node's
response is then fit with the unmodified official ``wp5_fit_imf.fit_one``.

Two posterior-marginalized mixtures are also fit, composed deterministically
from the per-node responses:

* ``posterior_9node`` — truth age marginalized over the WP4 posterior using
  the pipeline's own operational discretization (``age_posterior_nodes``: nine
  equiprobable split-normal nodes), each node snapped to its nearest native
  isochrone age.  This is the anti-tuning-compliant candidate fix.
* ``posterior_9node_sf1myr`` — the same posterior convolved with a 1 Myr
  uniform star-formation-duration spread (the WP7 branch the plan mandates).

Anti-tuning rule (brief section 4): the adopted repair_v4 fix must be the
posterior marginalization applied identically to all subgroups.  The
individual age nodes exist to *diagnose* the mechanism, never to pick a
best-scoring age.

No stored repair_v3 artifact is overwritten.  Node responses/curves are saved
as new snapshot parquets and hashed in the provenance record.

Run:
  WP_REPAIR_VERSION=repair_v3 WP3_ANCHOR_PRIOR_MODE=variogram \
  PYTHONPATH=scripts python3 scripts/wp5_age_conditional_scan.py
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import scipy
from scipy.stats import norm
from sklearn.isotonic import IsotonicRegression

import wp5_common as w
import wp5_fit_imf as F
import wp5_injections_repair as R
from wp3_repair_common import ANCHOR_PRIOR_MODE, REPAIR_VERSION, AnchorMap, load_template_library
from wp4_repair_common import age_posterior_nodes
from wp5_bump_shape_diagnostic import fold_interval

if REPAIR_VERSION != "repair_v3" or ANCHOR_PRIOR_MODE != "variogram":
    raise RuntimeError(
        "run with WP_REPAIR_VERSION=repair_v3 WP3_ANCHOR_PRIOR_MODE=variogram; "
        f"got {REPAIR_VERSION!r}/{ANCHOR_PRIOR_MODE!r}"
    )

SUBGROUP = "CygOB2-B"
FAMILY = "PARSEC"
RV = 3.1
BUMP_BIN = 2
# Native PARSEC isochrone ages (the truth generator snaps to them anyway).
# 2.818 = B's UMS MAP (the production truth; Monte-Carlo control node).
# 2.239/2.512 = younger probes below the posterior (2.512 ~ C's age).
# 3.162/3.548 = inside B's posterior (68%/90% upper tail).
# 3.981/4.467/5.012 = older probes in the PMS-indicator direction (B's PMS MAP
# is unmeasurable: n=2, grid-railed at 10 Myr; these stand in as direction
# probes; 4.467 also carries the 0.5% upper tail of the SF-spread mixture).
AGE_NODES = [2.239, 2.512, 2.818, 3.162, 3.548, 3.981, 4.467, 5.012]
CONTROL_NODE = 2.818


def node_tag(age: float) -> str:
    return f"{age:.3f}".replace(".", "p")


def snap_to_native(age: float, native_ages: np.ndarray) -> float:
    return float(native_ages[np.argmin(np.abs(native_ages - age))])


def split_normal_params(age_posterior: pd.DataFrame) -> tuple[float, float, float]:
    row = age_posterior[
        age_posterior["subgroup"].eq(SUBGROUP)
        & age_posterior["family"].eq(FAMILY)
        & age_posterior["R_V"].eq(RV)
        & age_posterior["f_bin"].eq(w.F_BINARY)
        & age_posterior["indicator"].eq("ums")
        & age_posterior["dmu"].eq(0.0)
    ].iloc[0]
    centre = float(row["age_map"])
    lower_scale = max(centre - float(row["age_lo68"]), 0.02) / abs(norm.ppf(0.16))
    upper_scale = max(float(row["age_hi68"]) - centre, 0.02) / norm.ppf(0.84)
    return centre, lower_scale, upper_scale


def nine_node_weights(
    age_posterior: pd.DataFrame, native_ages: np.ndarray
) -> dict[float, float]:
    """The pipeline's own posterior discretization, snapped to native ages."""
    nodes = age_posterior_nodes(age_posterior, SUBGROUP, FAMILY, RV)
    weights: dict[float, float] = {}
    for node in nodes:
        native = snap_to_native(float(node), native_ages)
        weights[native] = weights.get(native, 0.0) + 1.0 / len(nodes)
    return dict(sorted(weights.items()))


def sf_spread_weights(
    age_posterior: pd.DataFrame, native_ages: np.ndarray, duration_myr: float = 1.0
) -> dict[float, float]:
    """Split-normal posterior convolved with a U(+-duration/2) SF spread,
    integrated numerically and binned to nearest native isochrone age."""
    centre, lower_scale, upper_scale = split_normal_params(age_posterior)
    ages = np.linspace(1.0, 10.0, 18001)
    offsets = np.linspace(-duration_myr / 2.0, duration_myr / 2.0, 201)
    shifted = ages[None, :] - offsets[:, None]
    scale = np.where(shifted < centre, lower_scale, upper_scale)
    # split normal: common mode height, side-dependent scale
    density = np.exp(-0.5 * ((shifted - centre) / scale) ** 2) * (
        2.0 / (np.sqrt(2.0 * np.pi) * (lower_scale + upper_scale))
    )
    profile = density.mean(axis=0)
    profile /= np.trapezoid(profile, ages)
    nearest = np.array([snap_to_native(a, native_ages) for a in ages])
    weights: dict[float, float] = {}
    for native in np.unique(nearest):
        mask = nearest == native
        weight = float(np.trapezoid(np.where(mask, profile, 0.0), ages))
        # below one injection row per mass (1/400) a node cannot contribute
        if weight >= 1.0 / w.N_INJECT_PER_MASS:
            weights[float(native)] = weight
    total = sum(weights.values())
    return {k: v / total for k, v in sorted(weights.items())}


def largest_remainder(weights: dict[float, float], total: int) -> dict[float, int]:
    raw = {k: v * total for k, v in weights.items()}
    floors = {k: int(np.floor(v)) for k, v in raw.items()}
    shortfall = total - sum(floors.values())
    order = sorted(raw, key=lambda k: raw[k] - floors[k], reverse=True)
    for key in order[:shortfall]:
        floors[key] += 1
    return floors


def compose_mixture(
    node_responses: dict[float, pd.DataFrame], weights: dict[float, float]
) -> pd.DataFrame:
    """Deterministic posterior-weighted truth mixture: per true mass, take the
    first n_j rows from node j (rows within a mass are exchangeable draws)."""
    allocation = largest_remainder(weights, w.N_INJECT_PER_MASS)
    missing = [k for k, n in allocation.items() if n > 0 and k not in node_responses]
    if missing:
        raise RuntimeError(f"mixture needs unavailable age nodes {missing}")
    pieces = []
    for mass in w.MASS_GRID:
        for age, take in allocation.items():
            if take == 0:
                continue
            rows = node_responses[age]
            rows = rows[rows["true_primary_mass"].eq(mass)].head(take)
            if len(rows) != take:
                raise RuntimeError(f"node {age} has too few rows at mass {mass}")
            pieces.append(rows)
    return pd.concat(pieces, ignore_index=True)


def curve_from_response(response: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for mass in w.MASS_GRID:
        select = response[response["true_primary_mass"].eq(mass)]
        rows.append(
            {
                "subgroup": SUBGROUP,
                "family": FAMILY,
                "R_V": RV,
                "primary_mass": float(mass),
                "n_injected": int(len(select)),
                "recovery_fraction": float(
                    np.isfinite(select["recovered_mass"]).mean()
                ),
            }
        )
    curve = pd.DataFrame(rows)
    curve["recovery_isotonic"] = IsotonicRegression(
        y_min=0.0, y_max=1.0, increasing=True, out_of_bounds="clip"
    ).fit_transform(
        curve["primary_mass"], curve["recovery_fraction"], sample_weight=curve["n_injected"]
    )
    return curve


def fit_branch(
    masses: pd.DataFrame,
    curve: pd.DataFrame,
    response: pd.DataFrame,
    samples: np.ndarray,
    draw_cols: list[str],
) -> dict:
    out = {}
    rng = np.random.default_rng(w.SEED)
    for alpha in w.IMF_SLOPES:
        summary, bins, _ = F.fit_one(
            masses, curve, response, SUBGROUP, FAMILY, RV, alpha, rng,
            samples, draw_cols,
        )
        out[f"alpha_{alpha:.1f}"] = {
            "residuals": [float(r) for r in bins.sort_values("bin_index")["pearson_residual"]],
            "chi2_p": float(summary["poisson_chi_square_p"]),
            "trend_p": float(summary["residual_trend_p"]),
            "max_abs_residual": float(summary["max_abs_pearson_residual"]),
            "bump_bin_residual": float(
                bins.sort_values("bin_index")["pearson_residual"].iloc[BUMP_BIN]
            ),
            "k_median": float(summary["k_median"]),
            "residual_gate_pass": bool(summary["residual_gate_pass"]),
        }
    return out


def main() -> None:
    masses = pd.read_parquet(w.PROC / f"wp4_mass_posteriors_{REPAIR_VERSION}.parquet")
    sample_store = np.load(w.PROC / f"wp4_mass_posterior_samples_{REPAIR_VERSION}.npz")
    branch_samples = sample_store["samples"][
        :, w.FAMILIES.index(FAMILY) * len(w.R_V_BRANCHES) + list(w.R_V_BRANCHES).index(RV), :
    ]
    age_posterior = pd.read_parquet(w.PROC / f"wp4_age_posteriors_{REPAIR_VERSION}.parquet")
    iso_parsec = pd.read_parquet(w.PROC / "wp3_isochrones_parsec.parquet")
    native_ages = np.sort(iso_parsec["age_Myr"].unique())
    stored_norm = pd.read_parquet(w.PROC / f"wp5_imf_normalization_{REPAIR_VERSION}.parquet")
    stored_bins = pd.read_parquet(w.PROC / f"wp5_mass_function_bins_{REPAIR_VERSION}.parquet")

    # Snap the node list to the exact native grid values so node keys, mixture
    # weights and snapshot tags all refer to the same float.
    global AGE_NODES
    AGE_NODES = [snap_to_native(a, native_ages) for a in AGE_NODES]

    weights_posterior = nine_node_weights(age_posterior, native_ages)
    weights_sf = sf_spread_weights(age_posterior, native_ages, duration_myr=1.0)
    print("posterior 9-node weights:", weights_posterior, flush=True)
    print("posterior + 1 Myr SF spread weights:", weights_sf, flush=True)

    # --- heavy shared setup (identical to wp5_injections_repair.main) ---
    classifier = w.reconstruct_wp2_classifier()
    donor_pool, donor_model = R.build_donor_pool(classifier)
    donor_pool = R.augment_donor_pool(donor_pool)
    normal_points = R.sobol_normals(w.MEMBERSHIP_QMC_POINTS)
    validation = R.validate_qmc(classifier, normal_points)
    if validation["decision_agreement"] < 0.97:
        raise RuntimeError("QMC validation failed")
    posterior = np.load(w.PROC / f"wp3_extinction_posterior_{REPAIR_VERSION}.npz")
    posterior_ids = posterior["source_id"].astype("int64")
    posterior_cube = posterior["probability"]
    anchor_map = AnchorMap.from_frozen_wp3()
    _, template_magnitudes, template_weights = load_template_library()
    repair_provenance = json.loads(
        (w.PROVENANCE / "wp3_repair_execution.json").read_text(encoding="utf-8")
    )
    branch_sigma = float(
        repair_provenance["configuration"]["template_branch_uncertainty_calibration"][
            f"rv{RV:.1f}"
        ]["adopted_template_branch_sigma_mag"]
    )
    draw_cols = None

    original_loader = w.load_isochrone_at_age
    node_records = {}
    node_responses: dict[float, pd.DataFrame] = {}
    snapshots = {}
    try:
        for age_node in AGE_NODES:
            print(f"=== truth-age node {age_node} Myr ===", flush=True)

            def patched_loader(family: str, age_myr: float, _node=age_node):
                # Truth side only: force the node age.  The recovery side does
                # not call this function (wp4_repair_common has its own loader).
                return original_loader(family, _node)

            w.load_isochrone_at_age = patched_loader
            rng = np.random.default_rng(w.SEED)
            curve, response, summary = R.inject_curve(
                SUBGROUP, FAMILY, RV, classifier, donor_pool, donor_model,
                normal_points, rng, posterior_ids, posterior_cube, anchor_map,
                template_magnitudes, template_weights, branch_sigma,
                age_posterior,
            )
            w.load_isochrone_at_age = original_loader
            if draw_cols is None:
                draw_cols = sorted(
                    c for c in response.columns if c.startswith("recovered_mass_draw_")
                )
            tag = node_tag(age_node)
            response_path = (
                w.PROC / f"wp5_age_scan_B_response_age{tag}_{REPAIR_VERSION}.parquet"
            )
            curve_path = (
                w.PROC / f"wp5_age_scan_B_curve_age{tag}_{REPAIR_VERSION}.parquet"
            )
            response.to_parquet(response_path, index=False)
            curve.to_parquet(curve_path, index=False)
            snapshots[f"age_{tag}"] = {
                str(response_path.relative_to(w.ROOT)): w.sha256(response_path),
                str(curve_path.relative_to(w.ROOT)): w.sha256(curve_path),
            }
            node_responses[age_node] = response
            fits = fit_branch(masses, curve, response, branch_samples, draw_cols)
            node_records[f"{age_node:.3f}"] = {
                "truth_age_Myr": age_node,
                "native_isochrone_age_Myr": float(summary["age_isochrone_Myr"]),
                "posterior_weight_9node": weights_posterior.get(age_node, 0.0),
                "posterior_weight_sf1myr": weights_sf.get(age_node, 0.0),
                "in_wp4_central95": bool(2.750 <= age_node <= 3.421),
                "parsec_G0_fold_interval_Msun": fold_interval(iso_parsec, age_node),
                "injection_summary": {
                    k: v for k, v in summary.items()
                    if k in ["membership_pass", "mass_recovered", "member_donors",
                             "injected", "av_draw_median"]
                },
                "fits": fits,
            }
            print(json.dumps(node_records[f"{age_node:.3f}"]["fits"]["alpha_2.3"],
                             indent=2), flush=True)
    finally:
        w.load_isochrone_at_age = original_loader

    mixtures = {}
    for name, weights in [
        ("posterior_9node", weights_posterior),
        ("posterior_9node_sf1myr", weights_sf),
    ]:
        response = compose_mixture(node_responses, weights)
        curve = curve_from_response(response)
        response_path = (
            w.PROC / f"wp5_age_scan_B_response_{name}_{REPAIR_VERSION}.parquet"
        )
        response.to_parquet(response_path, index=False)
        snapshots[name] = {
            str(response_path.relative_to(w.ROOT)): w.sha256(response_path)
        }
        mixtures[name] = {
            "weights": {f"{k:.3f}": v for k, v in weights.items()},
            "rows_per_mass": largest_remainder(weights, w.N_INJECT_PER_MASS),
            "fits": fit_branch(masses, curve, response, branch_samples, draw_cols),
        }
        mixtures[name]["rows_per_mass"] = {
            f"{k:.3f}": v for k, v in mixtures[name]["rows_per_mass"].items()
        }
        print(f"=== mixture {name} ===", flush=True)
        print(json.dumps(mixtures[name]["fits"]["alpha_2.3"], indent=2), flush=True)

    stored_baseline = stored_norm[
        stored_norm.subgroup.eq(SUBGROUP) & stored_norm.family.eq(FAMILY)
        & stored_norm.R_V.eq(RV) & stored_norm.alpha.eq(2.3)
    ].iloc[0]
    stored_residuals = [
        float(r) for r in stored_bins[
            stored_bins.subgroup.eq(SUBGROUP) & stored_bins.family.eq(FAMILY)
            & stored_bins.R_V.eq(RV) & stored_bins.alpha.eq(2.3)
        ].sort_values("bin_index")["pearson_residual"]
    ]

    supported = [a for a in AGE_NODES if weights_posterior.get(a, 0.0) > 0.0]
    supported_pass = [
        a for a in supported
        if node_records[f"{a:.3f}"]["fits"]["alpha_2.3"]["residual_gate_pass"]
    ]
    mixture_pass = {
        name: mixtures[name]["fits"]["alpha_2.3"]["residual_gate_pass"]
        for name in mixtures
    }
    central95_flattens = [
        a for a in AGE_NODES
        if node_records[f"{a:.3f}"]["in_wp4_central95"]
        and node_records[f"{a:.3f}"]["fits"]["alpha_2.3"]["residual_gate_pass"]
    ]
    gate_g2 = {
        "criterion": (
            "B passes chi2_p>=0.01, trend_p>=0.05, max|res|<=3.0 at a truth age "
            "with non-negligible WP4 posterior support, or under the "
            "posterior-marginalized mixture"
        ),
        "nodes_with_posterior_support": [float(a) for a in supported],
        "supported_nodes_passing": [float(a) for a in supported_pass],
        "mixture_pass": mixture_pass,
        "central95_nodes_passing": [float(a) for a in central95_flattens],
        "mechanism_confirmed": bool(supported_pass or any(mixture_pass.values())),
        "verdict": (
            "PASS — age mechanism confirmed, proceed to step 3 (repair_v4 "
            "truth-side age marginalization for all subgroups)"
            if supported_pass or any(mixture_pass.values())
            else "FAIL — age excluded within the WP4 posterior, go to step 4 "
                 "contingencies"
        ),
    }

    inputs = [
        w.PROC / f"wp4_mass_posteriors_{REPAIR_VERSION}.parquet",
        w.PROC / f"wp4_mass_posterior_samples_{REPAIR_VERSION}.npz",
        w.PROC / f"wp4_age_posteriors_{REPAIR_VERSION}.parquet",
        w.PROC / f"wp3_extinction_{REPAIR_VERSION}.parquet",
        w.PROC / f"wp3_extinction_posterior_{REPAIR_VERSION}.npz",
        w.PROC / "wp3_isochrones_parsec.parquet",
        w.PROC / f"wp5_imf_normalization_{REPAIR_VERSION}.parquet",
        w.PROC / f"wp5_mass_function_bins_{REPAIR_VERSION}.parquet",
        w.PROVENANCE / "wp3_repair_execution.json",
    ]
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_age_conditional_scan.py",
        "status": "SUCCESS",
        "issue": "#1c step 2 — age-conditional refit of CygOB2-B (gate G2)",
        "repair_version": REPAIR_VERSION,
        "anchor_prior_mode": ANCHOR_PRIOR_MODE,
        "stored_artifacts_overwritten": False,
        "seed": w.SEED,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
        "method": {
            "branch": f"{SUBGROUP} / {FAMILY} / R_V={RV} (baseline)",
            "truth_age_monkeypatch": (
                "wp5_common.load_isochrone_at_age forced to the node age; used "
                "only by the truth generator.  Recovery side unchanged: WP3 "
                "posterior estimator + WP4 mass sampler marginalizing the "
                "production 9-node age posterior."
            ),
            "node_choice": (
                "native PARSEC isochrone ages (truth generator snaps to them): "
                "2.818 = UMS MAP control; 2.239/2.512 younger probes (and exact "
                "support for the 1 Myr SF-spread mixture); 3.162/3.548 inside "
                "the WP4 posterior; 3.981/5.012 PMS-direction probes.  B's PMS "
                "indicator MAP is unusable (n=2, grid-railed at 10 Myr, "
                "exclusion n_stars_below_15)."
            ),
            "per_node_rng": "fresh default_rng(SEED) per node for comparability",
            "fit": "unmodified wp5_fit_imf.fit_one, alpha in {2.0, 2.3, 2.6}",
            "mixture_composition": (
                "largest-remainder allocation of the 400 injections per true "
                "mass across nodes proportional to posterior weight; rows "
                "within a mass are exchangeable draws"
            ),
            "anti_tuning_rule": (
                "adopted fix must be marginalization over the existing WP4 age "
                "posterior applied identically to all subgroups; picking the "
                "single best-scoring age is forbidden (brief section 4)"
            ),
        },
        "qmc_validation": validation,
        "wp4_age_posterior_B_baseline": {
            "age_map": 2.818,
            "central68": [2.814, 3.246],
            "central95_from_lo90_hi90": [2.750, 3.421],
            "pms_indicator": "unmeasurable (n=2, grid_railed at 10 Myr)",
        },
        "stored_baseline_for_reference": {
            "residuals": stored_residuals,
            "chi2_p": float(stored_baseline.poisson_chi_square_p),
            "trend_p": float(stored_baseline.residual_trend_p),
            "max_abs_residual": float(stored_baseline.max_abs_pearson_residual),
        },
        "inputs": {str(p.relative_to(w.ROOT)): w.sha256(p) for p in inputs},
        "age_nodes": node_records,
        "mixtures": mixtures,
        "gate_G2": gate_g2,
        "response_snapshots": snapshots,
    }
    w.write_json(w.PROVENANCE / "wp5_age_scan_execution.json", record)
    print(json.dumps(gate_g2, indent=2))
    print("wrote provenance/wp5_age_scan_execution.json")


if __name__ == "__main__":
    main()
