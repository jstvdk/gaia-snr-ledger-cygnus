#!/usr/bin/env python3
"""Rerun preserved WP5 injections through the actual repaired WP3/WP4 estimators.

Selection, donor cloning, mass grid, injection count, WP2 classifier and gate
thresholds are inherited unchanged from ``wp5_injections.py``.  The repaired
path changes only the inverse estimator: every recovered synthetic source is
fit by the versioned WP3 full-posterior estimator and then by the repaired WP4
mass-posterior propagation.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.isotonic import IsotonicRegression

import wp5_common as w
from wp2_finalize_membership import TARGET_CENTER, measurement_covariances
from wp2_membership_pipeline import FEATURES
from wp3_common import BANDS, gaia_mag_error
from wp3_extinction_law import band_coefficients
from wp3_repair_common import (
    REPAIR_VERSION,
    AnchorMap,
    fit_extinction_posterior,
    load_template_library,
)
from wp4_repair_common import age_posterior_nodes, infer_mass_samples
from wp5_injections import (
    build_donor_pool,
    choose_observational_donors,
    sobol_normals,
    validate_qmc,
    wilson_interval,
)


MASS_POSTERIOR_DRAWS_INJECTION = 64


def augment_donor_pool(donor_pool: pd.DataFrame) -> pd.DataFrame:
    twomass = pd.read_parquet(
        w.PROC / "wp1_2mass_join.parquet",
        columns=[
            "source_id",
            "j_msigcom",
            "h_msigcom",
            "ks_msigcom",
            "has_j",
            "has_h",
            "has_ks",
        ],
    )
    return donor_pool.merge(twomass, on="source_id", how="left", validate="one_to_one")


def load_member_donors(subgroup: str, rv: float) -> pd.DataFrame:
    members = pd.read_parquet(w.PROC / "wp2_members.parquet")
    labels = pd.read_parquet(w.TABLES / "wp2_subgroup_labels.parquet")
    extinction = pd.read_parquet(
        w.PROC / f"wp3_extinction_{REPAIR_VERSION}.parquet"
    )
    base = members[
        members["membership_probability"].gt(0.5)
        & ~members["anchor_quality_exempt"]
    ].drop(columns=["subgroup"], errors="ignore")
    base = base.merge(labels[["source_id", "subgroup"]], on="source_id", how="inner")
    base = base.merge(
        extinction[
            ["source_id", f"av_rv{rv:.1f}", f"av_err_rv{rv:.1f}"]
        ],
        on="source_id",
        how="left",
        validate="one_to_one",
    )
    base = base[
        base["subgroup"].eq(subgroup)
        & base[f"av_rv{rv:.1f}"].notna()
        & base[["l_deg", "b_deg", "parallax_corrected", "pmra", "pmdec"]]
        .notna()
        .all(axis=1)
    ].reset_index(drop=True)
    if len(base) < 50:
        raise RuntimeError(f"too few repair injection donors for {subgroup}/R_V={rv}")
    return base


def draw_donor_extinction(
    source_id: np.ndarray,
    rv: float,
    posterior_ids: np.ndarray,
    posterior_cube: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    lookup = {int(value): index for index, value in enumerate(posterior_ids)}
    rv_index = list(w.R_V_BRANCHES).index(rv)
    output = np.empty(len(source_id), dtype=float)
    for index, value in enumerate(source_id):
        probability = posterior_cube[lookup[int(value)], rv_index].astype(float)
        probability /= probability.sum()
        output[index] = np.interp(
            rng.random(), np.cumsum(probability), np.linspace(0.0, 15.0, 301)
        )
    return output


def infer_repaired_inverse(
    subgroup: str,
    family: str,
    rv: float,
    source_seed: np.ndarray,
    longitude: np.ndarray,
    latitude: np.ndarray,
    observed_magnitude: dict[str, np.ndarray],
    observed_error: dict[str, np.ndarray],
    recovered: np.ndarray,
    anchor_map: AnchorMap,
    template_magnitudes: np.ndarray,
    template_weights: np.ndarray,
    branch_sigma: float,
    age_nodes: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n_row = len(longitude)
    av_median = np.full(n_row, np.nan)
    mass_median = np.full(n_row, np.nan)
    mass_q16 = np.full(n_row, np.nan)
    mass_q84 = np.full(n_row, np.nan)
    mass_draws = np.full((n_row, MASS_POSTERIOR_DRAWS_INJECTION), np.nan)
    prior_mean, prior_separation = anchor_map.evaluate(longitude, latitude, rv)
    # Per-star prior width, so injected synthetic stars are fit under exactly
    # the same prior their real neighbours received (ANCHOR_PRIOR_MODE).
    prior_width = anchor_map.prior_sigma_at(prior_separation, rv)
    selected = np.flatnonzero(recovered)
    for progress, index in enumerate(selected, start=1):
        magnitudes = np.array(
            [observed_magnitude[band][index] for band in BANDS], dtype=float
        )
        errors = np.array(
            [observed_error[band][index] for band in BANDS], dtype=float
        )
        probability, summary = fit_extinction_posterior(
            magnitudes,
            errors,
            rv,
            float(prior_mean[index]),
            float(prior_width[index]),
            template_magnitudes,
            template_weights,
            template_branch_sigma=branch_sigma,
        )
        if probability is None:
            continue
        av_median[index] = float(summary["av_q50"])
        row = pd.Series(
            {
                "source_id": int(source_seed[index]),
                **{
                    band: observed_magnitude[band][index]
                    for band in BANDS
                },
                **{
                    f"{band}_err": observed_error[band][index]
                    for band in BANDS
                },
            }
        )
        samples = infer_mass_samples(
            row,
            probability,
            family,
            rv,
            age_nodes,
            n_draws=MASS_POSTERIOR_DRAWS_INJECTION,
            seed_offset=800_000,
        )
        finite = samples[np.isfinite(samples)]
        mass_draws[index] = samples
        if len(finite):
            mass_q16[index], mass_median[index], mass_q84[index] = np.quantile(
                finite, [0.16, 0.50, 0.84]
            )
        if progress % 2500 == 0:
            print(
                f"  repaired inverse {progress}/{len(selected)} recovered sources",
                flush=True,
            )
    return av_median, mass_median, mass_q16, mass_q84, mass_draws


def inject_curve(
    subgroup: str,
    family: str,
    rv: float,
    classifier: w.WP2Classifier,
    donor_pool: pd.DataFrame,
    donor_model,
    normal_points: np.ndarray,
    rng: np.random.Generator,
    posterior_ids: np.ndarray,
    posterior_cube: np.ndarray,
    anchor_map: AnchorMap,
    template_magnitudes: np.ndarray,
    template_weights: np.ndarray,
    branch_sigma: float,
    age_posterior: pd.DataFrame,
    truth_age_override: float | None = None,
    interpolate_truth_age: bool = False,
    mass_grid: np.ndarray | None = None,
    truth_binary_fraction=None,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Inject one branch.

    ``truth_age_override`` replaces the upper-MS MAP age used to generate the
    injected *truth* photometry; the recovery side keeps marginalizing the full
    WP4 age posterior either way.  It is consumed before any RNG draw, so a
    node family generated with the same seed differs only in truth age.  The
    default None reproduces the frozen single-age repair_v1..v3 behaviour
    exactly.

    ``interpolate_truth_age`` builds the truth isochrone by interpolating
    between native ages instead of snapping to the nearest one (issue #13),
    which is what the recovery side already does.  It changes no RNG draw, so
    the default False likewise reproduces repair_v1..v5 exactly.

    ``mass_grid`` replaces ``wp5_common.MASS_GRID`` as the set of true primary
    masses to inject.  WP6 uses it to extend the response above the frozen
    18 Msun ceiling (see provenance/wp6_mass_extension_decision.json) without
    touching MASS_GRID itself, which the accepted repair_v6 fit filters on.
    The default None is the frozen grid.

    ``truth_binary_fraction`` is a callable f(mass) -> probability replacing the
    constant ``wp5_common.F_BINARY`` on the TRUTH side only (issue #15).  The
    recovery side is deliberately left alone: the bias under test is exactly
    the mismatch between the rate nature makes binaries at and the 0.40 the
    estimator assumes.  The per-star threshold consumes the same single
    ``rng.random(n_injected)`` draw as the constant version, so donor,
    extinction and QMC realizations are bit-identical to the parent node.  The
    default None reproduces the frozen behaviour exactly.
    """

    age_row = age_posterior[
        age_posterior["subgroup"].eq(subgroup)
        & age_posterior["family"].eq(family)
        & age_posterior["R_V"].eq(rv)
        & age_posterior["f_bin"].eq(w.F_BINARY)
        & age_posterior["indicator"].eq("ums")
        & age_posterior["dmu"].eq(0.0)
    ]
    if len(age_row) != 1:
        raise RuntimeError(
            f"missing repaired age for {subgroup}/{family}/R_V={rv}"
        )
    age = float(age_row["age_map"].iloc[0])
    if truth_age_override is not None:
        age = float(truth_age_override)
    if interpolate_truth_age:
        iso, native_age = w.load_isochrone_between_ages(family, age)
    else:
        iso, native_age = w.load_isochrone_at_age(family, age)
    member_donors = load_member_donors(subgroup, rv)
    p_weights = member_donors["membership_probability"].to_numpy(float, copy=True)
    p_weights /= p_weights.sum()

    grid = w.MASS_GRID if mass_grid is None else np.asarray(mass_grid, dtype=float)
    primary_mass = np.repeat(grid, w.N_INJECT_PER_MASS)
    n_injected = len(primary_mass)
    donor_index = rng.choice(
        len(member_donors), size=n_injected, replace=True, p=p_weights
    )
    member = member_donors.iloc[donor_index].reset_index(drop=True)
    binary_threshold = (
        np.full(n_injected, float(w.F_BINARY))
        if truth_binary_fraction is None
        else np.asarray(truth_binary_fraction(primary_mass), dtype=float)
    )
    if binary_threshold.shape != (n_injected,):
        raise RuntimeError(
            "truth_binary_fraction must return one probability per injected star"
        )
    binary = rng.random(n_injected) < binary_threshold
    q = rng.uniform(w.Q_MIN, 1.0, n_injected)
    secondary_mass = np.where(binary, q * primary_mass, 0.0)
    absolute = w.interpolate_photometry(iso, primary_mass, secondary_mass)
    av_true = draw_donor_extinction(
        member["source_id"].to_numpy("int64"),
        rv,
        posterior_ids,
        posterior_cube,
        rng,
    )
    coefficients = band_coefficients(rv)
    apparent = {
        band: absolute[band] + w.DIST_MODULUS + coefficients[band] * av_true
        for band in absolute
    }
    colour = apparent["BP"] - apparent["RP"]
    observational_index = choose_observational_donors(
        donor_pool,
        donor_model,
        member["l_deg"].to_numpy(float),
        member["b_deg"].to_numpy(float),
        apparent["G"],
        colour,
        rng,
    )
    observed = donor_pool.iloc[observational_index].reset_index(drop=True)

    observed_error: dict[str, np.ndarray] = {}
    observed_magnitude: dict[str, np.ndarray] = {}
    for band, donor_mag, donor_flux_error in [
        ("G", "phot_g_mean_mag", "phot_g_mean_flux_error"),
        ("BP", "phot_bp_mean_mag", "phot_bp_mean_flux_error"),
        ("RP", "phot_rp_mean_mag", "phot_rp_mean_flux_error"),
    ]:
        observed_error[band] = gaia_mag_error(
            observed[donor_mag].to_numpy(float),
            observed[donor_flux_error].to_numpy(float),
            band,
        )
        observed_magnitude[band] = apparent[band] + rng.normal(
            0.0, observed_error[band]
        )
    for band, error_column, availability_column in [
        ("J", "j_msigcom", "has_j"),
        ("H", "h_msigcom", "has_h"),
        ("Ks", "ks_msigcom", "has_ks"),
    ]:
        error = observed[error_column].to_numpy(float)
        available = observed[availability_column].fillna(False).to_numpy(bool)
        valid_error = available & np.isfinite(error) & (error > 0)
        observed_error[band] = np.where(valid_error, error, np.nan)
        magnitude = apparent[band] + rng.normal(
            0.0, np.where(valid_error, error, 0.0)
        )
        observed_magnitude[band] = np.where(valid_error, magnitude, np.nan)

    covariance_scaled, covariance_repairs = measurement_covariances(
        observed, classifier.scaler.scale_
    )
    covariance_physical = covariance_scaled.copy()
    covariance_physical *= classifier.scaler.scale_[None, :, None]
    covariance_physical *= classifier.scaler.scale_[None, None, :]
    cholesky_physical = np.linalg.cholesky(covariance_physical)
    central_true = member[FEATURES].to_numpy(float)
    central_measured = central_true + np.einsum(
        "ni,nji->nj", rng.standard_normal((n_injected, 3)), cholesky_physical
    )
    raw_parallax = (
        central_measured[:, 0] + observed["parallax_zero_point"].to_numpy(float)
    )
    query_pass = (
        (observed_magnitude["G"] < 19.0)
        & (raw_parallax > 0.35)
        & (raw_parallax < 1.10)
    )
    quality_pass = query_pass & observed["quality_pass"].to_numpy(bool)
    cmd_ready = quality_pass & observed["has_cmd_photometry"].to_numpy(bool)
    membership_probability = np.zeros(n_injected, dtype=float)
    evaluated = np.flatnonzero(cmd_ready)
    if len(evaluated):
        x_scaled = classifier.scaler.transform(central_measured[evaluated])
        offsets = np.column_stack(
            [
                (
                    member["l_deg"].to_numpy(float)[evaluated] - TARGET_CENTER[0]
                )
                * np.cos(np.deg2rad(TARGET_CENTER[1])),
                member["b_deg"].to_numpy(float)[evaluated] - TARGET_CENTER[1],
            ]
        )
        membership_probability[evaluated] = w.qmc_membership_probabilities(
            x_scaled,
            offsets,
            covariance_scaled[evaluated],
            classifier,
            normal_points,
        )
    membership_pass = cmd_ready & (
        membership_probability > w.MEMBERSHIP_THRESHOLD
    )
    recovered = membership_pass

    nodes = age_posterior_nodes(age_posterior, subgroup, family, rv)
    synthetic_id = (
        np.arange(n_injected, dtype=np.int64)
        + int(round(rv * 10)) * 10_000_000
        + (0 if family == "PARSEC" else 100_000_000)
        + w.SUBGROUPS.index(subgroup) * 1_000_000
    )
    (
        av_estimate,
        recovered_mass,
        recovered_q16,
        recovered_q84,
        recovered_draws,
    ) = infer_repaired_inverse(
        subgroup,
        family,
        rv,
        synthetic_id,
        member["l_deg"].to_numpy(float),
        member["b_deg"].to_numpy(float),
        observed_magnitude,
        observed_error,
        recovered,
        anchor_map,
        template_magnitudes,
        template_weights,
        branch_sigma,
        nodes,
    )
    recovered_mass[~recovered] = np.nan

    rows = []
    for mass in grid:
        select = primary_mass == mass
        n = int(select.sum())
        stages = {
            "n_query": int(query_pass[select].sum()),
            "n_quality": int(quality_pass[select].sum()),
            "n_cmd_ready": int(cmd_ready[select].sum()),
            "n_membership": int(membership_pass[select].sum()),
            "n_recovered": int(np.isfinite(recovered_mass[select]).sum()),
        }
        low, high = wilson_interval(np.array(stages["n_recovered"]), n)
        rows.append(
            {
                "subgroup": subgroup,
                "family": family,
                "R_V": rv,
                "age_fit_Myr": age,
                "age_isochrone_Myr": native_age,
                "primary_mass": mass,
                "n_injected": n,
                **stages,
                "query_fraction": stages["n_query"] / n,
                "quality_fraction": stages["n_quality"] / n,
                "cmd_ready_fraction": stages["n_cmd_ready"] / n,
                "membership_fraction": stages["n_membership"] / n,
                "recovery_fraction": stages["n_recovered"] / n,
                "recovery_wilson95_lo": float(low),
                "recovery_wilson95_hi": float(high),
            }
        )
    curve = pd.DataFrame(rows)
    curve["recovery_isotonic"] = IsotonicRegression(
        y_min=0.0, y_max=1.0, increasing=True, out_of_bounds="clip"
    ).fit_transform(
        curve["primary_mass"],
        curve["recovery_fraction"],
        sample_weight=curve["n_injected"],
    )
    response = pd.DataFrame(
        {
            "subgroup": subgroup,
            "family": family,
            "R_V": rv,
            "true_primary_mass": primary_mass,
            "is_unresolved_binary": binary,
            "mass_ratio": np.where(binary, q, np.nan),
            "av_true": av_true,
            "av_estimate": av_estimate,
            "G_observed": observed_magnitude["G"],
            "BP_RP_observed": observed_magnitude["BP"] - observed_magnitude["RP"],
            "has_complete_jhk": (
                np.isfinite(observed_magnitude["J"])
                & np.isfinite(observed_magnitude["H"])
                & np.isfinite(observed_magnitude["Ks"])
            ),
            "query_pass": query_pass,
            "quality_pass": quality_pass,
            "cmd_ready": cmd_ready,
            "membership_probability_qmc": membership_probability,
            "membership_pass": membership_pass,
            "recovered_mass": recovered_mass,
            "recovered_mass_q16": recovered_q16,
            "recovered_mass_q84": recovered_q84,
        }
    )
    for draw_index in range(MASS_POSTERIOR_DRAWS_INJECTION):
        response[f"recovered_mass_draw_{draw_index:02d}"] = recovered_draws[
            :, draw_index
        ]
    summary = {
        "subgroup": subgroup,
        "family": family,
        "R_V": rv,
        "age_fit_Myr": age,
        "age_isochrone_Myr": native_age,
        "injected": n_injected,
        "member_donors": int(len(member_donors)),
        "observational_donors_used": int(len(np.unique(observational_index))),
        "covariance_repairs": int(covariance_repairs),
        "binary_fraction_realized": float(binary.mean()),
        "binary_fraction_truth_model": (
            "constant F_BINARY"
            if truth_binary_fraction is None
            else "mass-dependent (issue #15)"
        ),
        "binary_fraction_truth_by_mass": {
            f"{mass:g}": round(
                float(binary_threshold[primary_mass == mass][0]), 4
            )
            for mass in grid
        },
        "av_draw_median": float(np.median(av_true)),
        "g_mag_range": [float(np.min(apparent["G"])), float(np.max(apparent["G"]))],
        "membership_pass": int(recovered.sum()),
        "mass_recovered": int(np.isfinite(recovered_mass).sum()),
        "inverse_mass_draws": MASS_POSTERIOR_DRAWS_INJECTION,
    }
    return curve, response, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--baseline-only",
        action="store_true",
        help="replace only PARSEC R_V=3.1 branches in existing repair outputs",
    )
    parser.add_argument(
        "--output-version",
        default=REPAIR_VERSION,
        help=(
            "suffix for the WP5 injection outputs and provenance record.  "
            "Upstream WP3/WP4 inputs always stay at REPAIR_VERSION; use this "
            "when a WP5-only re-run must not overwrite an earlier version."
        ),
    )
    args = parser.parse_args()
    output_version = args.output_version
    classifier = w.reconstruct_wp2_classifier()
    donor_pool, donor_model = build_donor_pool(classifier)
    donor_pool = augment_donor_pool(donor_pool)
    normal_points = sobol_normals(w.MEMBERSHIP_QMC_POINTS)
    validation = validate_qmc(classifier, normal_points)
    if validation["decision_agreement"] < 0.97:
        raise RuntimeError("repair injection QMC validation failed")

    posterior = np.load(
        w.PROC / f"wp3_extinction_posterior_{REPAIR_VERSION}.npz"
    )
    posterior_ids = posterior["source_id"].astype("int64")
    posterior_cube = posterior["probability"]
    anchor_map = AnchorMap.from_frozen_wp3()
    _, template_magnitudes, template_weights = load_template_library()
    repair_provenance = json.loads(
        (w.PROVENANCE / "wp3_repair_execution.json").read_text(encoding="utf-8")
    )
    branch_sigma = {
        rv: float(
            repair_provenance["configuration"][
                "template_branch_uncertainty_calibration"
            ][f"rv{rv:.1f}"]["adopted_template_branch_sigma_mag"]
        )
        for rv in w.R_V_BRANCHES
    }
    age_posterior = pd.read_parquet(
        w.PROC / f"wp4_age_posteriors_{REPAIR_VERSION}.parquet"
    )
    rng = np.random.default_rng(w.SEED)
    curves = []
    responses = []
    summaries = []
    families = ["PARSEC"] if args.baseline_only else w.FAMILIES
    rv_branches = [3.1] if args.baseline_only else w.R_V_BRANCHES
    for family in families:
        for rv in rv_branches:
            for subgroup in w.SUBGROUPS:
                print(f"repair inject {subgroup} {family} R_V={rv}", flush=True)
                curve, response, summary = inject_curve(
                    subgroup,
                    family,
                    rv,
                    classifier,
                    donor_pool,
                    donor_model,
                    normal_points,
                    rng,
                    posterior_ids,
                    posterior_cube,
                    anchor_map,
                    template_magnitudes,
                    template_weights,
                    branch_sigma[rv],
                    age_posterior,
                )
                curves.append(curve)
                responses.append(response)
                summaries.append(summary)

    curve_path = w.PROC / f"wp5_completeness_curves_{output_version}.parquet"
    response_path = w.PROC / f"wp5_injection_response_{output_version}.parquet"
    curve_result = pd.concat(curves, ignore_index=True)
    response_result = pd.concat(responses, ignore_index=True)
    if args.baseline_only:
        if not curve_path.exists() or not response_path.exists():
            raise RuntimeError("--baseline-only requires existing repair response outputs")
        old_curve = pd.read_parquet(curve_path)
        old_response = pd.read_parquet(response_path)
        keep_curve = ~(
            old_curve["family"].eq("PARSEC") & old_curve["R_V"].eq(3.1)
        )
        keep_response = ~(
            old_response["family"].eq("PARSEC") & old_response["R_V"].eq(3.1)
        )
        curve_result = pd.concat(
            [old_curve.loc[keep_curve], curve_result], ignore_index=True
        )
        response_result = pd.concat(
            [old_response.loc[keep_response], response_result], ignore_index=True
        )
    curve_result.to_parquet(curve_path, index=False)
    response_result.to_parquet(response_path, index=False)
    inputs = [
        w.PROC / "wp1_gaia_narrow.parquet",
        w.PROC / "wp1_2mass_join.parquet",
        w.PROC / "wp2_members.parquet",
        w.TABLES / "wp2_subgroup_labels.parquet",
        w.PROC / f"wp3_extinction_{REPAIR_VERSION}.parquet",
        w.PROC / f"wp3_extinction_posterior_{REPAIR_VERSION}.npz",
        w.PROC / f"wp4_mass_posteriors_{REPAIR_VERSION}.parquet",
        w.PROC / f"wp4_age_posteriors_{REPAIR_VERSION}.parquet",
        w.PROVENANCE / "wp2_membership_manifest.json",
        w.PROVENANCE / "wp3_repair_execution.json",
    ]
    record = {
        "repair_version": REPAIR_VERSION,
        "output_version": output_version,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_injections_repair.py",
        "status": "SUCCESS",
        "execution_scope": (
            "baseline_PARSEC_RV3.1_replacement"
            if args.baseline_only
            else "all_branches"
        ),
        "seed": w.SEED,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "sklearn": sklearn.__version__,
        },
        "preserved_configuration": {
            "mass_grid": w.MASS_GRID.tolist(),
            "parent_mass_range_Msun": [
                float(w.MASS_GRID.min()),
                float(w.MASS_GRID.max()),
            ],
            "observed_calibration_window_Msun": [
                w.CALIBRATION_NOMINAL_LO,
                w.CALIBRATION_HI,
            ],
            "injections_per_mass_branch": w.N_INJECT_PER_MASS,
            "membership_qmc_points": w.MEMBERSHIP_QMC_POINTS,
            "membership_threshold": w.MEMBERSHIP_THRESHOLD,
            "binary_fraction": w.F_BINARY,
            "minimum_mass_ratio": w.Q_MIN,
        },
        "repair_change": (
            "synthetic G/BP/RP/J/H/Ks and their cloned errors pass through the "
            "actual repair_v1 WP3 posterior estimator, then the repair_v1 WP4 "
            "posterior sampler (16 stratified draws for response throughput)"
        ),
        "qmc_validation": validation,
        "inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path) for path in inputs
        },
        "branches": summaries,
        "outputs": {
            str(path.relative_to(w.ROOT)): {
                "sha256": w.sha256(path),
                "bytes": path.stat().st_size,
                "rows": int(len(pd.read_parquet(path))),
            }
            for path in [curve_path, response_path]
        },
        "frozen_wp5_outputs_overwritten": False,
    }
    provenance_name = (
        "wp5_injections_repair_execution.json"
        if output_version == REPAIR_VERSION
        else f"wp5_injections_repair_execution_{output_version}.json"
    )
    w.write_json(w.PROVENANCE / provenance_name, record)
    print(
        f"wrote {curve_path.relative_to(w.ROOT)} and "
        f"{response_path.relative_to(w.ROOT)}"
    )


if __name__ == "__main__":
    main()
