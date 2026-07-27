#!/usr/bin/env python3
"""WP5 step 1: end-to-end catalogue-level synthetic-star injections.

This is deliberately a catalogue-level experiment.  Gaia epoch images and the
upstream AGIS source-detection pipeline are not public/re-runnable, so the
experiment cannot be an image-level artificial-star test.  Instead, synthetic
association stars receive real observational states cloned from nearby Gaia
DR3 sources in the one-degree field and are passed through the frozen query,
quality, mixture-membership, P>0.5, and mass-readiness logic.

Outputs
-------
data/processed/wp5_completeness_curves.parquet
provenance/wp5_injections_execution.json
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.special import ndtri
from scipy.stats import qmc
from sklearn.isotonic import IsotonicRegression
from sklearn.neighbors import NearestNeighbors

import wp5_common as w
from wp2_finalize_membership import (
    TARGET_CENTER,
    circle_mask,
    measurement_covariances,
)
from wp2_membership_pipeline import FEATURES
from wp3_common import gaia_mag_error
from wp3_extinction_law import band_coefficients
from wp4_masses import SIG_COL as MASS_SIG_COL
from wp4_masses import SIG_MG as MASS_SIG_MG


def wilson_interval(success: np.ndarray, total: int, z: float = 1.96):
    p = success / total
    denominator = 1.0 + z**2 / total
    center = (p + z**2 / (2.0 * total)) / denominator
    half = z * np.sqrt(p * (1.0 - p) / total + z**2 / (4.0 * total**2)) / denominator
    return center - half, center + half


def sobol_normals(n_points: int) -> np.ndarray:
    if n_points & (n_points - 1):
        raise ValueError("Sobol point count must be a power of two")
    sampler = qmc.Sobol(d=3, scramble=True, seed=w.SEED)
    uniform = sampler.random_base2(int(np.log2(n_points)))
    return ndtri(np.clip(uniform, 1e-9, 1.0 - 1e-9))


def build_donor_pool(classifier: w.WP2Classifier) -> tuple[pd.DataFrame, NearestNeighbors]:
    """Real-field donor pool with query, quality, colour, and band flags."""
    frame = classifier.frame.copy()
    target = circle_mask(frame, TARGET_CENTER)
    frame = frame.loc[target].reset_index(drop=True)
    tmass = pd.read_parquet(
        w.PROC / "wp1_2mass_join.parquet",
        columns=["source_id", "has_2mass_psc_match", "has_complete_jhk"],
    )
    frame = frame.merge(tmass, on="source_id", how="left", validate="one_to_one")
    frame["bp_available"] = (
        frame["phot_bp_mean_mag"].notna()
        & ~frame["phot_bp_mean_mag.mask"].fillna(True)
    )
    frame["rp_available"] = (
        frame["phot_rp_mean_mag"].notna()
        & ~frame["phot_rp_mean_mag.mask"].fillna(True)
    )
    frame["bp_rp"] = frame["phot_bp_mean_mag"] - frame["phot_rp_mean_mag"]
    frame["has_cmd_photometry"] = (
        frame["phot_g_mean_mag"].notna()
        & frame["bp_available"]
        & frame["rp_available"]
    )
    # Match location and magnitude first.  Colour is used when choosing among
    # the 24 neighbours so that red-source astrometric-solution failures remain
    # represented without making missing BP/RP impossible to draw.
    donor_features = np.column_stack(
        [
            (frame["l_deg"].to_numpy(float) - TARGET_CENTER[0]) / 0.08,
            (frame["b_deg"].to_numpy(float) - TARGET_CENTER[1]) / 0.08,
            frame["phot_g_mean_mag"].to_numpy(float) / 0.30,
        ]
    )
    model = NearestNeighbors(n_neighbors=24, algorithm="auto", n_jobs=-1).fit(
        donor_features
    )
    return frame, model


def choose_observational_donors(
    donor_pool: pd.DataFrame,
    donor_model: NearestNeighbors,
    l_deg: np.ndarray,
    b_deg: np.ndarray,
    g_mag: np.ndarray,
    bp_rp: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    features = np.column_stack(
        [
            (l_deg - TARGET_CENTER[0]) / 0.08,
            (b_deg - TARGET_CENTER[1]) / 0.08,
            g_mag / 0.30,
        ]
    )
    _, neighbours = donor_model.kneighbors(features)
    chosen = np.empty(len(features), dtype=int)
    donor_colour = donor_pool["bp_rp"].to_numpy(float)
    for index, candidates in enumerate(neighbours):
        colour_delta = donor_colour[candidates] - bp_rp[index]
        weights = np.exp(-0.5 * (colour_delta / 0.45) ** 2)
        # Missing-colour donors represent BP/RP acquisition failures and must
        # remain injectable.  Give them the median finite colour-match weight.
        finite = np.isfinite(weights)
        replacement = float(np.median(weights[finite])) if finite.any() else 1.0
        weights[~finite] = replacement
        if not np.isfinite(weights).all() or weights.sum() <= 0:
            weights = np.ones(len(candidates))
        chosen[index] = rng.choice(candidates, p=weights / weights.sum())
    return chosen


def load_member_donors(subgroup: str, rv: float) -> pd.DataFrame:
    members = pd.read_parquet(w.PROC / "wp2_members.parquet")
    labels = pd.read_parquet(w.TABLES / "wp2_subgroup_labels.parquet")
    extinction = pd.read_parquet(w.PROC / "wp3_extinction.parquet")
    base = members[
        members["membership_probability"].gt(0.5)
        & ~members["anchor_quality_exempt"]
    ].drop(columns=["subgroup"], errors="ignore")
    base = base.merge(labels[["source_id", "subgroup"]], on="source_id", how="inner")
    base = base.merge(
        extinction[
            [
                "source_id",
                f"av_rv{rv:.1f}",
                f"av_err_rv{rv:.1f}",
            ]
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
        raise RuntimeError(f"too few injection donors for {subgroup}/R_V={rv}")
    return base


def validate_qmc(classifier: w.WP2Classifier, normal_points: np.ndarray) -> dict:
    """Compare the 128-point injection classifier with frozen 10k WP2 labels."""
    members = pd.read_parquet(w.PROC / "wp2_members.parquet")
    test = members[
        members["membership_probability_astrometric"].notna()
        & members["quality_pass"]
    ].copy()
    covariance_columns = [
        "source_id",
        "parallax_pmra_corr",
        "parallax_pmdec_corr",
        "pmra_pmdec_corr",
    ]
    test = test.merge(
        classifier.analysis[covariance_columns],
        on="source_id",
        how="left",
        validate="one_to_one",
    )
    # Keep validation bounded but span the whole probability distribution.
    test = (
        test.assign(
            probability_bin=pd.cut(
                test["membership_probability_astrometric"],
                np.linspace(0.0, 1.0, 21),
                include_lowest=True,
            )
        )
        .groupby("probability_bin", observed=True, group_keys=False)
        .apply(lambda group: group.sample(min(len(group), 40), random_state=w.SEED))
        .reset_index(drop=True)
    )
    x = classifier.scaler.transform(test[FEATURES])
    covariance, repairs = measurement_covariances(test, classifier.scaler.scale_)
    offsets = np.column_stack(
        [
            (test["l_deg"].to_numpy(float) - TARGET_CENTER[0])
            * np.cos(np.deg2rad(TARGET_CENTER[1])),
            test["b_deg"].to_numpy(float) - TARGET_CENTER[1],
        ]
    )
    approximate = w.qmc_membership_probabilities(
        x, offsets, covariance, classifier, normal_points
    )
    exact = test["membership_probability_astrometric"].to_numpy(float)
    exact_label = exact > w.MEMBERSHIP_THRESHOLD
    approximate_label = approximate > w.MEMBERSHIP_THRESHOLD
    return {
        "rows": int(len(test)),
        "covariance_repairs": int(repairs),
        "decision_agreement": float(np.mean(exact_label == approximate_label)),
        "false_negative_count": int(np.sum(exact_label & ~approximate_label)),
        "false_positive_count": int(np.sum(~exact_label & approximate_label)),
        "median_abs_probability_difference": float(np.median(np.abs(exact - approximate))),
        "p95_abs_probability_difference": float(np.quantile(np.abs(exact - approximate), 0.95)),
    }


def inject_curve(
    subgroup: str,
    family: str,
    rv: float,
    classifier: w.WP2Classifier,
    donor_pool: pd.DataFrame,
    donor_model: NearestNeighbors,
    normal_points: np.ndarray,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    age = w.fitted_age(subgroup, family, rv)
    iso, native_age = w.load_isochrone_at_age(family, age)
    member_donors = load_member_donors(subgroup, rv)
    p_weights = member_donors["membership_probability"].to_numpy(float, copy=True)
    p_weights /= p_weights.sum()

    primary_mass = np.repeat(w.MASS_GRID, w.N_INJECT_PER_MASS)
    n_injected = len(primary_mass)
    donor_index = rng.choice(
        len(member_donors), size=n_injected, replace=True, p=p_weights
    )
    member = member_donors.iloc[donor_index].reset_index(drop=True)

    binary = rng.random(n_injected) < w.F_BINARY
    q = rng.uniform(w.Q_MIN, 1.0, n_injected)
    secondary_mass = np.where(binary, q * primary_mass, 0.0)
    absolute = w.interpolate_photometry(iso, primary_mass, secondary_mass)

    av = member[f"av_rv{rv:.1f}"].to_numpy(float)
    av_error = member[f"av_err_rv{rv:.1f}"].to_numpy(float)
    av_true = np.maximum(0.0, av + rng.normal(0.0, av_error))
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

    # Clone the real donor's magnitude-error state, then perturb the synthetic
    # central photometry.  Donors are already matched in G, colour, and sky.
    photometric_error = {}
    observed_mag = {}
    for band, donor_mag, donor_flux_error in [
        ("G", "phot_g_mean_mag", "phot_g_mean_flux_error"),
        ("BP", "phot_bp_mean_mag", "phot_bp_mean_flux_error"),
        ("RP", "phot_rp_mean_mag", "phot_rp_mean_flux_error"),
    ]:
        photometric_error[band] = gaia_mag_error(
            observed[donor_mag].to_numpy(float),
            observed[donor_flux_error].to_numpy(float),
            band,
        )
        observed_mag[band] = apparent[band] + rng.normal(
            0.0, photometric_error[band]
        )

    # One noisy central astrometric measurement, followed by the WP2
    # full-covariance posterior-odds integration around that measurement.
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

    # Frozen query: observed G and raw parallax.  The donor's Lindegren zero
    # point is added back before applying the raw-query parallax bounds.
    raw_parallax = (
        central_measured[:, 0] + observed["parallax_zero_point"].to_numpy(float)
    )
    query_pass = (
        (observed_mag["G"] < 19.0)
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

    # Run the synthetic photometry through the same WP4 inverse mass mapping.
    # This response is essential: the nearest-isochrone estimator and binaries
    # can move true masses between observed bins, so a scalar completeness curve
    # alone is not a sufficient forward model for the mass-function residuals.
    av_estimate = av_true + rng.normal(0.0, av_error)
    recovered_colour = (
        observed_mag["BP"]
        - observed_mag["RP"]
        - (coefficients["BP"] - coefficients["RP"]) * av_estimate
    )
    recovered_mg = (
        observed_mag["G"]
        - w.DIST_MODULUS
        - coefficients["G"] * av_estimate
    )
    iso_mass = iso["Mass"].to_numpy(float)
    iso_colour = (iso["BP0"] - iso["RP0"]).to_numpy(float)
    iso_mg = iso["G0"].to_numpy(float)
    recovered_mass = np.full(n_injected, np.nan)
    for start in range(0, n_injected, 512):
        stop = min(start + 512, n_injected)
        distance = (
            (
                recovered_colour[start:stop, None] - iso_colour[None, :]
            )
            / MASS_SIG_COL
        ) ** 2 + (
            (recovered_mg[start:stop, None] - iso_mg[None, :]) / MASS_SIG_MG
        ) ** 2
        nearest = np.argmin(distance, axis=1)
        recovered_mass[start:stop] = iso_mass[nearest]
    recovered_mass[~recovered] = np.nan

    rows = []
    for mass in w.MASS_GRID:
        select = primary_mass == mass
        n = int(select.sum())
        stages = {
            "n_query": int(query_pass[select].sum()),
            "n_quality": int(quality_pass[select].sum()),
            "n_cmd_ready": int(cmd_ready[select].sum()),
            "n_membership": int(membership_pass[select].sum()),
            "n_recovered": int(recovered[select].sum()),
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
            "G_observed": observed_mag["G"],
            "BP_RP_observed": observed_mag["BP"] - observed_mag["RP"],
            "query_pass": query_pass,
            "quality_pass": quality_pass,
            "cmd_ready": cmd_ready,
            "membership_probability_qmc": membership_probability,
            "membership_pass": membership_pass,
            "recovered_mass": recovered_mass,
        }
    )
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
        "av_draw_median": float(np.median(av_true)),
        "g_mag_range": [float(np.min(apparent["G"])), float(np.max(apparent["G"]))],
        "recovered": int(recovered.sum()),
    }
    return curve, response, summary


def main() -> None:
    w.PROC.mkdir(parents=True, exist_ok=True)
    w.PROVENANCE.mkdir(parents=True, exist_ok=True)
    classifier = w.reconstruct_wp2_classifier()
    donor_pool, donor_model = build_donor_pool(classifier)
    normal_points = sobol_normals(w.MEMBERSHIP_QMC_POINTS)
    validation = validate_qmc(classifier, normal_points)
    if validation["decision_agreement"] < 0.97:
        raise RuntimeError(
            "128-point QMC membership approximation failed validation: "
            f"{validation['decision_agreement']:.3f}"
        )

    rng = np.random.default_rng(w.SEED)
    curves = []
    responses = []
    summaries = []
    for family in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            for subgroup in w.SUBGROUPS:
                print(f"inject {subgroup} {family} R_V={rv}", flush=True)
                curve, response, summary = inject_curve(
                    subgroup,
                    family,
                    rv,
                    classifier,
                    donor_pool,
                    donor_model,
                    normal_points,
                    rng,
                )
                curves.append(curve)
                responses.append(response)
                summaries.append(summary)

    output = w.PROC / "wp5_completeness_curves.parquet"
    result = pd.concat(curves, ignore_index=True)
    result.to_parquet(output, index=False)
    response_output = w.PROC / "wp5_injection_response.parquet"
    response_result = pd.concat(responses, ignore_index=True)
    response_result.to_parquet(response_output, index=False)
    inputs = [
        w.PROC / "wp1_gaia_narrow.parquet",
        w.PROC / "wp1_2mass_join.parquet",
        w.PROC / "wp2_members.parquet",
        w.TABLES / "wp2_subgroup_labels.parquet",
        w.PROC / "wp3_extinction.parquet",
        w.PROC / "wp3_isochrones_parsec.parquet",
        w.PROC / "wp3_isochrones_mist.parquet",
        w.PROC / "wp4_age_posteriors.parquet",
        w.PROVENANCE / "wp2_membership_manifest.json",
    ]
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_injections.py",
        "status": "SUCCESS",
        "seed": w.SEED,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "sklearn": sklearn.__version__,
        },
        "inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path) for path in inputs
        },
        "experiment": {
            "kind": "catalogue-level synthetic-star injection into real Gaia DR3 field observational states",
            "not_image_level": True,
            "mass_grid_Msun": w.MASS_GRID.tolist(),
            "n_per_mass_subgroup_family_rv": w.N_INJECT_PER_MASS,
            "total_injected": int(result["n_injected"].sum()),
            "binary_fraction": w.F_BINARY,
            "mass_ratio_distribution": f"uniform[{w.Q_MIN},1]",
            "extinction": "spatial member resampling plus Gaussian A_V error, separately per R_V branch",
            "query": "G<19; 0.35<raw parallax<1.10 mas",
            "quality": "exact WP2 quality_pass cloned from nearby real-field observational donor",
            "membership": f"exact WP2 posterior-odds rule with {w.MEMBERSHIP_QMC_POINTS}-point Sobol-normal approximation; P>{w.MEMBERSHIP_THRESHOLD}",
            "mass_readiness": "G, BP, RP all available after observational-state cloning",
            "recovery": "query AND quality AND CMD-ready AND membership",
            "mass_response": (
                "recovered G/BP/RP and A_V are passed through the exact WP4 "
                "nearest-isochrone mass estimator; true-to-recovered migration "
                "is stored for the forward Poisson likelihood"
            ),
            "isotonic_role": "monotone summary used only to locate the completeness lower edge; raw fractions retained",
        },
        "wp2_classifier_reconstruction": classifier.reconstruction,
        "qmc_validation_against_frozen_wp2_10k": validation,
        "donor_pool": {
            "rows": int(len(donor_pool)),
            "field_center_l_b_deg": TARGET_CENTER,
            "radius_deg": 1.0,
            "nearest_neighbours": 24,
            "matching_scales": {"l_deg": 0.08, "b_deg": 0.08, "G_mag": 0.30, "BP_RP_mag": 0.45},
        },
        "branch_summaries": summaries,
        "output": {
            str(output.relative_to(w.ROOT)): {
                "sha256": w.sha256(output),
                "bytes": output.stat().st_size,
                "rows": int(len(result)),
            },
            str(response_output.relative_to(w.ROOT)): {
                "sha256": w.sha256(response_output),
                "bytes": response_output.stat().st_size,
                "rows": int(len(response_result)),
            },
        },
        "limitations": [
            "Gaia epoch images and AGIS are not re-runnable; catalogue acquisition is approximated by real-field observational-state cloning.",
            "The WP3 extinction distribution is observed-member conditioned and cannot reveal a completely Gaia-invisible high-A_V population.",
            "The P>0.5 handoff is corrected as part of recovery; upstream soft members lack frozen WP3/WP4 masses and are not retrofitted.",
        ],
    }
    w.write_json(w.PROVENANCE / "wp5_injections_execution.json", record)
    print(json.dumps({"validation": validation, "output_rows": len(result)}, indent=2))


if __name__ == "__main__":
    main()
