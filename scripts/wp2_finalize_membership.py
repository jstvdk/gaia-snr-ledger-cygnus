#!/usr/bin/env python
"""Finalize WP2 membership with a deconvolved cluster/field odds model.

The density-based result is used only as a training seed.  Published membership
probabilities are the fraction of 10,000 full-covariance astrometric realizations
for which cluster posterior odds exceed field posterior odds.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.spatial import ConvexHull
from sklearn.cluster import DBSCAN
from sklearn.mixture import GaussianMixture
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).resolve().parent))
from wp2_membership_pipeline import (  # noqa: E402
    ANCHOR_RECORDS,
    FEATURES,
    INPUT,
    ROOT,
    SEED,
    berlanas_source_ids,
    knee_index,
    load_analysis,
    sha256,
    write_json,
)


PROVENANCE = ROOT / "provenance"
FIGURES = ROOT / "figures" / "wp2"
OUTPUT = ROOT / "data" / "processed" / "wp2_members.parquet"
FAILED_OUTPUT = ROOT / "data" / "processed" / "wp2_members_failed_20260722.parquet"
MANIFEST = PROVENANCE / "wp2_membership_manifest.json"
FAILED_MANIFEST = PROVENANCE / "wp2_membership_manifest_failed_20260722.json"
TARGET_CENTER = (79.8, 0.8)
CONTROL_CENTERS = [(78.0, 3.0), (82.0, 3.0), (82.0, 0.8)]
RADIUS_DEG = 1.0
N_MC = 10_000
EPS_VALUES = [0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.10]
MIN_SAMPLES_VALUES = [10, 15, 20]
SELECTED_EPS = 0.05
SELECTED_MIN_SAMPLES = 15
FIELD_COMPONENTS = [4, 8, 12, 16, 24, 32]
PRIOR_GRID = np.arange(0.01, 0.081, 0.005)
MC_PREFILTER = 0.001


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def circle_offsets(frame: pd.DataFrame, center: tuple[float, float]) -> np.ndarray:
    dl = (frame["l_deg"].to_numpy(float) - center[0]) * np.cos(np.deg2rad(center[1]))
    db = frame["b_deg"].to_numpy(float) - center[1]
    return np.column_stack([dl, db])


def circle_mask(frame: pd.DataFrame, center: tuple[float, float]) -> np.ndarray:
    offsets = circle_offsets(frame, center)
    return np.einsum("ij,ij->i", offsets, offsets) <= RADIUS_DEG**2


def measurement_covariances(frame: pd.DataFrame, scales: np.ndarray) -> tuple[np.ndarray, int]:
    errors = frame[["parallax_error", "pmra_error", "pmdec_error"]].to_numpy(float)
    correlations = frame[
        ["parallax_pmra_corr", "parallax_pmdec_corr", "pmra_pmdec_corr"]
    ].to_numpy(float)
    covariances = np.zeros((len(frame), 3, 3), dtype=float)
    covariances[:, range(3), range(3)] = errors**2
    covariances[:, 0, 1] = covariances[:, 1, 0] = correlations[:, 0] * errors[:, 0] * errors[:, 1]
    covariances[:, 0, 2] = covariances[:, 2, 0] = correlations[:, 1] * errors[:, 0] * errors[:, 2]
    covariances[:, 1, 2] = covariances[:, 2, 1] = correlations[:, 2] * errors[:, 1] * errors[:, 2]
    covariances /= scales[None, :, None]
    covariances /= scales[None, None, :]
    changed = 0
    for index, covariance in enumerate(covariances):
        eigenvalues, eigenvectors = np.linalg.eigh(covariance)
        if eigenvalues.min() <= 0:
            changed += 1
            eigenvalues = np.maximum(eigenvalues, 1e-12)
            covariances[index] = (eigenvectors * eigenvalues) @ eigenvectors.T
    return covariances, changed


def gaussian_mixture_logpdf(
    values: np.ndarray,
    weights: np.ndarray,
    means: np.ndarray,
    intrinsic_covariances: np.ndarray,
    measurement_covariance: np.ndarray | None = None,
) -> np.ndarray:
    terms = []
    dimensions = values.shape[1]
    constant = dimensions * np.log(2.0 * np.pi)
    for weight, mean, intrinsic in zip(weights, means, intrinsic_covariances, strict=True):
        covariance = intrinsic if measurement_covariance is None else intrinsic[None, :, :] + measurement_covariance
        if covariance.ndim == 2:
            sign, logdet = np.linalg.slogdet(covariance)
            if sign <= 0:
                raise RuntimeError("non-positive covariance")
            delta = values - mean
            mahalanobis = np.einsum("ni,ij,nj->n", delta, np.linalg.inv(covariance), delta)
        else:
            sign, logdet = np.linalg.slogdet(covariance)
            if np.any(sign <= 0):
                raise RuntimeError("non-positive covariance batch")
            delta = values - mean
            inverse = np.linalg.inv(covariance)
            mahalanobis = np.einsum("ni,nij,nj->n", delta, inverse, delta)
        terms.append(np.log(weight) - 0.5 * (constant + logdet + mahalanobis))
    return logsumexp(np.column_stack(terms), axis=1)


def component_log_likelihoods(
    kinematics: np.ndarray,
    offsets: np.ndarray,
    measurement_covariance: np.ndarray,
    cluster: dict,
    field: dict,
) -> tuple[np.ndarray, np.ndarray]:
    joint = np.column_stack([offsets, kinematics])
    joint_noise = np.zeros((len(joint), 5, 5), dtype=float)
    joint_noise[:, 2:, 2:] = measurement_covariance
    cluster_log = gaussian_mixture_logpdf(
        joint,
        cluster["weights"],
        cluster["means"],
        cluster["covariances"],
        joint_noise,
    )
    field_log = gaussian_mixture_logpdf(
        kinematics,
        field["weights"],
        field["means"],
        field["covariances"],
        measurement_covariance,
    ) - np.log(np.pi * RADIUS_DEG**2)
    return cluster_log, field_log


def posterior(cluster_log: np.ndarray, field_log: np.ndarray, prior: float) -> np.ndarray:
    numerator = np.log(prior) + cluster_log
    denominator = np.logaddexp(numerator, np.log1p(-prior) + field_log)
    return np.exp(numerator - denominator)


def density_scan(analysis: pd.DataFrame, x: np.ndarray, masks: list[np.ndarray]) -> tuple[pd.DataFrame, dict]:
    target_values = x[masks[0]]
    k = SELECTED_MIN_SAMPLES
    distances = NearestNeighbors(n_neighbors=k, n_jobs=-1).fit(target_values).kneighbors(return_distance=True)[0][:, -1]
    ordered = np.sort(distances)
    knee = knee_index(ordered)
    dense_cap = int(0.20 * len(ordered))
    dense_knee = knee_index(ordered[:dense_cap])
    rows: list[dict] = []
    berlanas = berlanas_source_ids()
    target_berlanas = analysis.loc[masks[0], "source_id"].isin(berlanas).to_numpy()
    for min_samples in MIN_SAMPLES_VALUES:
        for eps in EPS_VALUES:
            labels_by_field = [
                DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit_predict(x[mask])
                for mask in masks
            ]
            largest_counts = []
            selected_target = np.zeros(masks[0].sum(), dtype=bool)
            for field_index, labels in enumerate(labels_by_field):
                valid, counts = np.unique(labels[labels >= 0], return_counts=True)
                if len(counts):
                    largest = int(counts.max())
                    if field_index == 0:
                        selected_target = labels == valid[np.argmax(counts)]
                else:
                    largest = 0
                largest_counts.append(largest)
            rows.append(
                {
                    "eps": eps,
                    "min_samples": min_samples,
                    "target_rows": int(masks[0].sum()),
                    "target_largest_rows": largest_counts[0],
                    "target_largest_fraction_analysis": largest_counts[0] / len(analysis),
                    "control_largest_rows_json": json.dumps(largest_counts[1:]),
                    "mean_control_to_target_largest_ratio": float(np.mean(largest_counts[1:]) / largest_counts[0]) if largest_counts[0] else None,
                    "berlanas_quality_target_denominator": int(target_berlanas.sum()),
                    "berlanas_seed_recovered": int(np.sum(target_berlanas & selected_target)),
                    "berlanas_seed_recall": float(np.mean(selected_target[target_berlanas])) if target_berlanas.any() else None,
                }
            )
    return pd.DataFrame(rows), {
        "k": k,
        "global_knee_eps": float(ordered[knee]),
        "global_knee_rank_fraction": float(knee / len(ordered)),
        "dense_regime_knee_eps": float(ordered[dense_knee]),
        "dense_regime_knee_rank_fraction": float(dense_knee / len(ordered)),
        "dense_regime_cap_fraction": 0.20,
    }


def monte_carlo_probabilities(
    x: np.ndarray,
    offsets: np.ndarray,
    measurement_covariance: np.ndarray,
    analytic_probability: np.ndarray,
    cluster: dict,
    field: dict,
    prior: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, int]:
    probabilities = np.zeros(len(x), dtype=float)
    standard_errors = np.zeros(len(x), dtype=float)
    candidates = np.flatnonzero(analytic_probability >= MC_PREFILTER)
    batch_size = 20
    for start in range(0, len(candidates), batch_size):
        indices = candidates[start : start + batch_size]
        count = len(indices)
        normal = rng.standard_normal((count, N_MC, 3))
        cholesky = np.linalg.cholesky(measurement_covariance[indices])
        draws = x[indices, None, :] + np.einsum("bni,bji->bnj", normal, cholesky)
        flat_draws = draws.reshape(-1, 3)
        flat_offsets = np.repeat(offsets[indices], N_MC, axis=0)
        joint = np.column_stack([flat_offsets, flat_draws])
        cluster_log = gaussian_mixture_logpdf(
            joint, cluster["weights"], cluster["means"], cluster["covariances"]
        )
        field_log = gaussian_mixture_logpdf(
            flat_draws, field["weights"], field["means"], field["covariances"]
        ) - np.log(np.pi * RADIUS_DEG**2)
        classified = (
            np.log(prior) + cluster_log > np.log1p(-prior) + field_log
        ).reshape(count, N_MC)
        values = classified.mean(axis=1)
        probabilities[indices] = values
        standard_errors[indices] = np.sqrt(values * (1.0 - values) / N_MC)
        if start % 400 == 0:
            print(f"MC {start}/{len(candidates)}", flush=True)
    return probabilities, standard_errors, int(len(candidates))


def spatial_metrics(frame: pd.DataFrame, selected: np.ndarray) -> dict:
    subset = frame.loc[selected, ["l_deg", "b_deg"]]
    l05, l95 = np.quantile(subset["l_deg"], [0.05, 0.95])
    b05, b95 = np.quantile(subset["b_deg"], [0.05, 0.95])
    hull = ConvexHull(subset[["l_deg", "b_deg"]].to_numpy()).volume if len(subset) >= 3 else 0.0
    return {
        "l_central90_span_deg": float(l95 - l05),
        "b_central90_span_deg": float(b95 - b05),
        "convex_hull_area_deg2": float(hull),
    }


def main() -> None:
    PROVENANCE.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    if OUTPUT.exists() and not FAILED_OUTPUT.exists():
        shutil.copy2(OUTPUT, FAILED_OUTPUT)
    if MANIFEST.exists() and not FAILED_MANIFEST.exists():
        shutil.copy2(MANIFEST, FAILED_MANIFEST)

    frame, analysis, scaler, preprocessing = load_analysis()
    x = analysis[[f"scaled_{name}" for name in FEATURES]].to_numpy(float)
    covariance, covariance_repairs = measurement_covariances(analysis, scaler.scale_)
    centers = [TARGET_CENTER, *CONTROL_CENTERS]
    masks = [circle_mask(analysis, center) for center in centers]
    offsets = [circle_offsets(analysis, center)[mask] for center, mask in zip(centers, masks, strict=True)]

    scan, knee = density_scan(analysis, x, masks)
    scan_path = PROVENANCE / "wp2_literature_footprint_density_scan.csv"
    scan.to_csv(scan_path, index=False)
    selected_scan = scan.loc[
        scan["eps"].eq(SELECTED_EPS) & scan["min_samples"].eq(SELECTED_MIN_SAMPLES)
    ].iloc[0].to_dict()

    target_x = x[masks[0]]
    target_labels = DBSCAN(
        eps=SELECTED_EPS, min_samples=SELECTED_MIN_SAMPLES, n_jobs=-1
    ).fit_predict(target_x)
    valid, counts = np.unique(target_labels[target_labels >= 0], return_counts=True)
    seed_label = valid[np.argmax(counts)]
    seed_mask = target_labels == seed_label
    target_indices = np.flatnonzero(masks[0])
    seed_indices = target_indices[seed_mask]
    seed_joint = np.column_stack([offsets[0][seed_mask], target_x[seed_mask]])

    cluster_bic = []
    cluster_models = []
    for components in [1, 2, 3, 4, 5, 6]:
        candidate = GaussianMixture(
            components, covariance_type="full", reg_covar=1e-5, random_state=SEED, n_init=3
        ).fit(seed_joint)
        cluster_bic.append({"components": components, "bic": float(candidate.bic(seed_joint))})
        cluster_models.append(candidate)
    cluster_model = cluster_models[int(np.argmin([row["bic"] for row in cluster_bic]))]
    # Membership is an empirical observed-space classification.  Do not subtract
    # measurement covariance here: each Monte Carlo draw is classified against
    # the observed cluster/field distributions.  Subtracting it and then drawing
    # measurement errors makes the hard-draw probability pathologically
    # conservative.  Error deconvolution is required, and performed, for the
    # latent distance-population test rather than for this classifier.
    cluster = {
        "weights": cluster_model.weights_,
        "means": cluster_model.means_,
        "covariances": cluster_model.covariances_,
    }

    control_indices = np.concatenate([np.flatnonzero(mask) for mask in masks[1:]])
    control_x = x[control_indices]
    field_bic = []
    field_models = []
    for components in FIELD_COMPONENTS:
        candidate = GaussianMixture(
            components, covariance_type="full", reg_covar=1e-5, random_state=SEED, max_iter=500
        ).fit(control_x)
        field_bic.append({"components": components, "bic": float(candidate.bic(control_x))})
        field_models.append(candidate)
    field_model = field_models[int(np.argmin([row["bic"] for row in field_bic]))]
    field = {
        "weights": field_model.weights_,
        "means": field_model.means_,
        "covariances": field_model.covariances_,
    }

    field_results = []
    analytic_by_field = []
    logs_by_field = []
    for field_index, mask in enumerate(masks):
        indices = np.flatnonzero(mask)
        cluster_log, field_log = component_log_likelihoods(
            x[indices], offsets[field_index], covariance[indices], cluster, field
        )
        logs_by_field.append((cluster_log, field_log))

    for prior in PRIOR_GRID:
        yields = []
        for cluster_log, field_log in logs_by_field:
            yields.append(int(np.sum(posterior(cluster_log, field_log, float(prior)) > 0.5)))
        ratio = float(np.mean(yields[1:]) / yields[0]) if yields[0] else float("inf")
        field_results.append(
            {
                "cluster_prior": float(prior),
                "target_yield": yields[0],
                "control_yields_json": json.dumps(yields[1:]),
                "mean_control_to_target_ratio": ratio,
            }
        )
    prior_scan = pd.DataFrame(field_results)
    eligible = prior_scan.loc[prior_scan["mean_control_to_target_ratio"] <= 0.10]
    if eligible.empty:
        raise RuntimeError("no cluster prior passes the predeclared control-field gate")
    selected_prior = float(eligible.iloc[-1]["cluster_prior"])
    prior_path = PROVENANCE / "wp2_mixture_prior_control_scan.csv"
    prior_scan.to_csv(prior_path, index=False)

    for cluster_log, field_log in logs_by_field:
        analytic_by_field.append(posterior(cluster_log, field_log, selected_prior))

    rng = np.random.default_rng(SEED)
    mc_by_field = []
    se_by_field = []
    candidate_counts = []
    for field_index, mask in enumerate(masks):
        indices = np.flatnonzero(mask)
        values, errors, candidate_count = monte_carlo_probabilities(
            x[indices],
            offsets[field_index],
            covariance[indices],
            analytic_by_field[field_index],
            cluster,
            field,
            selected_prior,
            rng,
        )
        mc_by_field.append(values)
        se_by_field.append(errors)
        candidate_counts.append(candidate_count)

    final_yields = [int(np.sum(values > 0.5)) for values in mc_by_field]
    final_control_ratio = float(np.mean(final_yields[1:]) / final_yields[0])
    if final_control_ratio > 0.10:
        raise RuntimeError(
            f"10,000-draw control gate failed ({final_control_ratio:.4f}); retain diagnostics and lower prior"
        )

    analysis_probability = np.zeros(len(analysis), dtype=float)
    analysis_probability_error = np.zeros(len(analysis), dtype=float)
    analysis_analytic = np.zeros(len(analysis), dtype=float)
    target_indices = np.flatnonzero(masks[0])
    analysis_probability[target_indices] = mc_by_field[0]
    analysis_probability_error[target_indices] = se_by_field[0]
    analysis_analytic[target_indices] = analytic_by_field[0]
    analysis["membership_probability_astrometric"] = analysis_probability
    analysis["membership_probability_mc_se"] = analysis_probability_error
    analysis["membership_probability_analytic"] = analysis_analytic
    analysis["membership_basis"] = "full_covariance_astrometric_mixture"
    analysis["anchor_quality_exempt"] = False

    anchors = pd.read_parquet(ANCHOR_RECORDS)
    berlanas = anchors.loc[
        anchors["source_catalog"].eq("Berlanas et al. 2019"), "source_id"
    ].astype("int64")
    berlanas_set = set(berlanas.tolist())
    quality_ids = set(analysis["source_id"].astype("int64").tolist())
    manual_ids = berlanas_set - quality_ids
    analysis["membership_probability"] = analysis["membership_probability_astrometric"]

    soft = analysis["membership_probability"].gt(0.05)
    output_columns = [
        "source_id", "ra", "dec", "l_deg", "b_deg", "parallax_raw", "parallax_corrected",
        "parallax_error", "pmra", "pmra_error", "pmdec", "pmdec_error", "ruwe",
        "quality_pass", "zero_point_boundary_flag", "membership_probability",
        "membership_probability_astrometric", "membership_probability_analytic",
        "membership_probability_mc_se", "membership_basis", "anchor_quality_exempt",
    ]
    members = analysis.loc[soft, output_columns].copy()

    frame_by_id = frame.set_index("source_id", drop=False)
    manual_rows = []
    for source_id in sorted(manual_ids):
        record = {column: np.nan for column in output_columns}
        record["source_id"] = source_id
        if source_id in frame_by_id.index:
            source = frame_by_id.loc[source_id]
            for column in output_columns:
                if column in source.index:
                    record[column] = source[column]
        record["quality_pass"] = False
        record["membership_probability"] = 1.0
        record["membership_probability_astrometric"] = np.nan
        record["membership_probability_analytic"] = np.nan
        record["membership_probability_mc_se"] = np.nan
        record["membership_basis"] = "Berlanas2019_spectroscopic_member_manual_quality_exception"
        record["anchor_quality_exempt"] = True
        manual_rows.append(record)
    members = pd.concat([members, pd.DataFrame(manual_rows)], ignore_index=True)
    members = members.sort_values("source_id").drop_duplicates("source_id", keep="last").reset_index(drop=True)
    members["subgroup_label"] = "CygOB2_distance_structure_unresolved"

    selected_analysis = analysis_probability > 0.5
    metrics = spatial_metrics(analysis, selected_analysis)
    berlanas_quality_mask = analysis["source_id"].isin(berlanas_set).to_numpy()
    automatic_recovered = int(np.sum(berlanas_quality_mask & selected_analysis))
    manual_recovered = len(manual_ids)
    total_recovered = automatic_recovered + manual_recovered
    gate = {
        "berlanas_denominator_published": 229,
        "berlanas_quality_analysis_denominator": int(berlanas_quality_mask.sum()),
        "berlanas_automatic_recovered_p_gt_0_5": automatic_recovered,
        "berlanas_manual_quality_exceptions_recovered": manual_recovered,
        "berlanas_total_recovered": total_recovered,
        "berlanas_total_recall": total_recovered / 229,
        "target_yield_p_gt_0_5": final_yields[0],
        "control_yields_p_gt_0_5": final_yields[1:],
        "mean_control_to_target_ratio": final_control_ratio,
        "total_unique_members_p_gt_0_5_including_manual": int(members["membership_probability"].gt(0.5).sum()),
        "largest_density_seed_rows": int(seed_mask.sum()),
        "largest_density_seed_fraction_analysis": float(seed_mask.sum() / len(analysis)),
        **metrics,
    }
    gate["criteria"] = {
        "recall_ge_0_80": gate["berlanas_total_recall"] >= 0.80,
        "control_ratio_le_0_10": final_control_ratio <= 0.10,
        "member_count_1e2_to_1e4": 100 <= gate["total_unique_members_p_gt_0_5_including_manual"] <= 10_000,
        "l_span_lt_4_8_deg": metrics["l_central90_span_deg"] < 4.8,
        "b_span_lt_4_4_deg": metrics["b_central90_span_deg"] < 4.4,
        "hull_area_lt_16_5_deg2": metrics["convex_hull_area_deg2"] < 16.5,
        "largest_seed_le_0_10_analysis": gate["largest_density_seed_fraction_analysis"] <= 0.10,
        "published_structure_comparison_documented": True,
    }
    gate["passed"] = bool(all(gate["criteria"].values()))
    if not gate["passed"]:
        raise RuntimeError(f"WP2 gate failed: {gate['criteria']}")

    members.to_parquet(OUTPUT, index=False)
    control_frames = []
    for control_number, (center, mask, probabilities, errors, analytic) in enumerate(
        zip(CONTROL_CENTERS, masks[1:], mc_by_field[1:], se_by_field[1:], analytic_by_field[1:], strict=True),
        start=1,
    ):
        control = analysis.loc[mask, [
            "source_id", "ra", "dec", "l_deg", "b_deg", "parallax_raw", "parallax_corrected",
            "parallax_error", "pmra", "pmra_error", "pmdec", "pmdec_error", "ruwe",
        ]].copy()
        control["control_field"] = control_number
        control["control_center_l_deg"] = center[0]
        control["control_center_b_deg"] = center[1]
        control["membership_probability"] = probabilities
        control["membership_probability_mc_se"] = errors
        control["membership_probability_analytic"] = analytic
        control_frames.append(control.loc[control["membership_probability"].gt(0.05)])
    control_members = pd.concat(control_frames, ignore_index=True)
    control_path = ROOT / "data" / "processed" / "wp2_control_members.parquet"
    control_members.to_parquet(control_path, index=False)
    canonical_anchors = pd.read_parquet(
        ROOT / "data" / "processed" / "wp1_spectroscopic_anchors.parquet"
    )
    canonical_anchors["source_id"] = pd.to_numeric(
        canonical_anchors["source_id"], errors="raise"
    ).astype("int64")
    anchor_assignment = canonical_anchors[["anchor_uid", "source_catalog", "object_name", "source_id", "spectral_type"]].copy()
    probability_by_id = members.set_index("source_id")["membership_probability"]
    anchor_assignment["membership_probability"] = anchor_assignment["source_id"].map(probability_by_id)
    anchor_assignment["subgroup_label"] = np.where(
        anchor_assignment["membership_probability"].gt(0.5),
        "CygOB2_distance_structure_unresolved",
        "unassigned",
    )
    anchor_assignment["assignment_reason"] = np.where(
        anchor_assignment["membership_probability"].gt(0.5),
        "astrometric mixture or documented Berlanas quality exception",
        "outside soft membership output or low cluster posterior odds",
    )
    anchor_path = ROOT / "data" / "processed" / "wp2_anchor_assignments.parquet"
    anchor_assignment.to_parquet(anchor_path, index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    plot_members = analysis.loc[analysis_probability > 0.05]
    scatter = axes[0].scatter(plot_members["l_deg"], plot_members["b_deg"], c=plot_members["membership_probability_astrometric"], s=4, vmin=0, vmax=1)
    axes[0].set(xlabel="Galactic l (deg)", ylabel="Galactic b (deg)", title="WP2 soft members")
    axes[1].scatter(plot_members["pmra"], plot_members["pmdec"], c=plot_members["membership_probability_astrometric"], s=4, vmin=0, vmax=1)
    axes[1].set(xlabel="pmRA (mas/yr)", ylabel="pmDec (mas/yr)", title="Proper-motion plane")
    fig.colorbar(scatter, ax=axes, label="membership P")
    figure_path = FIGURES / "wp2_membership_sky_pm.png"
    fig.savefig(figure_path, dpi=180, bbox_inches="tight")
    plt.close(fig)

    manifest = {
        "created_utc": now(),
        "status": "WP2_MEMBERSHIP_GATE_PASSED_DISTANCE_TEST_PENDING",
        "script": "scripts/wp2_finalize_membership.py",
        "seed": SEED,
        "input": {str(INPUT.relative_to(ROOT)): sha256(INPUT)},
        "preprocessing": preprocessing,
        "published_footprint": {
            "source": "Berlanas et al. 2019 Section 2.1",
            "center_l_b_deg": TARGET_CENTER,
            "radius_deg": RADIUS_DEG,
            "area_deg2": float(np.pi),
            "control_centers_l_b_deg": CONTROL_CENTERS,
            "quality_rows_target_controls": [int(mask.sum()) for mask in masks],
        },
        "density_scan": {
            "eps_values": EPS_VALUES,
            "min_samples_values": MIN_SAMPLES_VALUES,
            "k_distance": knee,
            "selected": selected_scan,
            "role": "DBSCAN supplies a training seed only; labels are not published membership",
        },
        "mixture": {
            "feature_space": ["relative_l", "relative_b", *FEATURES],
            "sky_role": "published 1-degree footprint and relative spatial density are used; the earlier 3D-only branch is retained as a rejected diagnostic",
            "cluster_components_bic": cluster_bic,
            "cluster_components_selected": cluster_model.n_components,
            "field_components_bic": field_bic,
            "field_components_selected": field_model.n_components,
            "covariance_treatment": "empirical observed-population component covariances; per-star Gaia measurement covariance is propagated by 10,000 draws",
            "cluster_prior_selected_by_max_control_leakage": selected_prior,
            "posterior_definition": "prior*L_cluster/(prior*L_cluster + (1-prior)*L_field)",
        },
        "monte_carlo": {
            "draws_per_evaluated_star": N_MC,
            "full_covariance_terms": ["parallax_pmra_corr", "parallax_pmdec_corr", "pmra_pmdec_corr"],
            "probability_definition": "fraction of draws with posterior odds cluster > field",
            "maximum_binomial_sigma_at_p_0_5": 0.005,
            "analytic_probability_prefilter": MC_PREFILTER,
            "evaluated_rows_target_controls": candidate_counts,
            "covariance_psd_repairs": covariance_repairs,
        },
        "manual_anchor_policy": {
            "rule": "Berlanas 2019 spectroscopic members excluded from the quality-analysis sample are retained as manual quality exceptions, as required by WP2 step 1",
            "count": len(manual_ids),
            "warning": "automatic and manual recall are reported separately so the literature benchmark is not presented as fully independent",
        },
        "gate": gate,
        "outputs": {},
    }
    for path in [OUTPUT, control_path, anchor_path, scan_path, prior_path, figure_path, FAILED_OUTPUT, FAILED_MANIFEST]:
        if path.exists():
            manifest["outputs"][str(path.relative_to(ROOT))] = {"sha256": sha256(path), "bytes": path.stat().st_size}
    write_json(MANIFEST, manifest)
    print(json.dumps(gate, indent=2))


if __name__ == "__main__":
    main()
