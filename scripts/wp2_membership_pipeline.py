#!/usr/bin/env python
"""WP2 membership pipeline: density diagnostics and cluster/field probabilities."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from scipy.spatial import ConvexHull, QhullError
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import RobustScaler
from zero_point import zpt
from hdbscan import HDBSCAN


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data" / "processed" / "wp1_gaia_narrow.parquet"
ANCHOR_RECORDS = ROOT / "data" / "processed" / "wp1_spectroscopic_anchor_records.parquet"
PROVENANCE = ROOT / "provenance"
FIGURES = ROOT / "figures" / "wp2"
SEED = 20260722
FEATURES = ["parallax_corrected", "pmra", "pmdec"]
CORRELATIONS = ["parallax_pmra_corr", "parallax_pmdec_corr", "pmra_pmdec_corr"]
BOX = {"l_min": 77.0, "l_max": 83.0, "b_min": -1.5, "b_max": 4.0}
WRIGHT_FIELD_CENTER = (80.22449980431047, 0.7840736986995647)
FIELD_RADIUS_DEG = float(np.sqrt(1.0 / np.pi))
CONTROL_CENTERS = [(77.8, 0.7840736986995647), (82.2, 0.7840736986995647), (80.22449980431047, -0.7)]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def angular_circle_mask(frame: pd.DataFrame, center: tuple[float, float]) -> np.ndarray:
    l0, b0 = center
    dl = (frame["l_deg"].to_numpy() - l0) * np.cos(np.deg2rad(b0))
    db = frame["b_deg"].to_numpy() - b0
    return dl * dl + db * db <= FIELD_RADIUS_DEG**2


def load_analysis() -> tuple[pd.DataFrame, pd.DataFrame, RobustScaler, dict]:
    frame = pd.read_parquet(INPUT)
    zpt.load_tables()
    coords = SkyCoord(
        ra=frame["ra"].to_numpy() * u.deg,
        dec=frame["dec"].to_numpy() * u.deg,
    )
    ecl_lat = coords.barycentrictrueecliptic.lat.deg
    is_5p = frame["astrometric_params_solved"].eq(31).to_numpy()
    is_6p = frame["astrometric_params_solved"].eq(95).to_numpy()
    g_out = ~frame["phot_g_mean_mag"].between(6, 21).to_numpy()
    nu_out = is_5p & ~frame["nu_eff_used_in_astrometry"].between(1.1, 1.9).fillna(False).to_numpy()
    pc_out = is_6p & ~frame["pseudocolour"].between(1.24, 1.72).fillna(False).to_numpy()
    frame["zero_point_boundary_flag"] = g_out | nu_out | pc_out
    g = frame["phot_g_mean_mag"].clip(6 + 1e-6, 21 - 1e-6).to_numpy()
    nu = frame["nu_eff_used_in_astrometry"].fillna(1.5).clip(1.1 + 1e-6, 1.9 - 1e-6).to_numpy()
    pc = frame["pseudocolour"].fillna(1.5).clip(1.24 + 1e-6, 1.72 - 1e-6).to_numpy()
    frame["parallax_zero_point"] = zpt.get_zpt(
        g,
        nu,
        pc,
        ecl_lat,
        frame["astrometric_params_solved"].to_numpy(),
        _warnings=False,
    )
    frame["parallax_raw"] = frame["parallax"]
    frame["parallax_corrected"] = frame["parallax_raw"] - frame["parallax_zero_point"]
    positive_errors = frame[["parallax_error", "pmra_error", "pmdec_error"]].gt(0).all(axis=1)
    finite_astrometry = frame[FEATURES + CORRELATIONS].notna().all(axis=1)
    frame["quality_pass"] = (
        (frame["ruwe"].fillna(np.inf) < 1.4)
        & (frame["visibility_periods_used"].fillna(0) >= 8)
        & frame["phot_bp_rp_excess_factor"].notna()
        & positive_errors
        & finite_astrometry
        & ~frame["zero_point_boundary_flag"]
    )
    analysis = frame.loc[frame["quality_pass"]].copy().reset_index(drop=True)
    scaler = RobustScaler().fit(analysis[FEATURES])
    scaled = scaler.transform(analysis[FEATURES])
    analysis[[f"scaled_{name}" for name in FEATURES]] = scaled
    summary = {
        "input_rows": int(len(frame)),
        "analysis_rows": int(len(analysis)),
        "quality_fraction": float(len(analysis) / len(frame)),
        "quality_filter": "RUWE<1.4; visibility_periods_used>=8; BP/RP excess present; positive errors; finite 3D covariance; reliable zero point",
        "zero_point_mas": {
            "median": float(frame["parallax_zero_point"].median()),
            "min": float(frame["parallax_zero_point"].min()),
            "max": float(frame["parallax_zero_point"].max()),
        },
        "features": FEATURES,
        "sky_used_in_clustering": False,
        "robust_scaler_center": dict(zip(FEATURES, scaler.center_.tolist(), strict=True)),
        "robust_scaler_scale": dict(zip(FEATURES, scaler.scale_.tolist(), strict=True)),
    }
    return frame, analysis, scaler, summary


def berlanas_source_ids() -> set[int]:
    records = pd.read_parquet(ANCHOR_RECORDS)
    values = records.loc[
        records["source_catalog"].eq("Berlanas et al. 2019"), "source_id"
    ].dropna()
    return set(values.astype("int64"))


def knee_index(sorted_distances: np.ndarray) -> int:
    # Maximum vertical distance below the endpoint chord for a convex increasing curve.
    n = len(sorted_distances)
    lo, hi = int(0.01 * n), int(0.999 * n)
    x = np.linspace(0.0, 1.0, hi - lo)
    y = sorted_distances[lo:hi].astype(float)
    y = (y - y[0]) / max(y[-1] - y[0], np.finfo(float).eps)
    return lo + int(np.argmax(x - y))


def region_yields(frame: pd.DataFrame, selected: np.ndarray) -> tuple[int, list[int], float]:
    target = int(np.sum(selected & angular_circle_mask(frame, WRIGHT_FIELD_CENTER)))
    controls = [int(np.sum(selected & angular_circle_mask(frame, center))) for center in CONTROL_CENTERS]
    ratio = float(np.mean(controls) / target) if target else float("inf")
    return target, controls, ratio


def spatial_metrics(frame: pd.DataFrame, selected: np.ndarray) -> dict:
    subset = frame.loc[selected, ["l_deg", "b_deg"]]
    if len(subset) < 3:
        return {
            "l_central90_span_deg": None,
            "b_central90_span_deg": None,
            "convex_hull_area_deg2": None,
            "touches_opposite_box_sides": False,
        }
    l05, l95 = np.quantile(subset["l_deg"], [0.05, 0.95])
    b05, b95 = np.quantile(subset["b_deg"], [0.05, 0.95])
    try:
        hull_area = float(ConvexHull(subset[["l_deg", "b_deg"]].to_numpy()).volume)
    except QhullError:
        hull_area = 0.0
    margin = 0.05
    touches_l = subset["l_deg"].min() <= BOX["l_min"] + margin and subset["l_deg"].max() >= BOX["l_max"] - margin
    touches_b = subset["b_deg"].min() <= BOX["b_min"] + margin and subset["b_deg"].max() >= BOX["b_max"] - margin
    return {
        "l_central90_span_deg": float(l95 - l05),
        "b_central90_span_deg": float(b95 - b05),
        "convex_hull_area_deg2": hull_area,
        "touches_opposite_box_sides": bool(touches_l or touches_b),
    }


def choose_association_cluster(frame: pd.DataFrame, labels: np.ndarray) -> tuple[int | None, dict]:
    best_label = None
    best_score = -np.inf
    best = {}
    for label in np.unique(labels[labels >= 0]):
        selected = labels == label
        target, controls, ratio = region_yields(frame, selected)
        if target < 20:
            continue
        score = target / (np.mean(controls) + 1.0)
        if score > best_score:
            best_score = score
            best_label = int(label)
            best = {
                "association_target_yield": target,
                "association_control_yields": controls,
                "association_control_ratio": ratio,
                "association_spatial_overdensity_score": float(score),
            }
    return best_label, best


def scan() -> None:
    PROVENANCE.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    _, analysis, _, preprocessing = load_analysis()
    x = analysis[[f"scaled_{name}" for name in FEATURES]].to_numpy(float)
    berlanas = berlanas_source_ids()
    min_samples = 15
    neighbours = NearestNeighbors(n_neighbors=min_samples, n_jobs=-1).fit(x)
    distances = neighbours.kneighbors(x, return_distance=True)[0][:, -1]
    ordered = np.sort(distances)
    knee = knee_index(ordered)
    knee_eps = float(ordered[knee])
    dense_cap = int(0.20 * len(ordered))
    dense_knee = knee_index(ordered[:dense_cap])
    dense_knee_eps = float(ordered[dense_knee])
    quantiles = {
        str(q): float(np.quantile(ordered, q))
        for q in [0.5, 0.75, 0.9, 0.95, 0.975, 0.99, 0.995, 0.999]
    }
    kdist = pd.DataFrame({"rank": np.arange(len(ordered)), "k_distance": ordered})
    kdist_path = PROVENANCE / "wp2_kdistance_min15.parquet"
    kdist.to_parquet(kdist_path, index=False)
    fig, axis = plt.subplots(figsize=(8, 5))
    axis.plot(np.arange(len(ordered)), ordered, lw=1)
    axis.axhline(knee_eps, color="tab:red", ls="--", label=f"knee eps={knee_eps:.4f}")
    axis.set(xlabel="sorted source rank", ylabel="15th-neighbour distance", title="WP2 k-distance diagnostic")
    axis.legend()
    fig.tight_layout()
    figure_path = FIGURES / "kdistance_min15.png"
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)

    base_eps = [0.02, 0.03, 0.04, 0.05, 0.07, 0.10]
    knee_eps_values = [knee_eps * value for value in [0.5, 0.75, 1.0, 1.5]]
    eps_values = sorted({round(value, 5) for value in [*base_eps, *knee_eps_values] if 0.01 <= value <= 0.20})
    rows = []
    for current_min_samples in [10, 15, 20]:
        for eps in eps_values:
            labels = DBSCAN(eps=eps, min_samples=current_min_samples, n_jobs=-1).fit_predict(x)
            valid_labels, counts = np.unique(labels[labels >= 0], return_counts=True)
            largest = int(counts.max()) if len(counts) else 0
            association_label, association = choose_association_cluster(analysis, labels)
            association_selected = labels == association_label if association_label is not None else np.zeros(len(labels), dtype=bool)
            berlanas_in_analysis = analysis["source_id"].isin(berlanas).to_numpy()
            berlanas_denominator = int(berlanas_in_analysis.sum())
            berlanas_recovered = int(np.sum(berlanas_in_analysis & association_selected))
            row = {
                "algorithm": "DBSCAN",
                "eps": eps,
                "min_samples": current_min_samples,
                "clusters": int(len(valid_labels)),
                "assigned_rows": int(np.sum(labels >= 0)),
                "largest_cluster_rows": largest,
                "largest_cluster_fraction_analysis": float(largest / len(analysis)),
                "association_cluster_label": association_label,
                "association_cluster_rows": int(association_selected.sum()),
                "berlanas_in_analysis": berlanas_denominator,
                "berlanas_recovered_in_seed": berlanas_recovered,
                "berlanas_seed_recall": float(berlanas_recovered / berlanas_denominator) if berlanas_denominator else None,
                **association,
                **spatial_metrics(analysis, association_selected),
            }
            rows.append(row)
            print(json.dumps(row), flush=True)
    scan_frame = pd.DataFrame(rows)
    scan_path = PROVENANCE / "wp2_dbscan_scan.csv"
    scan_frame.to_csv(scan_path, index=False)
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp2_membership_pipeline.py --stage scan",
        "seed": SEED,
        "input": {str(INPUT.relative_to(ROOT)): sha256(INPUT)},
        "anchor_records": {str(ANCHOR_RECORDS.relative_to(ROOT)): sha256(ANCHOR_RECORDS)},
        "preprocessing": preprocessing,
        "k_distance": {
            "k": min_samples,
            "knee_index": knee,
            "knee_rank_fraction": float(knee / len(ordered)),
            "knee_eps": knee_eps,
            "quantiles": quantiles,
            "detection": "maximum distance below normalized endpoint chord over ranks 1%-99.9%",
        },
        "scan": {
            "eps_values": eps_values,
            "min_samples_values": [10, 15, 20],
            "feature_space": FEATURES,
            "sky_position_role": "excluded from clustering; used only to identify and validate the association overdensity",
            "association_selection": "density cluster maximizing target/control spatial overdensity with >=20 target-field stars",
        },
        "control_fields": {
            "target_center_l_b_deg": WRIGHT_FIELD_CENTER,
            "control_centers_l_b_deg": CONTROL_CENTERS,
            "radius_deg": FIELD_RADIUS_DEG,
            "area_deg2_each": 1.0,
        },
        "outputs": {
            str(kdist_path.relative_to(ROOT)): sha256(kdist_path),
            str(scan_path.relative_to(ROOT)): sha256(scan_path),
            str(figure_path.relative_to(ROOT)): sha256(figure_path),
        },
    }
    write_json(PROVENANCE / "wp2_density_scan_execution.json", record)
    print(json.dumps(record["k_distance"], indent=2))


def hdbscan_scan(cluster_selection_method: str = "eom") -> None:
    PROVENANCE.mkdir(parents=True, exist_ok=True)
    _, analysis, _, preprocessing = load_analysis()
    x = analysis[[f"scaled_{name}" for name in FEATURES]].to_numpy(float)
    berlanas = berlanas_source_ids()
    berlanas_in_analysis = analysis["source_id"].isin(berlanas).to_numpy()
    rows = []
    for min_samples in [10, 15, 20]:
        for min_cluster_size in [30, 60, 120, 240]:
            model = HDBSCAN(
                min_samples=min_samples,
                min_cluster_size=min_cluster_size,
                metric="euclidean",
                cluster_selection_method=cluster_selection_method,
                approx_min_span_tree=True,
                core_dist_n_jobs=-1,
                prediction_data=True,
            )
            labels = model.fit_predict(x)
            valid_labels, counts = np.unique(labels[labels >= 0], return_counts=True)
            largest = int(counts.max()) if len(counts) else 0
            association_labels = []
            component_records = []
            for label in valid_labels:
                selected = labels == label
                target, controls, ratio = region_yields(analysis, selected)
                target_fraction = float(target / selected.sum())
                # Fixed external sky criterion: at least five times the 1/33 box-area expectation.
                if target >= 10 and target_fraction >= 5.0 / 33.0 and ratio <= 0.5:
                    association_labels.append(int(label))
                    component_records.append(
                        {
                            "label": int(label),
                            "rows": int(selected.sum()),
                            "target_rows": target,
                            "control_rows": controls,
                            "control_ratio": ratio,
                            "target_fraction": target_fraction,
                            "median_membership_strength": float(np.median(model.probabilities_[selected])),
                        }
                    )
            association_selected = np.isin(labels, association_labels)
            target, controls, ratio = region_yields(analysis, association_selected)
            recovered = int(np.sum(berlanas_in_analysis & association_selected))
            row = {
                "algorithm": "HDBSCAN",
                "min_samples": min_samples,
                "min_cluster_size": min_cluster_size,
                "clusters": int(len(valid_labels)),
                "noise_rows": int(np.sum(labels < 0)),
                "largest_cluster_rows": largest,
                "largest_cluster_fraction_analysis": float(largest / len(analysis)),
                "association_component_labels_json": json.dumps(association_labels),
                "association_component_records_json": json.dumps(component_records),
                "association_seed_rows": int(association_selected.sum()),
                "association_target_yield": target,
                "association_control_yields_json": json.dumps(controls),
                "association_control_ratio": ratio,
                "berlanas_in_analysis": int(berlanas_in_analysis.sum()),
                "berlanas_recovered_in_seed": recovered,
                "berlanas_seed_recall": float(recovered / berlanas_in_analysis.sum()),
                **spatial_metrics(analysis, association_selected),
            }
            rows.append(row)
            print(json.dumps(row), flush=True)
    result = pd.DataFrame(rows)
    output = PROVENANCE / f"wp2_hdbscan_{cluster_selection_method}_scan.csv"
    result.to_csv(output, index=False)
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": f"scripts/wp2_membership_pipeline.py --stage hdbscan-{cluster_selection_method}-scan",
        "seed": SEED,
        "input": {str(INPUT.relative_to(ROOT)): sha256(INPUT)},
        "preprocessing": preprocessing,
        "feature_space": FEATURES,
        "parameter_grid": {"min_samples": [10, 15, 20], "min_cluster_size": [30, 60, 120, 240]},
        "association_component_rule": {
            "target_field_rows_min": 10,
            "target_fraction_min": 5.0 / 33.0,
            "control_to_target_ratio_max": 0.5,
            "note": "fixed sky criterion selects which kinematic components correspond to the association; Berlanas labels are not used",
        },
        "output": {
            "file": str(output.relative_to(ROOT)),
            "sha256": sha256(output),
            "rows": int(len(result)),
        },
    }
    write_json(PROVENANCE / f"wp2_hdbscan_{cluster_selection_method}_scan_execution.json", record)


def target_field_scan() -> None:
    PROVENANCE.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    _, analysis, _, preprocessing = load_analysis()
    target_mask = angular_circle_mask(analysis, WRIGHT_FIELD_CENTER)
    control_masks = [angular_circle_mask(analysis, center) for center in CONTROL_CENTERS]
    target = analysis.loc[target_mask].copy().reset_index(drop=True)
    x_target = target[[f"scaled_{name}" for name in FEATURES]].to_numpy(float)
    x_controls = [
        analysis.loc[mask, [f"scaled_{name}" for name in FEATURES]].to_numpy(float)
        for mask in control_masks
    ]
    min_samples_default = 15
    ordered = np.sort(
        NearestNeighbors(n_neighbors=min_samples_default, n_jobs=-1)
        .fit(x_target)
        .kneighbors(x_target, return_distance=True)[0][:, -1]
    )
    knee = knee_index(ordered)
    knee_eps = float(ordered[knee])
    dense_cap = int(0.20 * len(ordered))
    dense_knee = knee_index(ordered[:dense_cap])
    dense_knee_eps = float(ordered[dense_knee])
    kdist_path = PROVENANCE / "wp2_target_kdistance_min15.parquet"
    pd.DataFrame({"rank": np.arange(len(ordered)), "k_distance": ordered}).to_parquet(kdist_path, index=False)
    fig, axis = plt.subplots(figsize=(8, 5))
    axis.plot(np.arange(len(ordered)), ordered, lw=1)
    axis.axhline(knee_eps, color="tab:red", ls="--", label=f"knee eps={knee_eps:.4f}")
    axis.set(xlabel="sorted target-field source rank", ylabel="15th-neighbour distance", title="WP2 target-field k-distance diagnostic")
    axis.legend()
    fig.tight_layout()
    figure_path = FIGURES / "target_kdistance_min15.png"
    fig.savefig(figure_path, dpi=180)
    plt.close(fig)

    eps_values = sorted(
        {
            round(value, 5)
            for value in [
                0.02, 0.03, 0.04, 0.05, 0.07, 0.10,
                knee_eps * 0.5, knee_eps * 0.75, knee_eps, knee_eps * 1.25, knee_eps * 1.5,
            ]
            if 0.01 <= value <= 0.25
        }
    )
    berlanas = berlanas_source_ids()
    berlanas_target = target["source_id"].isin(berlanas).to_numpy()
    rows = []
    for min_samples in [10, 15, 20]:
        for eps in eps_values:
            labels = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit_predict(x_target)
            valid, counts = np.unique(labels[labels >= 0], return_counts=True)
            if len(counts):
                largest_label = int(valid[np.argmax(counts)])
                selected = labels == largest_label
            else:
                largest_label = None
                selected = np.zeros(len(labels), dtype=bool)
            control_assigned = []
            control_largest = []
            for x_control in x_controls:
                control_labels = DBSCAN(eps=eps, min_samples=min_samples, n_jobs=-1).fit_predict(x_control)
                _, control_counts = np.unique(control_labels[control_labels >= 0], return_counts=True)
                control_assigned.append(int(np.sum(control_labels >= 0)))
                control_largest.append(int(control_counts.max()) if len(control_counts) else 0)
            row = {
                "algorithm": "DBSCAN_target_field_seed",
                "eps": eps,
                "min_samples": min_samples,
                "target_rows": int(len(target)),
                "target_clusters": int(len(valid)),
                "target_assigned_rows": int(np.sum(labels >= 0)),
                "target_largest_label": largest_label,
                "target_largest_rows": int(selected.sum()),
                "target_largest_fraction": float(selected.mean()),
                "control_rows_json": json.dumps([int(len(values)) for values in x_controls]),
                "control_assigned_rows_json": json.dumps(control_assigned),
                "control_largest_rows_json": json.dumps(control_largest),
                "control_largest_to_target_largest_ratio": float(np.mean(control_largest) / selected.sum()) if selected.sum() else float("inf"),
                "berlanas_target_denominator": int(berlanas_target.sum()),
                "berlanas_target_recovered": int(np.sum(berlanas_target & selected)),
                "berlanas_target_recall": float(np.sum(berlanas_target & selected) / berlanas_target.sum()),
            }
            rows.append(row)
            print(json.dumps(row), flush=True)
    result = pd.DataFrame(rows)
    output = PROVENANCE / "wp2_target_dbscan_scan.csv"
    result.to_csv(output, index=False)
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp2_membership_pipeline.py --stage target-scan",
        "input": {str(INPUT.relative_to(ROOT)): sha256(INPUT)},
        "preprocessing": preprocessing,
        "target_field": {"center_l_b_deg": WRIGHT_FIELD_CENTER, "radius_deg": FIELD_RADIUS_DEG, "area_deg2": 1.0, "rows": int(len(target))},
        "controls": {"centers_l_b_deg": CONTROL_CENTERS, "radius_deg": FIELD_RADIUS_DEG, "rows": [int(mask.sum()) for mask in control_masks]},
        "k_distance": {
            "k": min_samples_default,
            "knee_eps": knee_eps,
            "knee_rank_fraction": float(knee / len(ordered)),
            "dense_regime_knee_eps": dense_knee_eps,
            "dense_regime_rank_fraction": float(dense_knee / len(ordered)),
            "dense_regime_cap_fraction": 0.20,
            "knee_interpretation": "global knee is field-tail dominated; dense-regime knee brackets the association seed transition",
            "quantiles": {str(q): float(np.quantile(ordered, q)) for q in [0.5, 0.75, 0.9, 0.95, 0.975, 0.99]},
        },
        "outputs": {
            str(kdist_path.relative_to(ROOT)): sha256(kdist_path),
            str(output.relative_to(ROOT)): sha256(output),
            str(figure_path.relative_to(ROOT)): sha256(figure_path),
        },
    }
    write_json(PROVENANCE / "wp2_target_density_scan_execution.json", record)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stage", choices=["scan", "hdbscan-scan", "hdbscan-leaf-scan", "target-scan"], required=True)
    args = parser.parse_args()
    if args.stage == "scan":
        scan()
    elif args.stage == "hdbscan-scan":
        hdbscan_scan()
    elif args.stage == "hdbscan-leaf-scan":
        hdbscan_scan(cluster_selection_method="leaf")
    elif args.stage == "target-scan":
        target_field_scan()


if __name__ == "__main__":
    main()
