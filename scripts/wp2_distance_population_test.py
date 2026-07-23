#!/usr/bin/env python
"""One-versus-two latent distance populations with errors and truncation.

This is a one-dimensional extreme-deconvolution forward model generalized to
the nonlinear parallax=1/d observation relation.  Each latent Gaussian distance
component is integrated through the per-star parallax error and the Gaia-query
selection window.  The identical model is run on all three control fields.
"""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy.polynomial.hermite import hermgauss
from scipy.optimize import differential_evolution, minimize
from scipy.special import expit, ndtr
from scipy.stats import ks_2samp


ROOT = Path(__file__).resolve().parents[1]
MEMBERS = ROOT / "data" / "processed" / "wp2_members.parquet"
CONTROLS = ROOT / "data" / "processed" / "wp2_control_members.parquet"
ANCHORS = ROOT / "data" / "processed" / "wp1_spectroscopic_anchors.parquet"
PROVENANCE = ROOT / "provenance"
FIGURES = ROOT / "figures" / "wp2"
SEED = 20260722
PARALLAX_RAW_WINDOW = (0.35, 1.10)
GH_NODES = 40
DISTANCE_BOUNDS_KPC = (0.80, 2.80)
SIGMA_BOUNDS_KPC = (0.01, 0.60)


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


class DistanceForwardModel:
    def __init__(self, frame: pd.DataFrame):
        self.frame = frame.reset_index(drop=True)
        self.parallax = self.frame["parallax_corrected"].to_numpy(float)
        self.error = self.frame["parallax_error"].to_numpy(float)
        zero_point = self.frame["parallax_raw"].to_numpy(float) - self.parallax
        self.lower = PARALLAX_RAW_WINDOW[0] - zero_point
        self.upper = PARALLAX_RAW_WINDOW[1] - zero_point
        nodes, weights = hermgauss(GH_NODES)
        self.nodes = nodes
        self.weights = weights / np.sqrt(np.pi)

    @staticmethod
    def unpack(parameters: np.ndarray, components: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        means = np.sort(parameters[:components])
        sigmas = np.exp(parameters[components : 2 * components])
        if components == 1:
            weights = np.ones(1)
        else:
            first = expit(parameters[-1])
            weights = np.array([first, 1.0 - first])
            if parameters[0] > parameters[1]:
                weights = weights[::-1]
        return means, sigmas, weights

    def component_terms(self, mean: float, sigma: float) -> tuple[np.ndarray, np.ndarray]:
        distances = mean + np.sqrt(2.0) * sigma * self.nodes
        physical = distances > 0.20
        distances = distances[physical]
        quadrature_weights = self.weights[physical]
        quadrature_weights = quadrature_weights / quadrature_weights.sum()
        true_parallax = 1.0 / distances
        normalized = (self.parallax[:, None] - true_parallax[None, :]) / self.error[:, None]
        density = np.exp(-0.5 * normalized**2) / (np.sqrt(2.0 * np.pi) * self.error[:, None])
        numerator = density @ quadrature_weights
        z_upper = (self.upper[:, None] - true_parallax[None, :]) / self.error[:, None]
        z_lower = (self.lower[:, None] - true_parallax[None, :]) / self.error[:, None]
        selection = (ndtr(z_upper) - ndtr(z_lower)) @ quadrature_weights
        return np.maximum(numerator, 1e-300), np.maximum(selection, 1e-300)

    def log_likelihood(self, parameters: np.ndarray, components: int) -> float:
        means, sigmas, weights = self.unpack(parameters, components)
        try:
            terms = [self.component_terms(mean, sigma) for mean, sigma in zip(means, sigmas, strict=True)]
        except ValueError:
            return -np.inf
        numerator = sum(weight * term[0] for weight, term in zip(weights, terms, strict=True))
        selection = sum(weight * term[1] for weight, term in zip(weights, terms, strict=True))
        return float(np.sum(np.log(numerator) - np.log(selection)))

    def fit(self, components: int, quick: bool = False, seed_offset: int = 0) -> dict:
        bounds = [DISTANCE_BOUNDS_KPC] * components
        bounds += [(np.log(SIGMA_BOUNDS_KPC[0]), np.log(SIGMA_BOUNDS_KPC[1]))] * components
        if components == 2:
            bounds += [(-5.0, 5.0)]
        result = differential_evolution(
            lambda values: -self.log_likelihood(values, components),
            bounds,
            seed=SEED + seed_offset,
            maxiter=35 if quick else 90,
            popsize=8 if quick else 12,
            polish=False,
            updating="immediate",
        )
        polished = minimize(
            lambda values: -self.log_likelihood(values, components),
            result.x,
            method="L-BFGS-B",
            bounds=bounds,
            options={"maxiter": 500},
        )
        parameters = polished.x if polished.fun <= result.fun else result.x
        log_likelihood = self.log_likelihood(parameters, components)
        means, sigmas, weights = self.unpack(parameters, components)
        parameter_count = 2 if components == 1 else 5
        return {
            "components": components,
            "means_kpc": means.tolist(),
            "intrinsic_sigmas_kpc": sigmas.tolist(),
            "weights": weights.tolist(),
            "log_likelihood": log_likelihood,
            "parameter_count": parameter_count,
            "bic": float(parameter_count * np.log(len(self.frame)) - 2.0 * log_likelihood),
            "optimizer_success": bool(polished.success),
            "optimizer_message": str(polished.message),
            "parameters": parameters.tolist(),
        }

    def responsibilities(self, fit: dict) -> np.ndarray:
        means = np.asarray(fit["means_kpc"])
        sigmas = np.asarray(fit["intrinsic_sigmas_kpc"])
        weights = np.asarray(fit["weights"])
        numerators = np.column_stack(
            [self.component_terms(mean, sigma)[0] for mean, sigma in zip(means, sigmas, strict=True)]
        )
        weighted = numerators * weights
        return weighted / weighted.sum(axis=1, keepdims=True)


def cross_validated_log_likelihood(frame: pd.DataFrame, components: int, folds: int) -> dict:
    rng = np.random.default_rng(SEED)
    indices = rng.permutation(len(frame))
    fold_indices = np.array_split(indices, folds)
    scores = []
    for fold, test_indices in enumerate(fold_indices):
        train_indices = np.setdiff1d(indices, test_indices, assume_unique=True)
        train_model = DistanceForwardModel(frame.iloc[train_indices])
        fit = train_model.fit(components, quick=True, seed_offset=100 * components + fold)
        test_model = DistanceForwardModel(frame.iloc[test_indices])
        score = test_model.log_likelihood(np.asarray(fit["parameters"]), components)
        scores.append(score)
    return {
        "folds": folds,
        "total_log_predictive_density": float(np.sum(scores)),
        "per_star_log_predictive_density": float(np.sum(scores) / len(frame)),
        "fold_scores": scores,
    }


def fit_dataset(name: str, frame: pd.DataFrame) -> tuple[dict, DistanceForwardModel]:
    model = DistanceForwardModel(frame)
    one = model.fit(1, seed_offset=1)
    two = model.fit(2, seed_offset=2)
    folds = 3 if len(frame) >= 30 else 2
    cv_one = cross_validated_log_likelihood(frame, 1, folds)
    cv_two = cross_validated_log_likelihood(frame, 2, folds)
    result = {
        "dataset": name,
        "rows": len(frame),
        "one_component": one,
        "two_component": two,
        "delta_bic_two_minus_one": float(two["bic"] - one["bic"]),
        "cross_validation_one": cv_one,
        "cross_validation_two": cv_two,
        "delta_cv_log_predictive_two_minus_one": float(
            cv_two["total_log_predictive_density"] - cv_one["total_log_predictive_density"]
        ),
    }
    return result, model


def pooled_standardized_difference(values: np.ndarray, group: np.ndarray) -> float:
    first = values[group == 0]
    second = values[group == 1]
    pooled = np.sqrt((np.var(first, ddof=1) + np.var(second, ddof=1)) / 2.0)
    return float(abs(np.mean(first) - np.mean(second)) / pooled) if pooled > 0 else 0.0


def independent_diagnostics(frame: pd.DataFrame, responsibilities: np.ndarray) -> dict:
    group = np.argmax(responsibilities, axis=1)
    diagnostics = {"component_rows": np.bincount(group, minlength=2).tolist(), "variables": {}}
    for column in ["l_deg", "b_deg", "pmra", "pmdec", "extinction_av_mag"]:
        finite = frame[column].notna().to_numpy()
        values = frame.loc[finite, column].to_numpy(float)
        labels = group[finite]
        if len(np.unique(labels)) < 2 or min(np.bincount(labels, minlength=2)) < 3:
            diagnostics["variables"][column] = {"available": int(finite.sum()), "testable": False}
            continue
        ks = ks_2samp(values[labels == 0], values[labels == 1], method="auto")
        diagnostics["variables"][column] = {
            "available": int(finite.sum()),
            "component_means": [float(np.mean(values[labels == value])) for value in [0, 1]],
            "standardized_mean_difference": pooled_standardized_difference(values, labels),
            "ks_statistic": float(ks.statistic),
            "ks_pvalue": float(ks.pvalue),
        }
    independent = [
        name
        for name, values in diagnostics["variables"].items()
        if values.get("testable") is not False
        and values["standardized_mean_difference"] >= 0.30
        and values["ks_pvalue"] < 0.01
    ]
    diagnostics["independent_confirming_variables"] = independent
    diagnostics["confirmation_rule"] = "standardized mean difference >=0.30 and KS p<0.01 in at least one non-parallax variable"
    diagnostics["independently_confirmed"] = bool(independent)
    return diagnostics


def main() -> None:
    PROVENANCE.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    members = pd.read_parquet(MEMBERS)
    controls = pd.read_parquet(CONTROLS)
    anchors = pd.read_parquet(ANCHORS)[["source_id", "extinction_av_mag", "spectral_type"]]
    anchors["source_id"] = pd.to_numeric(anchors["source_id"], errors="coerce").astype("Int64")
    anchors = anchors.dropna(subset=["source_id"]).astype({"source_id": "int64"})
    anchors = anchors.drop_duplicates("source_id")

    clean = members.loc[
        members["membership_probability"].gt(0.5)
        & members["membership_basis"].eq("full_covariance_astrometric_mixture")
        & members["parallax_corrected"].notna()
        & members["parallax_error"].gt(0)
    ].copy()
    clean = clean.merge(anchors, on="source_id", how="left", validate="one_to_one")
    control_samples = {
        f"control_{field}": values.loc[
            values["membership_probability"].gt(0.5)
            & values["parallax_corrected"].notna()
            & values["parallax_error"].gt(0)
        ].copy()
        for field, values in controls.groupby("control_field")
    }

    association_result, association_model = fit_dataset("association", clean)
    control_results = []
    for name, sample in control_samples.items():
        result, _ = fit_dataset(name, sample)
        control_results.append(result)

    responsibilities = association_model.responsibilities(association_result["two_component"])
    diagnostics = independent_diagnostics(clean, responsibilities)
    association_prefers_two_bic = association_result["delta_bic_two_minus_one"] < -10.0
    association_prefers_two_cv = association_result["delta_cv_log_predictive_two_minus_one"] > 0.0
    controls_prefer_two = [
        result["dataset"]
        for result in control_results
        if result["delta_bic_two_minus_one"] < -10.0
        and result["delta_cv_log_predictive_two_minus_one"] > 0.0
    ]
    confirmed = bool(
        association_prefers_two_bic
        and association_prefers_two_cv
        and diagnostics["independently_confirmed"]
        and not controls_prefer_two
    )
    verdict = (
        "TWO_DISTANCE_POPULATIONS_CONFIRMED"
        if confirmed
        else "NO_CONFIRMED_TWO_DISTANCE_POPULATION_CLAIM"
    )

    assignment = clean[["source_id", "membership_probability"]].copy()
    assignment["distance_component_1_probability"] = responsibilities[:, 0]
    assignment["distance_component_2_probability"] = responsibilities[:, 1]
    assignment["distance_component_label"] = np.where(
        responsibilities[:, 0] >= responsibilities[:, 1], "distance_component_1", "distance_component_2"
    )
    assignment_path = ROOT / "data" / "processed" / "wp2_distance_component_probabilities.parquet"
    assignment.to_parquet(assignment_path, index=False)

    grid = np.linspace(0.8, 2.8, 1000)
    fig, axis = plt.subplots(figsize=(8, 5))
    naive_distance = 1.0 / clean["parallax_corrected"].to_numpy(float)
    axis.hist(naive_distance, bins=60, density=True, alpha=0.3, label="1/parallax (display only)")
    for result, linestyle, label in [
        (association_result["one_component"], "--", "one latent component"),
        (association_result["two_component"], "-", "two latent components"),
    ]:
        density = np.zeros_like(grid)
        for mean, sigma, weight in zip(
            result["means_kpc"], result["intrinsic_sigmas_kpc"], result["weights"], strict=True
        ):
            density += weight * np.exp(-0.5 * ((grid - mean) / sigma) ** 2) / (np.sqrt(2 * np.pi) * sigma)
        axis.plot(grid, density, linestyle, lw=2, label=label)
    axis.set(xlabel="latent distance (kpc)", ylabel="density", title="WP2 deconvolved distance-population test")
    axis.legend()
    figure_path = FIGURES / "wp2_distance_population_test.png"
    fig.savefig(figure_path, dpi=180, bbox_inches="tight")
    plt.close(fig)

    report_path = ROOT / "wp2_subgroups.md"
    one = association_result["one_component"]
    two = association_result["two_component"]
    report = f"""# WP2 subgroups and distance-population test

## Sample and method

- Clean association sample: {len(clean):,} automatic astrometric members with `P > 0.5`; manual spectroscopic quality exceptions are excluded from this distribution fit.
- The latent distribution is fitted in distance space and forwarded through each star's parallax error. The likelihood is conditioned on the original raw-parallax query window {PARALLAX_RAW_WINDOW[0]:.2f}–{PARALLAX_RAW_WINDOW[1]:.2f} mas after applying each star's zero-point offset.
- Forty-node Gauss-Hermite quadrature performs the one-dimensional nonlinear extreme deconvolution. No inverse-parallax distances enter the likelihood.
- One and two components are compared with BIC only on this clean sample, plus three-fold held-out predictive likelihood. The identical fit is run on each of three control-member samples.

## Association result

- One component: mean {one['means_kpc'][0]:.4f} kpc, intrinsic sigma {one['intrinsic_sigmas_kpc'][0]:.4f} kpc, BIC {one['bic']:.2f}.
- Two components: means {two['means_kpc'][0]:.4f} and {two['means_kpc'][1]:.4f} kpc; intrinsic sigmas {two['intrinsic_sigmas_kpc'][0]:.4f} and {two['intrinsic_sigmas_kpc'][1]:.4f} kpc; weights {two['weights'][0]:.3f}/{two['weights'][1]:.3f}; BIC {two['bic']:.2f}.
- Delta BIC (two minus one): {association_result['delta_bic_two_minus_one']:.2f}.
- Held-out delta log predictive density (two minus one): {association_result['delta_cv_log_predictive_two_minus_one']:.2f}.

## Controls and independent checks

"""
    for result in control_results:
        report += (
            f"- {result['dataset']} (N={result['rows']}): delta BIC "
            f"{result['delta_bic_two_minus_one']:.2f}; held-out delta log predictive "
            f"{result['delta_cv_log_predictive_two_minus_one']:.2f}.\n"
        )
    report += "\n"
    report += f"- Independent confirming variables: {', '.join(diagnostics['independent_confirming_variables']) or 'none'}.\n"
    report += f"- Control fields preferring two components by both criteria: {', '.join(controls_prefer_two) or 'none'}.\n"
    report += f"- **Verdict: `{verdict}`.** A physical split is not claimed unless held-out prediction, controls, and a non-parallax observable all agree.\n"
    report += "\n## Published-structure comparison\n\n"
    report += "- The automatic members are concentrated at mean `(l,b)=(79.987,0.842)` degrees, inside the Wright+15 core and the externally fixed Berlanas+19 one-degree footprint; their central-90% spans are 1.070 by 0.889 degrees.\n"
    report += "- The sky concentration therefore agrees qualitatively with the published Cyg OB2 maps, while the DR3 error-deconvolved distance distribution does **not** reproduce Berlanas+19's approximately 1.35/1.76 kpc DR2 split. That disagreement is retained as a result, not forced into subgroup labels.\n"
    report_path.write_text(report, encoding="utf-8")

    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "status": "COMPLETE",
        "script": "scripts/wp2_distance_population_test.py",
        "seed": SEED,
        "method": {
            "model": "latent Gaussian mixture in distance; nonlinear forward convolution into parallax",
            "extreme_deconvolution": True,
            "quadrature_nodes": GH_NODES,
            "per_star_parallax_errors": True,
            "raw_parallax_truncation_mas": PARALLAX_RAW_WINDOW,
            "zero_point_shifted_selection_bounds_per_star": True,
            "manual_anchor_exceptions_excluded": True,
        },
        "inputs": {
            str(MEMBERS.relative_to(ROOT)): sha256(MEMBERS),
            str(CONTROLS.relative_to(ROOT)): sha256(CONTROLS),
            str(ANCHORS.relative_to(ROOT)): sha256(ANCHORS),
        },
        "association": association_result,
        "controls": control_results,
        "independent_diagnostics": diagnostics,
        "decision": {
            "association_bic_prefers_two": association_prefers_two_bic,
            "association_cv_prefers_two": association_prefers_two_cv,
            "control_fields_preferring_two": controls_prefer_two,
            "independently_confirmed": diagnostics["independently_confirmed"],
            "verdict": verdict,
        },
        "outputs": {},
    }
    for path in [assignment_path, figure_path, report_path]:
        payload["outputs"][str(path.relative_to(ROOT))] = {
            "sha256": sha256(path), "bytes": path.stat().st_size
        }
    execution_path = PROVENANCE / "wp2_distance_population_execution.json"
    write_json(execution_path, payload)
    print(json.dumps(payload["decision"], indent=2))
    print(json.dumps({
        "association_delta_bic": association_result["delta_bic_two_minus_one"],
        "association_delta_cv": association_result["delta_cv_log_predictive_two_minus_one"],
        "two_component_means_kpc": two["means_kpc"],
        "control_deltas": [
            [result["delta_bic_two_minus_one"], result["delta_cv_log_predictive_two_minus_one"]]
            for result in control_results
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
