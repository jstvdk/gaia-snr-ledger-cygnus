#!/usr/bin/env python3
"""Shared WP5 machinery: branch metadata, injections, and IMF integrals.

WP5 consumes the frozen WP2--WP4 products.  It never edits an upstream
catalogue.  The completeness experiment treats every operation that can remove
a source from the mass catalogue as part of the selection function:

* the frozen Gaia query limits (G < 19 and 0.35 < raw parallax < 1.10 mas);
* the exact WP2 astrometric-quality filter;
* the exact WP2 cluster-versus-field classifier, reconstructed with the frozen
  seed and hyperparameters;
* the published P > 0.5 handoff used by WP3/WP4;
* availability of G, BP, and RP needed for a WP4 CMD mass.

Synthetic photometry is generated from the same PARSEC/MIST grids and fitted
subgroup ages used by WP4.  Extinction is resampled spatially from the matching
WP3 subgroup and R_V branch.  Gaia observational states are cloned from nearby,
similar-colour real sources in the one-degree association field.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from sklearn.cluster import DBSCAN
from sklearn.mixture import GaussianMixture

from wp2_finalize_membership import (
    CONTROL_CENTERS,
    FIELD_COMPONENTS,
    RADIUS_DEG,
    SELECTED_EPS,
    SELECTED_MIN_SAMPLES,
    TARGET_CENTER,
    circle_mask,
    circle_offsets,
    measurement_covariances,
)
from wp2_membership_pipeline import FEATURES, SEED as WP2_SEED, load_analysis
from wp3_extinction_law import R_V_BRANCHES, band_coefficients
from wp4_common import DIST_MODULUS, SUBGROUPS


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
PROVENANCE = ROOT / "provenance"
FIGURES = ROOT / "figures" / "wp5"
TABLES = ROOT / "tables"
NOTEBOOKS = ROOT / "notebooks"

SEED = 20260723
FAMILIES = ["PARSEC", "MIST"]
IMF_SLOPES = [2.0, 2.3, 2.6]
# Extend below the nominal 2-Msun calibration edge so the forward response can
# account for lower-mass stars scattered into the observed 2--8-Msun window by
# unresolved binaries and the WP3/WP4 inverse-mass uncertainty.
MASS_GRID = np.round(np.arange(0.50, 8.0001, 0.25), 2)
N_INJECT_PER_MASS = 400
F_BINARY = 0.40
Q_MIN = 0.10
MEMBERSHIP_QMC_POINTS = 128
MEMBERSHIP_THRESHOLD = 0.50
COMPLETENESS_TARGET = 0.95
CALIBRATION_NOMINAL_LO = 2.0
CALIBRATION_HI = 8.0
N_IMF_BINS = 6
N_POSTERIOR_DRAWS = 10_000
TOTAL_MASS_RANGE = (0.08, 120.0)
LOW_MASS_BREAK = 0.5
LOW_MASS_SLOPE = 1.3
LITERATURE_MASS_MSUN = 16_500.0
LITERATURE_FACTOR_GATE = 2.0


def sha256(path: Path | str) -> str:
    path = Path(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def branch_tag(family: str, rv: float) -> str:
    return f"{family}_rv{rv:.1f}"


def mass_column(family: str, rv: float) -> str:
    return f"mass_{family}_rv{rv:.1f}"


def fitted_age(subgroup: str, family: str, rv: float) -> float:
    """WP4 upper-MS MAP at f_bin=0.4 and the baseline distance."""
    posterior = pd.read_parquet(PROC / "wp4_age_posteriors.parquet")
    row = posterior[
        posterior["subgroup"].eq(subgroup)
        & posterior["family"].eq(family)
        & posterior["R_V"].eq(rv)
        & posterior["f_bin"].eq(F_BINARY)
        & posterior["indicator"].eq("ums")
        & posterior["dmu"].eq(0.0)
    ]
    if len(row) != 1:
        raise RuntimeError(
            f"expected one WP4 age row for {subgroup}/{family}/R_V={rv}; got {len(row)}"
        )
    return float(row["age_map"].iloc[0])


def load_isochrone_at_age(family: str, age_myr: float) -> tuple[pd.DataFrame, float]:
    path = PROC / f"wp3_isochrones_{family.lower()}.parquet"
    iso = pd.read_parquet(path)
    ages = np.sort(iso["age_Myr"].unique())
    nearest_age = float(ages[np.argmin(np.abs(ages - age_myr))])
    out = iso[np.isclose(iso["age_Myr"], nearest_age)].copy()
    # Below 8 Msun at 1--10 Myr the relevant locus is PMS/MS.  Removing later
    # phases also prevents interpolation across high-mass evolutionary loops.
    if family == "PARSEC":
        out = out[out["label"] <= 1]
    else:
        out = out[out["phase"] <= 1]
    out = out.sort_values("Mini").drop_duplicates("Mini", keep="first")
    required = ["Mini", "G0", "BP0", "RP0", "J0", "H0", "Ks0"]
    out = out.dropna(subset=required)
    return out, nearest_age


def interpolate_photometry(
    iso: pd.DataFrame,
    primary_mass: np.ndarray,
    secondary_mass: np.ndarray | None = None,
) -> dict[str, np.ndarray]:
    """Interpolate single or unresolved-binary absolute magnitudes."""
    result: dict[str, np.ndarray] = {}
    mini = iso["Mini"].to_numpy(float)
    for band in ["G", "BP", "RP", "J", "H", "Ks"]:
        primary = np.interp(primary_mass, mini, iso[f"{band}0"].to_numpy(float))
        if secondary_mass is None:
            result[band] = primary
            continue
        secondary = np.interp(
            np.clip(secondary_mass, mini.min(), mini.max()),
            mini,
            iso[f"{band}0"].to_numpy(float),
        )
        combined = -2.5 * np.log10(
            10.0 ** (-0.4 * primary) + 10.0 ** (-0.4 * secondary)
        )
        result[band] = np.where(secondary_mass > 0, combined, primary)
    return result


@dataclass
class WP2Classifier:
    frame: pd.DataFrame
    analysis: pd.DataFrame
    scaler: object
    cluster: dict[str, np.ndarray]
    field: dict[str, np.ndarray]
    prior: float
    reconstruction: dict


def reconstruct_wp2_classifier() -> WP2Classifier:
    """Rebuild the frozen deterministic WP2 observed-space mixture model."""
    frame, analysis, scaler, preprocessing = load_analysis()
    x = analysis[[f"scaled_{name}" for name in FEATURES]].to_numpy(float)
    centers = [TARGET_CENTER, *CONTROL_CENTERS]
    masks = [circle_mask(analysis, center) for center in centers]
    offsets = [
        circle_offsets(analysis, center)[mask]
        for center, mask in zip(centers, masks, strict=True)
    ]

    target_x = x[masks[0]]
    labels = DBSCAN(
        eps=SELECTED_EPS,
        min_samples=SELECTED_MIN_SAMPLES,
        n_jobs=-1,
    ).fit_predict(target_x)
    valid, counts = np.unique(labels[labels >= 0], return_counts=True)
    seed_label = valid[np.argmax(counts)]
    seed_mask = labels == seed_label
    seed_joint = np.column_stack([offsets[0][seed_mask], target_x[seed_mask]])

    cluster_models = [
        GaussianMixture(
            n_components,
            covariance_type="full",
            reg_covar=1e-5,
            random_state=WP2_SEED,
            n_init=3,
        ).fit(seed_joint)
        for n_components in [1, 2, 3, 4, 5, 6]
    ]
    cluster_bic = [float(model.bic(seed_joint)) for model in cluster_models]
    cluster_model = cluster_models[int(np.argmin(cluster_bic))]

    control_indices = np.concatenate([np.flatnonzero(mask) for mask in masks[1:]])
    control_x = x[control_indices]
    field_models = [
        GaussianMixture(
            n_components,
            covariance_type="full",
            reg_covar=1e-5,
            random_state=WP2_SEED,
            max_iter=500,
        ).fit(control_x)
        for n_components in FIELD_COMPONENTS
    ]
    field_bic = [float(model.bic(control_x)) for model in field_models]
    field_model = field_models[int(np.argmin(field_bic))]

    manifest = json.loads(
        (PROVENANCE / "wp2_membership_manifest.json").read_text(encoding="utf-8")
    )
    prior = float(
        manifest["mixture"]["cluster_prior_selected_by_max_control_leakage"]
    )
    expected_cluster = int(manifest["mixture"]["cluster_components_selected"])
    expected_field = int(manifest["mixture"]["field_components_selected"])
    if cluster_model.n_components != expected_cluster or field_model.n_components != expected_field:
        raise RuntimeError("reconstructed WP2 mixture does not match frozen manifest")

    reconstruction = {
        "preprocessing": preprocessing,
        "target_rows": int(masks[0].sum()),
        "control_rows": [int(mask.sum()) for mask in masks[1:]],
        "cluster_components": int(cluster_model.n_components),
        "field_components": int(field_model.n_components),
        "cluster_bic": cluster_bic,
        "field_bic": field_bic,
        "prior": prior,
        "matches_wp2_manifest": True,
    }
    return WP2Classifier(
        frame=frame,
        analysis=analysis,
        scaler=scaler,
        cluster={
            "weights": cluster_model.weights_,
            "means": cluster_model.means_,
            "covariances": cluster_model.covariances_,
        },
        field={
            "weights": field_model.weights_,
            "means": field_model.means_,
            "covariances": field_model.covariances_,
        },
        prior=prior,
        reconstruction=reconstruction,
    )


def _fixed_gmm_logpdf(values: np.ndarray, model: dict[str, np.ndarray]) -> np.ndarray:
    """Fast log density for draws with no per-row covariance."""
    dimensions = values.shape[1]
    constant = dimensions * np.log(2.0 * np.pi)
    terms = []
    for weight, mean, covariance in zip(
        model["weights"], model["means"], model["covariances"], strict=True
    ):
        sign, logdet = np.linalg.slogdet(covariance)
        if sign <= 0:
            raise RuntimeError("non-positive GMM covariance")
        delta = values - mean
        inverse = np.linalg.inv(covariance)
        mahalanobis = np.einsum("ni,ij,nj->n", delta, inverse, delta)
        terms.append(np.log(weight) - 0.5 * (constant + logdet + mahalanobis))
    return logsumexp(np.column_stack(terms), axis=1)


def qmc_membership_probabilities(
    x_scaled: np.ndarray,
    offsets: np.ndarray,
    measurement_covariance: np.ndarray,
    classifier: WP2Classifier,
    normal_points: np.ndarray,
    batch_size: int = 256,
) -> np.ndarray:
    """Approximate the exact WP2 10k-draw probability with fixed QMC points.

    This is the same hard posterior-odds classification as WP2.  The only
    approximation is replacing 10,000 pseudorandom draws per injected source
    with a deterministic 128-point Sobol normal design; its decision agreement
    is measured against the frozen WP2 probabilities before use.
    """
    output = np.zeros(len(x_scaled), dtype=float)
    n_points = len(normal_points)
    for start in range(0, len(x_scaled), batch_size):
        stop = min(start + batch_size, len(x_scaled))
        central = x_scaled[start:stop]
        covariance = measurement_covariance[start:stop]
        # PSD repair mirrors WP2; donor covariances should already be valid.
        cholesky = np.linalg.cholesky(covariance)
        draws = central[:, None, :] + np.einsum(
            "ni,bji->bnj", normal_points, cholesky
        )
        flat = draws.reshape(-1, 3)
        flat_offsets = np.repeat(offsets[start:stop], n_points, axis=0)
        joint = np.column_stack([flat_offsets, flat])
        cluster_log = _fixed_gmm_logpdf(joint, classifier.cluster)
        field_log = _fixed_gmm_logpdf(flat, classifier.field) - np.log(
            np.pi * RADIUS_DEG**2
        )
        selected = (
            np.log(classifier.prior) + cluster_log
            > np.log1p(-classifier.prior) + field_log
        ).reshape(stop - start, n_points)
        output[start:stop] = selected.mean(axis=1)
    return output


def power_integral(lo: float, hi: float, exponent: float) -> float:
    """Integral of m**exponent dm."""
    if np.isclose(exponent, -1.0):
        return float(np.log(hi / lo))
    return float((hi ** (exponent + 1.0) - lo ** (exponent + 1.0)) / (exponent + 1.0))


def imf_number_integral(lo: float, hi: float, alpha: float) -> float:
    return power_integral(lo, hi, -alpha)


def primary_system_mass_per_k(alpha: float) -> float:
    """Total primary-system mass for a high-mass normalization k."""
    lo, hi = TOTAL_MASS_RANGE
    low_coefficient = LOW_MASS_BREAK ** (LOW_MASS_SLOPE - alpha)
    low = low_coefficient * power_integral(lo, LOW_MASS_BREAK, 1.0 - LOW_MASS_SLOPE)
    high = power_integral(LOW_MASS_BREAK, hi, 1.0 - alpha)
    return low + high


def companion_mass_per_k(alpha: float, n_grid: int = 20_000) -> float:
    """Expected unresolved-companion mass for f_bin=0.4, q~U[0.1,1].

    A companion is counted only when q*m >= 0.08 Msun.  The convention matches
    the WP4 binary forward model and keeps k explicitly a *primary-system*
    normalization.  Both primary-system and multiplicity-adjusted masses are
    reported.
    """
    lo, hi = TOTAL_MASS_RANGE
    mass = np.geomspace(lo, hi, n_grid)
    high_density = mass ** (-alpha)
    low_coefficient = LOW_MASS_BREAK ** (LOW_MASS_SLOPE - alpha)
    density = np.where(
        mass < LOW_MASS_BREAK,
        low_coefficient * mass ** (-LOW_MASS_SLOPE),
        high_density,
    )
    q_lower = np.maximum(Q_MIN, lo / mass)
    mean_q_visible = np.where(
        q_lower < 1.0,
        (1.0 - q_lower**2) / (2.0 * (1.0 - Q_MIN)),
        0.0,
    )
    companion = F_BINARY * mass * mean_q_visible
    return float(np.trapezoid(density * companion, mass))
