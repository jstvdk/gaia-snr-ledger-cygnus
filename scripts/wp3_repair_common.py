#!/usr/bin/env python3
"""Shared estimator for the versioned WP3 extinction repair.

The repair differs from the frozen WP3 estimator in two deliberate ways:

* a 0.03-mag calibration/model floor is added to every available band;
* the 149 spectroscopic extinction anchors define a local Gaussian prior whose
  width is measured by leave-one-out eight-neighbour residuals.

The fitted object is a probability density on a fixed A_V grid.  Keeping this
density, rather than only its mean and a local curvature error, preserves the
PMS/ZAMS modes needed by the repaired WP4 mass inference and WP5 injections.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.signal import find_peaks
from sklearn.neighbors import NearestNeighbors

from wp3_common import BANDS, DIST_MODULUS
from wp3_extinction_law import R_V_BRANCHES, band_coefficients


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"

REPAIR_VERSION = "repair_v1"
PHOTOMETRIC_FLOOR_MAG = 0.03
ANCHOR_NEIGHBOURS = 8
AV_GRID = np.linspace(0.0, 15.0, 301)
AV_GRID_STEP = float(AV_GRID[1] - AV_GRID[0])
MIN_PRIOR_SIGMA_MAG = 0.30
MAX_TEMPLATES_PER_FAMILY_AGE = 64


def _template_quadrature_weights(frame: pd.DataFrame) -> np.ndarray:
    """Sampling-density correction with equal family and age-cell weight.

    Isochrone files have different and nonuniform native mass samplings.
    Treating every row equally therefore creates an accidental prior.  Within
    each family/age cell we integrate with d(log M), then give every age cell
    equal weight and both model families equal total weight.
    """
    weight = np.zeros(len(frame), dtype=float)
    family_count = frame["family"].nunique()
    for family, family_frame in frame.groupby("family"):
        age_groups = list(family_frame.groupby("logAge"))
        for _, group in age_groups:
            order = np.argsort(group["Mini"].to_numpy(float))
            index = group.index.to_numpy()[order]
            log_mass = np.log(np.clip(group.loc[index, "Mini"].to_numpy(float), 1e-4, None))
            if len(log_mass) == 1:
                cell = np.ones(1)
            else:
                cell = np.gradient(log_mass)
                cell = np.clip(np.abs(cell), 1e-8, None)
            cell /= cell.sum()
            weight[index] = cell / (family_count * len(age_groups))
    weight /= weight.sum()
    return weight


def load_template_library() -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    columns = ["logAge", "age_Myr", "Mini", "Mass", "logTe", "family"]
    columns += [f"{band}0" for band in BANDS]
    frames = []
    for family in ["PARSEC", "MIST"]:
        frame = pd.read_parquet(PROC / f"wp3_isochrones_{family.lower()}.parquet")
        frame = frame[[column for column in columns if column in frame.columns]].copy()
        frame["family"] = family
        frames.append(frame)
    templates = pd.concat(frames, ignore_index=True)
    templates = templates.dropna(subset=[f"{band}0" for band in BANDS]).reset_index(drop=True)
    sampled = []
    for _, group in templates.groupby(["family", "logAge"], sort=False):
        group = group.sort_values("Mini")
        if len(group) > MAX_TEMPLATES_PER_FAMILY_AGE:
            index = np.unique(
                np.round(
                    np.linspace(0, len(group) - 1, MAX_TEMPLATES_PER_FAMILY_AGE)
                ).astype(int)
            )
            group = group.iloc[index]
        sampled.append(group)
    templates = pd.concat(sampled, ignore_index=True)
    magnitudes = np.column_stack(
        [templates[f"{band}0"].to_numpy(float) for band in BANDS]
    )
    prior_weight = _template_quadrature_weights(templates)
    return templates, magnitudes, prior_weight


@dataclass
class AnchorMap:
    anchors: pd.DataFrame
    neighbour_model: NearestNeighbors
    cos_b0: float
    prior_sigma: dict[float, float]
    variogram: dict[float, dict[str, float]]

    @classmethod
    def from_frozen_wp3(cls) -> "AnchorMap":
        extinction = pd.read_parquet(PROC / "wp3_extinction.parquet")
        anchors = extinction[
            extinction["av_method"].eq("intrinsic_color_spectroscopic")
            & extinction["l_deg"].notna()
            & extinction["b_deg"].notna()
        ].copy()
        cos_b0 = float(np.cos(np.deg2rad(anchors["b_deg"].median())))
        coordinates = np.column_stack(
            [anchors["l_deg"].to_numpy(float) * cos_b0, anchors["b_deg"].to_numpy(float)]
        )
        model = NearestNeighbors(n_neighbors=ANCHOR_NEIGHBOURS + 1).fit(coordinates)
        distances, indices = model.kneighbors(coordinates)
        prior_sigma: dict[float, float] = {}
        variogram: dict[float, dict[str, float]] = {}
        for rv in R_V_BRANCHES:
            values = anchors[f"av_rv{rv:.1f}"].to_numpy(float)
            local = np.nanmedian(values[indices[:, 1:]], axis=1)
            residual = values - local
            centre = float(np.nanmedian(residual))
            robust_sigma = float(1.4826 * np.nanmedian(np.abs(residual - centre)))
            sigma = max(MIN_PRIOR_SIGMA_MAG, robust_sigma)
            prior_sigma[float(rv)] = sigma
            variogram[float(rv)] = {
                "n_anchors": int(np.isfinite(values).sum()),
                "n_neighbours": ANCHOR_NEIGHBOURS,
                "median_eighth_neighbour_separation_deg": float(np.nanmedian(distances[:, -1])),
                "leave_one_out_residual_median_mag": centre,
                "leave_one_out_robust_sigma_mag": robust_sigma,
                "leave_one_out_rms_mag": float(np.sqrt(np.nanmean(residual**2))),
                "adopted_prior_sigma_mag": sigma,
            }
        return cls(anchors, model, cos_b0, prior_sigma, variogram)

    def evaluate(
        self, longitude: np.ndarray, latitude: np.ndarray, rv: float
    ) -> tuple[np.ndarray, np.ndarray]:
        longitude = np.asarray(longitude, dtype=float)
        latitude = np.asarray(latitude, dtype=float)
        centre = np.full(len(longitude), self.anchors[f"av_rv{rv:.1f}"].median(), dtype=float)
        separation = np.full(len(longitude), np.nan, dtype=float)
        valid = np.isfinite(longitude) & np.isfinite(latitude)
        if valid.any():
            coordinates = np.column_stack(
                [longitude[valid] * self.cos_b0, latitude[valid]]
            )
            distances, indices = self.neighbour_model.kneighbors(
                coordinates, n_neighbors=ANCHOR_NEIGHBOURS
            )
            anchor_values = self.anchors[f"av_rv{rv:.1f}"].to_numpy(float)
            centre[valid] = np.nanmedian(anchor_values[indices], axis=1)
            separation[valid] = distances[:, -1]
        return centre, separation


def _grid_summary(probability: np.ndarray) -> dict[str, float | int]:
    probability = np.asarray(probability, dtype=float)
    probability = probability / probability.sum()
    cdf = np.cumsum(probability)
    quantile = lambda q: float(np.interp(q, cdf, AV_GRID))
    mean = float(np.sum(probability * AV_GRID))
    sd = float(np.sqrt(np.sum(probability * (AV_GRID - mean) ** 2)))
    peaks, properties = find_peaks(probability, prominence=0.02 * probability.max())
    if len(peaks) == 0:
        peaks = np.array([int(np.argmax(probability))])
        prominences = np.array([probability[peaks[0]]])
    else:
        prominences = properties["prominences"]
    order = np.argsort(probability[peaks])[::-1]
    peaks = peaks[order]
    prominences = prominences[order]
    return {
        "av_mean": mean,
        "av_sd": sd,
        "av_q16": quantile(0.16),
        "av_q50": quantile(0.50),
        "av_q84": quantile(0.84),
        "av_map": float(AV_GRID[int(np.argmax(probability))]),
        "n_modes": int(len(peaks)),
        "mode1_av": float(AV_GRID[peaks[0]]),
        "mode1_probability_density": float(probability[peaks[0]] / AV_GRID_STEP),
        "mode2_av": float(AV_GRID[peaks[1]]) if len(peaks) > 1 else np.nan,
        "mode2_relative_height": (
            float(probability[peaks[1]] / probability[peaks[0]])
            if len(peaks) > 1
            else 0.0
        ),
        "peak_prominence_sum": float(np.sum(prominences)),
    }


def gaussian_grid(mean: float, sigma: float) -> np.ndarray:
    sigma = max(float(sigma), AV_GRID_STEP / 2.0)
    probability = np.exp(-0.5 * ((AV_GRID - mean) / sigma) ** 2)
    probability /= probability.sum()
    return probability


def fit_extinction_posterior(
    magnitudes: np.ndarray,
    errors: np.ndarray,
    rv: float,
    prior_mean: float,
    prior_sigma: float,
    template_magnitudes: np.ndarray,
    template_weights: np.ndarray,
    template_branch_sigma: float = 0.0,
) -> tuple[np.ndarray | None, dict[str, float | int]]:
    """Fit one star and return its complete gridded A_V marginal posterior."""
    magnitudes = np.asarray(magnitudes, dtype=float)
    errors = np.asarray(errors, dtype=float)
    available = np.isfinite(magnitudes) & np.isfinite(errors) & (errors > 0)
    if int(available.sum()) < 2:
        return None, {"n_bands_fit": int(available.sum())}

    coefficient = np.array([band_coefficients(rv)[band] for band in BANDS])
    effective_error = np.sqrt(errors**2 + PHOTOMETRIC_FLOOR_MAG**2)
    weight = np.zeros_like(effective_error)
    weight[available] = 1.0 / effective_error[available] ** 2
    y = np.where(available, magnitudes - DIST_MODULUS, 0.0)
    denominator = float(np.sum(weight * coefficient**2))
    weighted_k = weight * coefficient
    numerator = (
        np.sum(y * weighted_k)
        - np.sum(template_magnitudes * weighted_k[None, :], axis=1)
    )
    av_conditional = np.clip(numerator / denominator, AV_GRID[0], AV_GRID[-1])
    residual = y[None, :] - template_magnitudes - av_conditional[:, None] * coefficient
    chi2 = np.sum(weight[None, :] * residual**2, axis=1)
    log_weight = -0.5 * (chi2 - np.nanmin(chi2)) + np.log(template_weights)
    mixture_weight = np.exp(log_weight - np.max(log_weight))
    mixture_weight /= mixture_weight.sum()

    edges = np.r_[
        AV_GRID[0] - AV_GRID_STEP / 2.0,
        AV_GRID[:-1] + AV_GRID_STEP / 2.0,
        AV_GRID[-1] + AV_GRID_STEP / 2.0,
    ]
    probability, _ = np.histogram(
        av_conditional, bins=edges, weights=mixture_weight
    )
    conditional_sigma = np.sqrt(1.0 / denominator)
    probability = gaussian_filter1d(
        probability.astype(float),
        sigma=max(conditional_sigma / AV_GRID_STEP, 0.5),
        mode="constant",
    )
    if template_branch_sigma > 0:
        probability = gaussian_filter1d(
            probability,
            sigma=template_branch_sigma / AV_GRID_STEP,
            mode="constant",
        )
    spatial_prior = np.exp(-0.5 * ((AV_GRID - prior_mean) / prior_sigma) ** 2)
    probability *= spatial_prior
    if not np.isfinite(probability).all() or probability.sum() <= 0:
        probability = gaussian_grid(prior_mean, prior_sigma)
    else:
        probability /= probability.sum()
    summary = _grid_summary(probability)
    summary.update(
        {
            "n_bands_fit": int(available.sum()),
            "conditional_av_sigma": float(conditional_sigma),
            "chi2_min": float(np.nanmin(chi2)),
            "effective_error_floor_mag": PHOTOMETRIC_FLOOR_MAG,
            "prior_mean": float(prior_mean),
            "prior_sigma": float(prior_sigma),
            "template_branch_sigma": float(template_branch_sigma),
        }
    )
    return probability, summary
