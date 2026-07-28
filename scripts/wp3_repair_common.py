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

import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from sklearn.neighbors import NearestNeighbors

from wp3_common import BANDS, DIST_MODULUS
from wp3_extinction_law import R_V_BRANCHES, band_coefficients


ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"

# Both knobs are environment-overridable so the whole repair chain (WP3 -> WP4
# -> WP5) can be re-run at a new version without editing any script.  Defaults
# reproduce repair_v1 exactly.
REPAIR_VERSION = os.environ.get("WP_REPAIR_VERSION", "repair_v1")

# "global"    -- one prior width per R_V branch, from leave-one-out at the
#                anchor density.  The repair_v1 behaviour.
# "variogram" -- per-star width read off the fitted anchor variogram at that
#                star's own distance to the eighth-nearest anchor.  The global
#                width is calibrated at a median anchor separation of 0.071 deg
#                but applied to CygOB2-B members sitting 0.377 deg away, where
#                the A_V field has decorrelated; it understates the true
#                predictive width there by 2.3x and collapses that subgroup's
#                differential extinction 9.5x.  See
#                provenance/wp3_anchor_prior_diagnostic_execution.json and
#                reports/WP3_ANCHOR_PRIOR_REPAIR_repair_v3.md.
ANCHOR_PRIOR_MODE = os.environ.get("WP3_ANCHOR_PRIOR_MODE", "global")

PHOTOMETRIC_FLOOR_MAG = 0.03
ANCHOR_NEIGHBOURS = 8
AV_GRID = np.linspace(0.0, 15.0, 301)
AV_GRID_STEP = float(AV_GRID[1] - AV_GRID[0])
MIN_PRIOR_SIGMA_MAG = 0.30
# Tikhonov jitter for the 8x8 kriging solve; the fitted nugget is ~0 so the
# covariance matrix is near-singular when two anchors nearly coincide.
KRIGING_JITTER = 1e-6
MAX_TEMPLATES_PER_FAMILY_AGE = 64
# Separation bins and minimum pair count for the anchor variogram, matched to
# scripts/wp3_anchor_prior_diagnostic.py so the two always agree.
VARIOGRAM_BINS_DEG = np.array(
    [0.0, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.70, 1.00, 2.00]
)
VARIOGRAM_MIN_PAIRS = 20


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


def exponential_variogram(
    separation: np.ndarray, nugget: float, sill: float, correlation_range: float
) -> np.ndarray:
    """Standard exponential variogram, returned as a sigma (not a variance)."""
    variance = nugget**2 + (sill**2 - nugget**2) * (
        1.0 - np.exp(-3.0 * np.asarray(separation, dtype=float) / correlation_range)
    )
    return np.sqrt(np.maximum(variance, 1e-12))


def _fit_anchor_variogram(
    pair_separation: np.ndarray, value: np.ndarray, n_anchor: int
) -> dict[str, float]:
    """Fit the exponential variogram of the anchor A_V field.

    Gives the predictive width of the anchor-median prior as a function of how
    far a star sits from the anchors that inform it.
    """
    pair_difference = (value[:, None] - value[None, :])[
        np.triu_indices(n_anchor, 1)
    ]
    finite = np.isfinite(pair_difference) & np.isfinite(pair_separation)
    separation, difference = pair_separation[finite], pair_difference[finite]
    centres, sigmas, counts = [], [], []
    for low, high in zip(VARIOGRAM_BINS_DEG[:-1], VARIOGRAM_BINS_DEG[1:], strict=True):
        inside = (separation >= low) & (separation < high)
        if int(inside.sum()) < VARIOGRAM_MIN_PAIRS:
            continue
        residual = difference[inside]
        centres.append(float(np.median(separation[inside])))
        # A difference of two independent draws carries sqrt(2) x the
        # single-star scatter, so divide it back out.
        sigmas.append(
            float(1.4826 * np.median(np.abs(residual - np.median(residual))))
            / np.sqrt(2.0)
        )
        counts.append(int(inside.sum()))
    parameters, _ = curve_fit(
        exponential_variogram,
        np.array(centres),
        np.array(sigmas),
        p0=[0.35, 1.05, 0.20],
        sigma=1.0 / np.sqrt(np.array(counts, dtype=float)),
        bounds=([0.0, 0.0, 0.01], [3.0, 5.0, 5.0]),
        maxfev=20000,
    )
    return {
        "nugget_mag": float(parameters[0]),
        "sill_mag": float(parameters[1]),
        "correlation_range_deg": float(parameters[2]),
        "n_bins_fitted": len(centres),
    }


@dataclass
class AnchorMap:
    anchors: pd.DataFrame
    neighbour_model: NearestNeighbors
    cos_b0: float
    prior_sigma: dict[float, float]
    variogram: dict[float, dict[str, float]]
    variogram_model: dict[float, dict[str, float]] = field(default_factory=dict)

    def prior_sigma_at(self, separation: np.ndarray, rv: float) -> np.ndarray:
        """Prior width for each star, honouring ANCHOR_PRIOR_MODE.

        In "global" mode this is the repair_v1 scalar broadcast to every star.
        In "variogram" mode it is the fitted anchor variogram evaluated at each
        star's own eighth-nearest-anchor separation: stars sitting inside the
        anchor field keep a tight prior, stars extrapolated beyond the
        correlation length correctly receive one near the sill, which lets their
        own photometry rather than distant anchors set their A_V.
        """
        separation = np.asarray(separation, dtype=float)
        scalar = float(self.prior_sigma[float(rv)])
        if ANCHOR_PRIOR_MODE == "global":
            return np.full(separation.shape, scalar, dtype=float)
        if ANCHOR_PRIOR_MODE not in ("variogram", "kriging"):
            raise ValueError(f"unknown ANCHOR_PRIOR_MODE {ANCHOR_PRIOR_MODE!r}")
        model = self.variogram_model[float(rv)]
        sigma = exponential_variogram(
            np.where(np.isfinite(separation), separation, np.inf),
            model["nugget_mag"],
            model["sill_mag"],
            model["correlation_range_deg"],
        )
        # A star with no usable position has no spatial information at all, so
        # it gets the sill rather than the tight calibrated width.
        sigma = np.where(np.isfinite(separation), sigma, model["sill_mag"])
        return np.clip(sigma, MIN_PRIOR_SIGMA_MAG, None)

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
        variogram_model: dict[float, dict[str, float]] = {}
        pair_separation = np.sqrt(
            (coordinates[:, 0][:, None] - coordinates[:, 0][None, :]) ** 2
            + (coordinates[:, 1][:, None] - coordinates[:, 1][None, :]) ** 2
        )[np.triu_indices(len(coordinates), 1)]
        for rv in R_V_BRANCHES:
            variogram_model[float(rv)] = _fit_anchor_variogram(
                pair_separation,
                anchors[f"av_rv{rv:.1f}"].to_numpy(float),
                len(coordinates),
            )
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
        return cls(anchors, model, cos_b0, prior_sigma, variogram, variogram_model)

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
            if ANCHOR_PRIOR_MODE == "kriging":
                centre[valid] = self._kriged_mean(
                    coordinates, distances, indices, anchor_values, rv
                )
            else:
                centre[valid] = np.nanmedian(anchor_values[indices], axis=1)
            separation[valid] = distances[:, -1]
        return centre, separation

    def _kriged_mean(
        self,
        coordinates: np.ndarray,
        distances: np.ndarray,
        indices: np.ndarray,
        anchor_values: np.ndarray,
        rv: float,
    ) -> np.ndarray:
        """Simple-kriging prior mean using the already-fitted variogram.

        The plain neighbour median gives the eight nearest anchors weights
        summing to one no matter how far away they are, so a star extrapolated
        far beyond the correlation range is still centred on distant anchors at
        full strength -- even though ``prior_sigma_at`` simultaneously widens
        that star's uncertainty to near the sill, i.e. declares those same
        anchors nearly uninformative.  The first and second moments of the
        prior therefore disagree.

        Simple kriging removes the inconsistency using the same fitted
        exponential model and introduces no new parameter: weights fall off
        with separation and the deficit is taken up by the field mean, so the
        estimate degrades gracefully to "no local information" instead of
        asserting a distant measurement.  Where anchors are genuinely local the
        weights sum to one and the estimate reduces to the neighbour average,
        leaving well-anchored stars untouched.

        The prior *width* is deliberately left at the unconditional variogram
        sigma rather than the (smaller) kriging variance: this change is meant
        to correct where the prior sits, not to make it more confident.
        """
        model = self.variogram_model[float(rv)]
        sill = float(model["sill_mag"])
        correlation_range = float(model["correlation_range_deg"])
        field_mean = float(np.nanmedian(anchor_values))
        anchor_coordinates = np.column_stack(
            [
                self.anchors["l_deg"].to_numpy(float) * self.cos_b0,
                self.anchors["b_deg"].to_numpy(float),
            ]
        )

        def covariance(separation: np.ndarray) -> np.ndarray:
            return sill**2 * np.exp(-3.0 * separation / correlation_range)

        output = np.empty(len(coordinates), dtype=float)
        for row in range(len(coordinates)):
            neighbour = indices[row]
            values = anchor_values[neighbour]
            usable = np.isfinite(values)
            if usable.sum() == 0:
                output[row] = field_mean
                continue
            neighbour = neighbour[usable]
            values = values[usable]
            offsets = anchor_coordinates[neighbour]
            pair = np.sqrt(
                ((offsets[:, None, :] - offsets[None, :, :]) ** 2).sum(axis=-1)
            )
            matrix = covariance(pair) + np.eye(len(neighbour)) * KRIGING_JITTER
            try:
                weights = np.linalg.solve(matrix, covariance(distances[row][usable]))
            except np.linalg.LinAlgError:
                output[row] = float(np.nanmedian(values))
                continue
            output[row] = field_mean + float(weights @ (values - field_mean))
        return output


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
