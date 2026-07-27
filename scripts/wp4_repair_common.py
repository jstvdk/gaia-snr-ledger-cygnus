#!/usr/bin/env python3
"""Posterior-propagation machinery for the versioned WP4 mass repair."""
from __future__ import annotations

from functools import lru_cache

import numpy as np
import pandas as pd
from scipy.stats import norm

from wp3_extinction_law import band_coefficients
from wp3_repair_common import AV_GRID, BANDS, PHOTOMETRIC_FLOOR_MAG, PROC
from wp4_common import DIST_MODULUS


MASS_POSTERIOR_DRAWS = 256
MASS_GRID_SIZE = 512
MASS_GRID = np.geomspace(0.10, 120.0, MASS_GRID_SIZE)
F_BINARY = 0.40
Q_MIN = 0.10
N_Q_COMPONENTS = 9
N_AGE_NODES = 9
# Class C/D: median robust per-band residual (0.35--0.39 mag) for 107
# non-extreme spectroscopic-HRD anchors against their repaired-age PARSEC SED.
# The larger RMS tail is handled by the explicit unresolved-binary mixture.
MODEL_BAND_SCALE = 0.38
RANDOM_SEED = 20260724


def _flux_add(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    return -2.5 * np.log10(
        10.0 ** (-0.4 * first) + 10.0 ** (-0.4 * second)
    )


@lru_cache(maxsize=2)
def _isochrones(family: str) -> pd.DataFrame:
    frame = pd.read_parquet(PROC / f"wp3_isochrones_{family.lower()}.parquet")
    if family == "PARSEC":
        frame = frame[frame["label"] <= 1]
    else:
        frame = frame[frame["phase"] <= 1]
    return frame


def age_posterior_nodes(
    posterior_rows: pd.DataFrame, subgroup: str, family: str, rv: float
) -> np.ndarray:
    row = posterior_rows[
        posterior_rows["subgroup"].eq(subgroup)
        & posterior_rows["family"].eq(family)
        & posterior_rows["R_V"].eq(rv)
        & posterior_rows["f_bin"].eq(F_BINARY)
        & posterior_rows["indicator"].eq("ums")
        & posterior_rows["dmu"].eq(0.0)
    ]
    if len(row) != 1:
        fallback = posterior_rows[
            posterior_rows["family"].eq(family)
            & posterior_rows["R_V"].eq(rv)
            & posterior_rows["f_bin"].eq(F_BINARY)
            & posterior_rows["indicator"].eq("ums")
            & posterior_rows["dmu"].eq(0.0)
        ]
        if len(fallback) == 0:
            raise RuntimeError(f"missing age posterior for {subgroup}/{family}/R_V={rv}")
        centre = float(fallback["age_map"].median())
        lo = float(fallback["age_lo68"].median())
        hi = float(fallback["age_hi68"].median())
    else:
        centre = float(row["age_map"].iloc[0])
        lo = float(row["age_lo68"].iloc[0])
        hi = float(row["age_hi68"].iloc[0])
    probabilities = (np.arange(N_AGE_NODES) + 0.5) / N_AGE_NODES
    z = norm.ppf(probabilities)
    lower_scale = max(centre - lo, 0.02) / abs(norm.ppf(0.16))
    upper_scale = max(hi - centre, 0.02) / norm.ppf(0.84)
    nodes = centre + np.where(z < 0, lower_scale * z, upper_scale * z)
    return np.clip(nodes, 1.0, 10.0)


def _interpolate_age_sequence(
    family: str, age_myr: float
) -> tuple[np.ndarray, dict[str, np.ndarray]]:
    frame = _isochrones(family)
    ages = np.sort(frame["age_Myr"].unique())
    upper_index = int(np.searchsorted(ages, age_myr, side="left"))
    upper_index = min(max(upper_index, 0), len(ages) - 1)
    lower_index = max(upper_index - 1, 0)
    lower_age = float(ages[lower_index])
    upper_age = float(ages[upper_index])
    if np.isclose(lower_age, upper_age):
        fraction = 0.0
    else:
        fraction = float((age_myr - lower_age) / (upper_age - lower_age))

    age_frames = []
    for native_age in [lower_age, upper_age]:
        native = frame[np.isclose(frame["age_Myr"], native_age)].copy()
        native = native.sort_values("Mini").drop_duplicates("Mini", keep="first")
        native = native.dropna(
            subset=["Mini"] + [f"{band}0" for band in BANDS]
        )
        age_frames.append(native)
    minimum = max(float(value["Mini"].min()) for value in age_frames)
    maximum = min(float(value["Mini"].max()) for value in age_frames)
    mass = MASS_GRID[(MASS_GRID >= minimum) & (MASS_GRID <= maximum)]
    magnitudes: dict[str, np.ndarray] = {}
    for band in BANDS:
        values = []
        for native in age_frames:
            values.append(
                np.interp(
                    mass,
                    native["Mini"].to_numpy(float),
                    native[f"{band}0"].to_numpy(float),
                )
            )
        magnitudes[band] = (1.0 - fraction) * values[0] + fraction * values[1]
    return mass, magnitudes


@lru_cache(maxsize=512)
def model_sequence(
    family: str, age_myr_rounded: float, q_component: int
) -> tuple[np.ndarray, np.ndarray]:
    mass, magnitude = _interpolate_age_sequence(family, float(age_myr_rounded))
    if q_component == 0:
        combined = magnitude
    else:
        q_grid = np.linspace(Q_MIN, 1.0, N_Q_COMPONENTS)
        q = float(q_grid[q_component - 1])
        secondary_mass = q * mass
        secondary = {}
        for band in BANDS:
            secondary[band] = np.interp(
                np.clip(secondary_mass, mass.min(), mass.max()),
                mass,
                magnitude[band],
            )
        combined = {
            band: _flux_add(magnitude[band], secondary[band])
            for band in BANDS
        }
    return np.log(mass), np.column_stack([combined[band] for band in BANDS])


def draw_from_grid(
    probability: np.ndarray, rng: np.random.Generator, n_draws: int
) -> np.ndarray:
    probability = np.asarray(probability, dtype=float)
    if not np.isfinite(probability).all() or probability.sum() <= 0:
        return np.full(n_draws, np.nan)
    probability /= probability.sum()
    cdf = np.cumsum(probability)
    uniform = (np.arange(n_draws) + rng.random(n_draws)) / n_draws
    rng.shuffle(uniform)
    return np.interp(uniform, cdf, AV_GRID)


def _draw_mass_likelihood(
    observed_absolute_magnitude: np.ndarray,
    log_mass: np.ndarray,
    model_absolute_magnitude: np.ndarray,
    uniform: np.ndarray,
) -> np.ndarray:
    available = np.isfinite(observed_absolute_magnitude)
    observed = np.where(available, observed_absolute_magnitude, 0.0)
    residual = (
        observed[:, None, :] - model_absolute_magnitude[None, :, :]
    ) / MODEL_BAND_SCALE
    residual = np.where(available[:, None, :], residual, 0.0)
    distance = np.sum(residual**2, axis=2)
    distance -= np.min(distance, axis=1)[:, None]
    probability = np.exp(-0.5 * np.clip(distance, 0.0, 1400.0))
    # MASS_GRID is logarithmic, so equal cell weights are a log-uniform
    # (scale-invariant) mass prior rather than an IMF prior.
    probability /= probability.sum(axis=1)[:, None]
    cdf = np.cumsum(probability, axis=1)
    index = np.sum(cdf < uniform[:, None], axis=1)
    index = np.clip(index, 0, len(log_mass) - 1)
    result = log_mass[index].copy()
    if len(log_mass) > 1:
        step = float(np.median(np.diff(log_mass)))
        result += (uniform - 0.5) * step
        result = np.clip(result, log_mass[0], log_mass[-1])
    return np.exp(result)


def infer_mass_samples(
    row: pd.Series,
    av_probability: np.ndarray,
    family: str,
    rv: float,
    age_nodes: np.ndarray,
    n_draws: int = MASS_POSTERIOR_DRAWS,
    seed_offset: int = 0,
) -> np.ndarray:
    """Propagate A_V, photometry, age, and unresolved-binary uncertainty."""
    required = np.array([row.get(band, np.nan) for band in BANDS])
    errors = np.array([row.get(f"{band}_err", np.nan) for band in BANDS])
    if not np.isfinite(required[:3]).all():
        return np.full(n_draws, np.nan)
    available = np.isfinite(required) & np.isfinite(errors) & (errors > 0)
    effective_error = np.full(len(BANDS), np.nan)
    effective_error[available] = np.sqrt(
        errors[available] ** 2 + PHOTOMETRIC_FLOOR_MAG**2
    )
    source_seed = int(row["source_id"]) % 2_000_000_000
    rng = np.random.default_rng(RANDOM_SEED + source_seed + seed_offset)
    av = draw_from_grid(av_probability, rng, n_draws)
    if not np.isfinite(av).all():
        return np.full(n_draws, np.nan)
    noisy = np.full((len(BANDS), n_draws), np.nan)
    noisy[available] = (
        required[available, None]
        + effective_error[available, None]
        * rng.standard_normal((int(available.sum()), n_draws))
    )
    coefficient = np.array(
        [band_coefficients(rv)[band] for band in BANDS]
    )
    observed_absolute = (
        noisy.T - DIST_MODULUS - av[:, None] * coefficient[None, :]
    )

    age_index = np.arange(n_draws) % len(age_nodes)
    rng.shuffle(age_index)
    binary = rng.random(n_draws) < F_BINARY
    q_index = np.zeros(n_draws, dtype=int)
    q_index[binary] = rng.integers(1, N_Q_COMPONENTS + 1, size=int(binary.sum()))
    samples = np.full(n_draws, np.nan)
    for age_choice in np.unique(age_index):
        rounded_age = float(np.round(age_nodes[age_choice], 6))
        for q_choice in np.unique(q_index[age_index == age_choice]):
            selection = (age_index == age_choice) & (q_index == q_choice)
            log_mass, model_absolute = model_sequence(
                family, rounded_age, int(q_choice)
            )
            samples[selection] = _draw_mass_likelihood(
                observed_absolute[selection],
                log_mass,
                model_absolute,
                rng.random(int(selection.sum())),
            )
    return samples
