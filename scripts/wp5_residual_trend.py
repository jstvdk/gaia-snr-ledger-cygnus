#!/usr/bin/env python3
"""The WP5 residual-trend statistic: incumbent rank test and its replacement.

Pre-declared in CUTS_AND_THRESHOLDS.md section 14 (2026-07-28) as the fix for
open issue #11 -- the incumbent Spearman rank test is not stable against
Monte-Carlo resampling of an unchanged model, because it responds to the
residuals' rank order rather than their magnitude and because its p-value is
confined to a lattice with no achievable value between 0.042 and 0.072.

Both statistics are implemented here so every comparison in the validation is
made on identical inputs.

Replacement statistic (section 14.2):

    b = sum (x_i - xbar) res_i / sum (x_i - xbar)^2
    T = b * sqrt(sum (x_i - xbar)^2)

with x = log10(bin geometric centre) and res the Pearson residuals, which are
already standardized so no variance need be estimated from 6 points.

Null (section 14.3): parametric bootstrap.  Simulate Poisson counts from the
fitted expectations, refit the normalization by the same Jeffreys rule used in
production, recompute the residuals and hence T, and read the p-value off that
distribution.  Exact for the actual bin counts.
"""
from __future__ import annotations

import numpy as np
from scipy import stats

BOOTSTRAP_DRAWS = 20_000
BOOTSTRAP_SEED = 20260728
TREND_THRESHOLD = 0.05


def slope_statistic(log_mass: np.ndarray, residual: np.ndarray) -> float:
    """T of section 14.2: standardized-residual least-squares slope."""
    centred = log_mass - log_mass.mean()
    sxx = float(np.sum(centred**2))
    if sxx <= 0:
        return 0.0
    slope = float(np.sum(centred * residual) / sxx)
    return slope * np.sqrt(sxx)


def incumbent_trend_p(log_mass: np.ndarray, residual: np.ndarray) -> float:
    """The Spearman rank test currently in wp5_fit_imf.fit_one."""
    result = stats.spearmanr(log_mass, residual)
    return float(result.pvalue) if np.isfinite(result.pvalue) else 1.0


def jeffreys_k(total_counts: float, total_rate: float) -> float:
    """Posterior median of Gamma(N + 1/2, rate=R), the production k rule."""
    if total_rate <= 0:
        return 0.0
    return float(stats.gamma.ppf(0.5, a=total_counts + 0.5, scale=1.0 / total_rate))


def bootstrap_null(
    expected: np.ndarray,
    rate_per_k: np.ndarray,
    log_mass: np.ndarray,
    draws: int = BOOTSTRAP_DRAWS,
    seed: int = BOOTSTRAP_SEED,
) -> np.ndarray:
    """Null distribution of T under Poisson counts about ``expected``.

    The normalization is refitted on every simulated dataset, so the constraint
    that k was estimated from the same counts is carried exactly rather than
    approximated.
    """
    rng = np.random.default_rng(seed)
    simulated = rng.poisson(np.maximum(expected, 1e-12), size=(draws, len(expected)))
    total_rate = float(np.sum(rate_per_k))
    totals = simulated.sum(axis=1).astype(float)
    # Vectorized form of [jeffreys_k(v, total_rate) for v in totals]; scipy's
    # ppf is elementwise so this is numerically identical, only faster.
    k = (
        stats.gamma.ppf(0.5, a=totals + 0.5, scale=1.0 / total_rate)
        if total_rate > 0
        else np.zeros_like(totals)
    )
    lam = k[:, None] * rate_per_k[None, :]
    lam = np.maximum(lam, 1e-12)
    residual = (simulated - lam) / np.sqrt(lam)
    centred = log_mass - log_mass.mean()
    sxx = float(np.sum(centred**2))
    return (residual @ centred) / np.sqrt(sxx)


def replacement_trend_p(
    residual: np.ndarray,
    expected: np.ndarray,
    rate_per_k: np.ndarray,
    log_mass: np.ndarray,
    draws: int = BOOTSTRAP_DRAWS,
    seed: int = BOOTSTRAP_SEED,
    null: np.ndarray | None = None,
) -> tuple[float, float]:
    """Return (p_value, T_observed) for the replacement statistic."""
    observed = slope_statistic(log_mass, residual)
    if null is None:
        null = bootstrap_null(expected, rate_per_k, log_mass, draws, seed)
    p_value = float(np.mean(np.abs(null) >= abs(observed)))
    # never report an impossible zero from a finite bootstrap
    p_value = max(p_value, 1.0 / (len(null) + 1))
    return p_value, observed


def gate_pass(chi2_p: float, trend_p: float, max_abs_residual: float) -> bool:
    """The unchanged three-way conjunction."""
    return bool(chi2_p >= 0.01 and trend_p >= TREND_THRESHOLD and max_abs_residual <= 3.0)
