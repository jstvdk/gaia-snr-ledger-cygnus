#!/usr/bin/env python3
"""Issue #11: criteria R3a/R3b/R3c of CUTS_AND_THRESHOLDS.md section 14.6.

R3 as originally written demanded that two independent Monte-Carlo
realizations of an identical model always return the same gate verdict.  That
is unsatisfiable by any statistic: a cell with pass-probability pi disagrees
across independent realizations with probability 2*pi*(1-pi), so only a
degenerate test (pi in {0,1} everywhere) could comply.  Section 14.6 withdraws
it and replaces it with three criteria that test the defect that is actually
present.

R3a  no forbidden region: the widest gap between achievable p-values inside
     [0.01, 0.20], measured on null data simulated from every cell's own
     fitted lambda, must be <= 0.005.  The incumbent rank test has no
     achievable value between 0.0416 and 0.0724 -- a gap straddling the
     threshold, so a cell landing there has no verdict the data can support.

R3b  flip rate consistent with calibration: on the four identical-model pairs
     of issue #11 the observed flip count must lie inside the 95% interval of
     the Poisson-binomial implied by each cell's own pi.  Too few flips fails
     as well as too many.

R3c  indeterminacy declared: every cell carries pi, and cells with
     0.05 < pi < 0.95 are labelled indeterminate rather than counted as clean
     passes or failures.

pi is computed from the fit's own model of injection noise -- the Dirichlet
posterior over the response matrix's category counts, the same one already
propagated into every published k.  M response replicates are drawn, the
normalization is refitted on each by the production Jeffreys rule against the
UNCHANGED observed counts, and the full three-way gate is re-evaluated.

Outputs:
  tables/wp5_verdict_stability.csv
  provenance/wp5_verdict_stability_execution.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp5_verdict_stability.py [--versions ...]
"""
from __future__ import annotations

import argparse
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy import stats

import wp5_common as w
import wp5_joint_age_fit as J
import wp5_residual_trend as T

RESPONSE_REPLICATES = 400
REPLICATE_SEED = 20260729
NULL_SIMULATIONS = 20_000
NULL_SEED = 20260730
INDETERMINATE_LO = 0.05
INDETERMINATE_HI = 0.95
GAP_WINDOW = (0.01, 0.20)
GAP_LIMIT = 0.005

# The four identical-model cell pairs of issue #11: same subgroup, family, R_V
# and alpha, same single truth age, different injection realization.
IDENTICAL_MODEL_PAIRS = [
    ("CygOB2-A", "MIST", 3.1, 2.0),
    ("CygOB2-C", "MIST", 3.1, 2.0),
    ("CygOB2-C", "MIST", 3.1, 2.3),
    ("CygOB2-A", "MIST", 3.5, 2.6),
]

NODE_VERSIONS = {"repair_v4", "repair_v5", "repair_v6"}


def bins_path(version: str):
    return w.PROC / f"wp5_mass_function_bins_{version}.parquet"


def normalization_path(version: str):
    return w.PROC / f"wp5_imf_normalization_{version}.parquet"


class ResponseCache:
    """Category counts per (branch, alpha), built once and reused."""

    def __init__(self, version: str, age_posterior: pd.DataFrame, upstream: str):
        self.version = version
        self.age_posterior = age_posterior
        self.upstream = upstream
        self.native = {f: J.native_isochrone_ages(f) for f in w.FAMILIES}
        self.legacy = None
        if version not in NODE_VERSIONS:
            path = w.PROC / f"wp5_injection_response_{version}.parquet"
            if not path.exists():
                raise RuntimeError(f"no injection response on disk for {version}")
            self.legacy = pd.read_parquet(path)
        # Only the branch in hand is held.  The full node set is ~440 MB for a
        # 162-node version, which is enough to be killed when an injection run
        # is in flight at the same time.
        self._cache: dict[tuple, dict[float, pd.DataFrame]] = {}

    def _scan_snapshot(self, subgroup: str, family: str, rv: float, age: float):
        """repair_v4 reused the gate-G2 scan's CygOB2-B baseline nodes.

        ``J.reusable_scan_snapshot`` gates on the WP3_REPAIR_VERSION env var,
        which here names the version being *fitted*, not the upstream of the
        version being *scored*.  Resolve against this cache's own upstream.
        """
        if self.upstream != "repair_v3":
            return None
        if subgroup != "CygOB2-B" or family != "PARSEC" or not np.isclose(rv, 3.1):
            return None
        tag = f"{age:.3f}".replace(".", "p")
        response = w.PROC / f"wp5_age_scan_B_response_age{tag}_repair_v3.parquet"
        return response if response.exists() else None

    def responses(self, subgroup: str, family: str, rv: float) -> dict[float, pd.DataFrame]:
        key = (subgroup, family, rv)
        if key in self._cache:
            return self._cache[key]
        self._cache.clear()
        if self.legacy is not None:
            frame = self.legacy[
                self.legacy.subgroup.eq(subgroup)
                & self.legacy.family.eq(family)
                & self.legacy.R_V.eq(rv)
            ]
            out = {0.0: frame}
        else:
            prior = J.truth_age_nodes(
                self.age_posterior, subgroup, family, rv, self.native[family],
                snap=not J.uses_age_interpolation(self.version),
            )
            out = {}
            for age in prior:
                path = J.node_response_path(subgroup, family, rv, age, self.version)
                if not path.exists():
                    reuse = self._scan_snapshot(subgroup, family, rv, age)
                    if reuse is None:
                        raise RuntimeError(
                            f"missing node response {path.name} for {self.version}"
                        )
                    path = reuse
                out[age] = pd.read_parquet(path)
        self._cache[key] = out
        return out


def node_posterior_weights(
    normalization: pd.DataFrame, subgroup: str, family: str, rv: float, alpha: float,
    ages: list[float],
) -> np.ndarray:
    row = normalization[
        normalization.subgroup.eq(subgroup)
        & normalization.family.eq(family)
        & normalization.R_V.eq(rv)
        & normalization.alpha.eq(alpha)
    ]
    if len(row) != 1 or "truth_age_posterior_weights" not in normalization.columns:
        return np.full(len(ages), 1.0 / len(ages))
    weights = np.asarray(row["truth_age_posterior_weights"].iloc[0], dtype=float)
    if len(weights) != len(ages):
        return np.full(len(ages), 1.0 / len(ages))
    return weights / weights.sum()


def replicate_gate(
    category_counts: list[np.ndarray],
    true_masses: list[np.ndarray],
    node_weights: np.ndarray,
    alpha: float,
    weighted_counts: np.ndarray,
    log_mass: np.ndarray,
    null: np.ndarray,
    rng: np.random.Generator,
    replicates: int,
) -> dict:
    """Pass-probability pi under the response's own Dirichlet uncertainty."""
    weighted_n = float(weighted_counts.sum())
    rate = np.zeros((replicates, w.N_IMF_BINS), dtype=float)
    for counts, masses, weight in zip(category_counts, true_masses, node_weights):
        _, imf_weight = J._rate_bins(masses, counts, alpha)
        node_rate = np.zeros((replicates, w.N_IMF_BINS), dtype=float)
        for index in range(len(masses)):
            draw = rng.dirichlet(counts[index].astype(float) + 0.5, size=replicates)
            node_rate += draw[:, : w.N_IMF_BINS] * imf_weight[index]
        rate += weight * node_rate

    total_rate = rate.sum(axis=1)
    k = stats.gamma.ppf(0.5, a=weighted_n + 0.5, scale=1.0 / np.maximum(total_rate, 1e-300))
    expected = np.maximum(k[:, None] * rate, np.finfo(float).eps)
    residual = (weighted_counts[None, :] - expected) / np.sqrt(expected)

    chi_square = np.sum(residual**2, axis=1)
    chi_p = stats.chi2.sf(chi_square, max(w.N_IMF_BINS - 1, 1))
    max_abs = np.max(np.abs(residual), axis=1)
    centred = log_mass - log_mass.mean()
    statistic = (residual @ centred) / np.sqrt(float(np.sum(centred**2)))
    trend_p = np.array(
        [max(float(np.mean(np.abs(null) >= abs(value))), 1.0 / (len(null) + 1))
         for value in statistic]
    )
    incumbent_p = np.array(
        [T.incumbent_trend_p(log_mass, row) for row in residual]
    )
    passes = (chi_p >= 0.01) & (trend_p >= T.TREND_THRESHOLD) & (max_abs <= 3.0)
    passes_incumbent = (
        (chi_p >= 0.01) & (incumbent_p >= T.TREND_THRESHOLD) & (max_abs <= 3.0)
    )
    return {
        "pi": float(passes.mean()),
        "pi_incumbent": float(passes_incumbent.mean()),
        "replicate_trend_p_lo68": float(np.quantile(trend_p, 0.16)),
        "replicate_trend_p_hi68": float(np.quantile(trend_p, 0.84)),
        "replicate_max_abs_lo68": float(np.quantile(max_abs, 0.16)),
        "replicate_max_abs_hi68": float(np.quantile(max_abs, 0.84)),
    }


def score_version(version: str, replicates: int) -> pd.DataFrame:
    frame = pd.read_parquet(bins_path(version))
    normalization = pd.read_parquet(normalization_path(version))
    upstream = "repair_v5" if version in {"repair_v5", "repair_v6"} else "repair_v3"
    age_posterior = pd.read_parquet(
        w.PROC / f"wp4_age_posteriors_{upstream}.parquet"
    )
    cache = ResponseCache(version, age_posterior, upstream)
    rng = np.random.default_rng(REPLICATE_SEED)

    rows = []
    for (subgroup, family, rv), branch in frame.groupby(
        ["subgroup", "family", "R_V"]
    ):
        responses = cache.responses(subgroup, family, float(rv))
        ages = list(responses)
        draw_columns = sorted(
            c for c in responses[ages[0]].columns
            if c.startswith("recovered_mass_draw_")
        )
        # The calibration window is set by the PRIOR-weighted mixture curve, so
        # the bin edges — and hence the category counts — are identical across
        # the three alphas of a branch.  Build them once.
        category_cache: dict[bytes, tuple[list, list]] = {}
        for alpha, cell in branch.groupby("alpha"):
            ordered = cell.sort_values("bin_index")
            bin_edges = np.append(
                ordered["mass_lo"].to_numpy(float), ordered["mass_hi"].iloc[-1]
            )
            weighted_counts = ordered["membership_weighted_count"].to_numpy(float)
            residual = ordered["pearson_residual"].to_numpy(float)
            log_mass = np.log10(ordered["mass_geometric_center"].to_numpy(float))
            expected = ordered["expected_count_at_k_median"].to_numpy(float)
            rate = ordered[
                "completeness_weighted_imf_integral_per_k"
            ].to_numpy(float)

            null = T.bootstrap_null(expected, rate, log_mass)
            trend_p, statistic = T.replacement_trend_p(
                residual, expected, rate, log_mass, null=null
            )
            chi_square = float(np.sum(residual**2))
            chi_p = float(stats.chi2.sf(chi_square, max(w.N_IMF_BINS - 1, 1)))
            max_abs = float(np.max(np.abs(residual)))
            incumbent_p = T.incumbent_trend_p(log_mass, residual)

            key = bin_edges.tobytes()
            if key not in category_cache:
                counts, masses = [], []
                for age in ages:
                    true_mass, category = J._category_counts(
                        responses[age], bin_edges, draw_columns
                    )
                    counts.append(category)
                    masses.append(true_mass)
                category_cache[key] = (counts, masses)
            counts, masses = category_cache[key]
            weights = node_posterior_weights(
                normalization, subgroup, family, float(rv), float(alpha), ages
            )
            stability = replicate_gate(
                counts, masses, weights, float(alpha), weighted_counts,
                log_mass, null, rng, replicates,
            )
            rows.append(
                {
                    "version": version, "subgroup": subgroup, "family": family,
                    "R_V": float(rv), "alpha": float(alpha),
                    "chi2_p": chi_p, "max_abs_residual": max_abs,
                    "incumbent_trend_p": incumbent_p,
                    "replacement_trend_p": trend_p,
                    "replacement_T": statistic,
                    "gate_incumbent": T.gate_pass(chi_p, incumbent_p, max_abs),
                    "gate_replacement": T.gate_pass(chi_p, trend_p, max_abs),
                    **stability,
                    "indeterminate": bool(
                        INDETERMINATE_LO < stability["pi"] < INDETERMINATE_HI
                    ),
                    "indeterminate_incumbent": bool(
                        INDETERMINATE_LO < stability["pi_incumbent"] < INDETERMINATE_HI
                    ),
                }
            )
        print(f"  {version} {subgroup} {family} R_V={rv} done", flush=True)
    return pd.DataFrame(rows)


def measure_r3a(frame: pd.DataFrame) -> dict:
    """Widest gap between achievable p-values in [0.01, 0.20], on null data."""
    rng = np.random.default_rng(NULL_SEED)
    baseline = frame[
        frame.family.eq("PARSEC") & frame.R_V.eq(3.1) & frame.alpha.eq(2.3)
    ]
    incumbent_values, replacement_values = [], []
    for row in baseline.itertuples():
        cell = pd.read_parquet(bins_path(row.version))
        ordered = cell[
            cell.subgroup.eq(row.subgroup) & cell.family.eq(row.family)
            & cell.R_V.eq(row.R_V) & cell.alpha.eq(row.alpha)
        ].sort_values("bin_index")
        expected = ordered["expected_count_at_k_median"].to_numpy(float)
        rate = ordered["completeness_weighted_imf_integral_per_k"].to_numpy(float)
        log_mass = np.log10(ordered["mass_geometric_center"].to_numpy(float))
        null = T.bootstrap_null(expected, rate, log_mass)
        simulated = rng.poisson(np.maximum(expected, 1e-12), size=(NULL_SIMULATIONS, len(expected)))
        totals = simulated.sum(axis=1).astype(float)
        k = stats.gamma.ppf(0.5, a=totals + 0.5, scale=1.0 / float(np.sum(rate)))
        lam = np.maximum(k[:, None] * rate[None, :], 1e-12)
        residual = (simulated - lam) / np.sqrt(lam)
        centred = log_mass - log_mass.mean()
        statistic = (residual @ centred) / np.sqrt(float(np.sum(centred**2)))
        replacement_values.append(
            np.array([np.mean(np.abs(null) >= abs(v)) for v in statistic[:4000]])
        )
        incumbent_values.append(
            np.array([T.incumbent_trend_p(log_mass, r) for r in residual[:4000]])
        )

    def widest_gap(values: np.ndarray) -> dict:
        unique = np.unique(values)
        window = unique[(unique >= GAP_WINDOW[0]) & (unique <= GAP_WINDOW[1])]
        if len(window) < 2:
            return {"widest_gap": float(GAP_WINDOW[1] - GAP_WINDOW[0]),
                    "gap_containing_threshold": None, "n_achievable": int(len(window))}
        gaps = np.diff(window)
        below = window[window <= T.TREND_THRESHOLD]
        above = window[window > T.TREND_THRESHOLD]
        straddle = (
            float(above.min() - below.max())
            if len(below) and len(above)
            else None
        )
        return {
            "widest_gap": float(gaps.max()),
            "gap_containing_threshold": straddle,
            "gap_containing_threshold_bounds": (
                [float(below.max()), float(above.min())]
                if len(below) and len(above) else None
            ),
            "n_achievable": int(len(window)),
        }

    incumbent = widest_gap(np.concatenate(incumbent_values))
    replacement = widest_gap(np.concatenate(replacement_values))
    return {
        "criterion": (
            "widest gap between achievable p-values inside [0.01, 0.20] on null "
            f"data simulated from each cell's own fitted lambda, limit {GAP_LIMIT}"
        ),
        "null_simulations_per_cell": NULL_SIMULATIONS,
        "p_values_retained_per_cell": 4000,
        "cells": int(len(baseline)),
        "incumbent": incumbent,
        "replacement": replacement,
        "R3a_pass": bool(replacement["widest_gap"] <= GAP_LIMIT),
        "incumbent_would_pass": bool(incumbent["widest_gap"] <= GAP_LIMIT),
    }


def measure_r3b(frame: pd.DataFrame) -> dict:
    """Observed flips on the four identical-model pairs vs the implied rate."""
    observed_incumbent = {
        ("CygOB2-A", "MIST", 3.1, 2.0): True,
        ("CygOB2-C", "MIST", 3.1, 2.0): True,
        ("CygOB2-C", "MIST", 3.1, 2.3): True,
        ("CygOB2-A", "MIST", 3.5, 2.6): True,
    }
    # Under the replacement, 3 of the 4 pairs became consistent; only
    # CygOB2-C MIST R_V=3.1 alpha=2.3 still flips (recorded in
    # provenance/wp5_trend_replacement_validation_execution.json).
    observed_replacement = {
        ("CygOB2-A", "MIST", 3.1, 2.0): False,
        ("CygOB2-C", "MIST", 3.1, 2.0): False,
        ("CygOB2-C", "MIST", 3.1, 2.3): True,
        ("CygOB2-A", "MIST", 3.5, 2.6): False,
    }
    pairs = []
    for key in IDENTICAL_MODEL_PAIRS:
        subgroup, family, rv, alpha = key
        selection = frame[
            frame.subgroup.eq(subgroup) & frame.family.eq(family)
            & frame.R_V.eq(rv) & frame.alpha.eq(alpha)
            & frame.version.isin(["repair_v3", "repair_v4"])
        ]
        if len(selection) != 2:
            continue
        pi_replacement = float(selection["pi"].mean())
        pi_incumbent = float(selection["pi_incumbent"].mean())
        pairs.append(
            {
                "subgroup": subgroup, "family": family, "R_V": rv, "alpha": alpha,
                "pi_replacement_mean_over_realizations": round(pi_replacement, 4),
                "pi_incumbent_mean_over_realizations": round(pi_incumbent, 4),
                "implied_flip_probability_replacement": round(
                    2 * pi_replacement * (1 - pi_replacement), 4
                ),
                "implied_flip_probability_incumbent": round(
                    2 * pi_incumbent * (1 - pi_incumbent), 4
                ),
                "observed_flip_incumbent": observed_incumbent[key],
                "observed_flip_replacement": observed_replacement[key],
            }
        )

    def poisson_binomial_interval(probabilities: list[float]) -> tuple[int, int, float]:
        distribution = np.array([1.0])
        for probability in probabilities:
            distribution = np.convolve(distribution, [1 - probability, probability])
        cumulative = np.cumsum(distribution)
        lo = int(np.searchsorted(cumulative, 0.025))
        hi = int(np.searchsorted(cumulative, 0.975))
        return lo, hi, float(np.sum(np.arange(len(distribution)) * distribution))

    result = {
        "criterion": (
            "observed verdict flips on the four identical-model pairs must lie "
            "inside the 95% Poisson-binomial interval implied by each cell's pi; "
            "too few flips fails as well as too many"
        ),
        "pairs": pairs,
    }
    for label, key, observed in (
        ("replacement", "implied_flip_probability_replacement", observed_replacement),
        ("incumbent", "implied_flip_probability_incumbent", observed_incumbent),
    ):
        probabilities = [pair[key] for pair in pairs]
        lo, hi, mean = poisson_binomial_interval(probabilities)
        count = sum(observed[k] for k in IDENTICAL_MODEL_PAIRS if any(
            (p["subgroup"], p["family"], p["R_V"], p["alpha"]) == k for p in pairs
        ))
        result[label] = {
            "observed_flips": int(count),
            "expected_flips": round(mean, 3),
            "interval_95": [lo, hi],
            "consistent": bool(lo <= count <= hi),
        }
    if len(pairs) != len(IDENTICAL_MODEL_PAIRS):
        result["R3b_pass"] = None
        result["not_evaluated"] = (
            f"only {len(pairs)} of {len(IDENTICAL_MODEL_PAIRS)} identical-model "
            "pairs are available; R3b needs both repair_v3 and repair_v4 scored"
        )
        return result
    result["R3b_pass"] = bool(result["replacement"]["consistent"])
    return result


def classify_regressions(frame: pd.DataFrame, before: str, after: str) -> dict:
    """Gate-G3 no-regression clause under section 14.7.

    A CygOB2-A or CygOB2-C cell moving pass -> fail blocks acceptance under the
    strict per-branch reading.  Section 14.7 refines *what counts as a move*: a
    cell indeterminate in BOTH versions contributes to neither a regression nor
    an improvement, symmetrically, because its verdict is a coin flip on the
    injection realization.  Both readings are reported.
    """
    keys = ["subgroup", "family", "R_V", "alpha"]
    old = frame[frame.version.eq(before)].set_index(keys)
    new = frame[frame.version.eq(after)].set_index(keys)
    shared = old.index.intersection(new.index)
    strict, refined, improvements = [], [], []
    for key in shared:
        subgroup = key[0]
        if subgroup == "CygOB2-B":
            continue
        was, now = old.loc[key], new.loc[key]
        both_indeterminate = bool(was.indeterminate and now.indeterminate)
        entry = {
            "subgroup": subgroup, "family": key[1], "R_V": key[2], "alpha": key[3],
            f"{before}_trend_p": round(float(was.replacement_trend_p), 4),
            f"{after}_trend_p": round(float(now.replacement_trend_p), 4),
            f"{before}_max_abs_residual": round(float(was.max_abs_residual), 3),
            f"{after}_max_abs_residual": round(float(now.max_abs_residual), 3),
            f"{before}_pi": round(float(was.pi), 3),
            f"{after}_pi": round(float(now.pi), 3),
            "indeterminate_in_both": both_indeterminate,
        }
        if bool(was.gate_replacement) and not bool(now.gate_replacement):
            strict.append(entry)
            if not both_indeterminate:
                refined.append(entry)
        elif not bool(was.gate_replacement) and bool(now.gate_replacement):
            improvements.append({**entry, "counted": not both_indeterminate})
    return {
        "comparison": f"{before} -> {after}",
        "clause": "CygOB2-A and CygOB2-C must not regress (strict per-branch reading)",
        "strict_reading": {
            "regressions": strict,
            "count": len(strict),
            "passes": len(strict) == 0,
        },
        "refined_reading_section_14_7": {
            "regressions": refined,
            "count": len(refined),
            "passes": len(refined) == 0,
            "discounted_as_indeterminate_in_both": len(strict) - len(refined),
        },
        "improvements": improvements,
        "readings_agree": bool((len(strict) == 0) == (len(refined) == 0)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--versions", nargs="+",
        default=["repair_v3", "repair_v4", "repair_v5"],
    )
    parser.add_argument("--replicates", type=int, default=RESPONSE_REPLICATES)
    args = parser.parse_args()

    tables = []
    for version in args.versions:
        if not bins_path(version).exists():
            print(f"(skipping {version}: not on disk)")
            continue
        print(f"scoring {version} ...", flush=True)
        tables.append(score_version(version, args.replicates))
    frame = pd.concat(tables, ignore_index=True)
    out_csv = w.TABLES / "wp5_verdict_stability.csv"
    frame.to_csv(out_csv, index=False)

    r3a = measure_r3a(frame)
    r3b = measure_r3b(frame)

    per_version = []
    for version, block in frame.groupby("version", sort=False):
        indeterminate = int(block.indeterminate.sum())
        per_version.append(
            {
                "version": version,
                "cells": int(len(block)),
                "grid_replacement": int(block.gate_replacement.sum()),
                "grid_incumbent": int(block.gate_incumbent.sum()),
                "determinate_pass": int(
                    (block.gate_replacement & ~block.indeterminate).sum()
                ),
                "determinate_fail": int(
                    (~block.gate_replacement & ~block.indeterminate).sum()
                ),
                "indeterminate": indeterminate,
                "indeterminate_fraction": round(indeterminate / len(block), 4),
                "exceeds_25pc_indeterminate": bool(indeterminate / len(block) > 0.25),
                "median_pi": round(float(block["pi"].median()), 4),
                "indeterminate_incumbent": int(block.indeterminate_incumbent.sum()),
            }
        )

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_verdict_stability.py",
        "status": "SUCCESS",
        "issue": "#11 — criteria R3a/R3b/R3c",
        "predeclaration": "CUTS_AND_THRESHOLDS.md section 14.6, written before this run",
        "withdrawn_criterion": (
            "R3 as written in section 14.4 required identical verdicts on both "
            "realizations of an identical model.  A cell with pass-probability "
            "pi disagrees with probability 2*pi*(1-pi), so only a degenerate "
            "test could comply; the criterion was unsatisfiable by any "
            "statistic and is withdrawn, not weakened."
        ),
        "pi_definition": (
            "fraction of M response replicates, drawn from the fit's own "
            "Dirichlet posterior over the response matrix category counts, that "
            "pass the unchanged three-way gate.  The observed counts are held "
            "fixed; only the injection realization varies.  No new parameter "
            "and no new noise model."
        ),
        "response_replicates": args.replicates,
        "replicate_seed": REPLICATE_SEED,
        "null_seed": NULL_SEED,
        "indeterminate_band": [INDETERMINATE_LO, INDETERMINATE_HI],
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "inputs": {
            str(bins_path(v).relative_to(w.ROOT)): w.sha256(bins_path(v))
            for v in args.versions if bins_path(v).exists()
        },
        "R3a": r3a,
        "R3b": r3b,
        "R3c": {
            "criterion": (
                "every cell carries pi; cells with 0.05 < pi < 0.95 are "
                "labelled indeterminate and reported separately, never counted "
                "as a clean pass or a clean failure"
            ),
            "implemented": True,
            "per_version": per_version,
            # Section 14.6 requires only that pi be computed and indeterminate
            # cells be labelled and reported.  The 25% clause is section
            # 14.7(5), a *reporting obligation* on the grid, not an acceptance
            # criterion for the statistic -- do not conflate them.
            "R3c_pass": True,
            "section_14_7_5_reporting_obligation_triggered": bool(
                any(entry["exceeds_25pc_indeterminate"] for entry in per_version)
            ),
            "finding_if_triggered": (
                "more than 25% of the 54-cell grid has a verdict that is a coin "
                "flip under the injection experiment's own uncertainty.  The "
                "branch grid is underpowered at N_INJECT_PER_MASS = 400 and its "
                "per-cell verdicts must not be read as independent evidence.  "
                "This is a finding about the grid, not about the statistic."
            ),
        },
        "R3_overall_pass": bool(r3a["R3a_pass"] and bool(r3b.get("R3b_pass"))),
        "gate_G3_no_regression_clause": [
            classify_regressions(frame, before, after)
            for before, after in zip(args.versions[:-1], args.versions[1:])
        ],
        "outputs": {str(out_csv.relative_to(w.ROOT)): w.sha256(out_csv)},
    }
    w.write_json(w.PROVENANCE / "wp5_verdict_stability_execution.json", record)

    print("\nR3a  widest achievable-p gap in [0.01, 0.20]")
    print(f"  incumbent   {r3a['incumbent']['widest_gap']:.4f}"
          f"   straddling 0.05: {r3a['incumbent']['gap_containing_threshold']}")
    print(f"  replacement {r3a['replacement']['widest_gap']:.4f}"
          f"   straddling 0.05: {r3a['replacement']['gap_containing_threshold']}")
    print(f"  -> R3a {'PASS' if r3a['R3a_pass'] else 'FAIL'} (limit {GAP_LIMIT})")
    print("\nR3b  flips on the four identical-model pairs")
    for label in ("incumbent", "replacement"):
        block = r3b[label]
        print(f"  {label:12s} observed {block['observed_flips']}, expected "
              f"{block['expected_flips']:.2f}, 95% interval {block['interval_95']}"
              f"  -> {'consistent' if block['consistent'] else 'INCONSISTENT'}")
    print(f"  -> R3b {'PASS' if r3b['R3b_pass'] else 'FAIL'}")
    print("\nR3c  indeterminate cells")
    for entry in per_version:
        print(f"  {entry['version']:10s} grid {entry['grid_replacement']:2d}/54 = "
              f"{entry['determinate_pass']:2d} determinate pass + "
              f"{entry['indeterminate']:2d} indeterminate "
              f"({entry['indeterminate_fraction']:.0%})")
    if record["R3c"]["section_14_7_5_reporting_obligation_triggered"]:
        print("  -> section 14.7(5) triggered: >25% indeterminate, reported as a finding")
    print("\ngate G3 no-regression clause (A and C only)")
    for block in record["gate_G3_no_regression_clause"]:
        strict = block["strict_reading"]
        refined = block["refined_reading_section_14_7"]
        print(f"  {block['comparison']}: strict {strict['count']} regression(s) -> "
              f"{'PASS' if strict['passes'] else 'BLOCK'};  refined "
              f"{refined['count']} -> {'PASS' if refined['passes'] else 'BLOCK'} "
              f"({refined['discounted_as_indeterminate_in_both']} discounted)")
    print(f"\nR3 overall: {'PASS' if record['R3_overall_pass'] else 'FAIL'}")
    print("wrote provenance/wp5_verdict_stability_execution.json")


if __name__ == "__main__":
    main()
