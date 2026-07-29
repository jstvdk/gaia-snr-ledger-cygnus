#!/usr/bin/env python3
"""Joint age--k WP5 fit: truth-side age marginalized under the WP4 posterior.

Step 3 machinery for issue #1c.  The gate-G2 scan
(provenance/wp5_age_scan_execution.json) showed that WP5's injection truth
model generates photometry at a single age, that CygOB2-B's residual is an
artifact of that single age sitting on the PMS/Henyey isochrone fold, and that
plain marginalization is insufficient because B's WP4 age posterior is
bottom-heavy.  The brief's section 4 authorizes one alternative: a joint
age--k fit with the WP4 posterior as prior, applied identically to every
subgroup.  This module implements it.

Node rule (no new parameter anywhere).  The truth-side age nodes are the
*recovery* side's own posterior discretization --
``wp4_repair_common.age_posterior_nodes``, nine equiprobable split-normal
nodes -- snapped to the native isochrone ages the truth generator already
snaps to, with prior weights equal to the summed node counts.  The truth and
recovery sides therefore discretize the same WP4 posterior the same way.

Fit.  For each node j the forward rates r_ij per unit k are built exactly as
in ``wp5_fit_imf.fit_one``.  With the Jeffreys prior p(k) ~ k^(-1/2) the
normalization integrates out analytically, giving the node marginal
likelihood

    log ML_j = sum_i n_i log r_ij - (N + 1/2) log R_j + const,

with n_i the membership-weighted observed counts (identical across nodes, so
the constant cancels).  Posterior node weights are w_j ~ prior_j * ML_j; the
reported rates are the mixture sum_j w_j r_ij, and the k posterior draws a
node per draw before its Dirichlet response draw, so response and age
uncertainty both propagate into k.

Anti-tuning compliance: no gate statistic enters the weights -- only the
Poisson likelihood of the observed counts under the WP4 prior -- and the same
machinery runs for all three subgroups.  With a single node the whole
computation reduces to ``fit_one`` exactly, RNG stream included; that is
asserted by ``check_single_node_equivalence``.

Library module.  ``scripts/wp5_joint_fit_baseline_check.py`` is the step-3b
runner; the step-3c adoption chain drives it over the full grid.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import wp5_common as w
import wp5_fit_imf as F
from wp4_repair_common import age_posterior_nodes


# Versions whose truth model interpolates the isochrone between native ages
# and therefore does not snap the posterior nodes (issue #13).  Declared once
# here so the injection runner and the fit runner cannot disagree about which
# node set a version owns -- they address the same files by age.
AGE_INTERPOLATED_VERSIONS = frozenset({"repair_v6", "repair_v7"})


def uses_age_interpolation(version: str) -> bool:
    return version in AGE_INTERPOLATED_VERSIONS


def native_isochrone_ages(family: str) -> np.ndarray:
    frame = pd.read_parquet(w.PROC / f"wp3_isochrones_{family.lower()}.parquet")
    return np.sort(frame["age_Myr"].unique())


def snap_to_native(age: float, native_ages: np.ndarray) -> float:
    return float(native_ages[int(np.argmin(np.abs(native_ages - age)))])


def truth_age_nodes(
    age_posterior: pd.DataFrame,
    subgroup: str,
    family: str,
    rv: float,
    native_ages: np.ndarray | None = None,
    snap: bool = True,
) -> dict[float, float]:
    """Truth age -> WP4 prior weight, uniform rule for every branch.

    ``snap=True`` is the repair_v4 rule: the nine posterior quantiles are
    snapped to native isochrone ages and their weights summed.  Because the
    native grid is coarse this is a *discontinuous* function of the posterior —
    measured at 109x amplification in
    ``provenance/wp5_node_rule_continuity_execution.json`` — so an arbitrarily
    small WP4 shift can delete a node carrying tens of percent of the weight.

    ``snap=False`` (issue #13, repair_v6) keeps the quantiles where the
    posterior puts them, each with weight 1/9, and leaves it to
    ``wp5_common.load_isochrone_between_ages`` to build the truth isochrone
    there.  Node positions are then affine in the posterior and the rule is
    Lipschitz-1: exactly as continuous as a rule can be.  No new parameter
    either way — the node count is still ``wp4_repair_common.N_AGE_NODES``.
    """
    nodes = age_posterior_nodes(age_posterior, subgroup, family, rv)
    if not snap:
        return {float(node): 1.0 / len(nodes) for node in np.sort(nodes)}
    if native_ages is None:
        native_ages = native_isochrone_ages(family)
    weights: dict[float, float] = {}
    for node in nodes:
        native = snap_to_native(float(node), native_ages)
        weights[native] = weights.get(native, 0.0) + 1.0 / len(nodes)
    return dict(sorted(weights.items()))


def mixture_curve(
    curves: dict[float, pd.DataFrame], weights: dict[float, float]
) -> pd.DataFrame:
    """Posterior-weighted completeness curve, on the shared primary-mass grid."""
    reference = curves[next(iter(weights))].sort_values("primary_mass")
    mass = reference["primary_mass"].to_numpy(float)
    blended = np.zeros(len(mass), dtype=float)
    for age, weight in weights.items():
        node = curves[age].sort_values("primary_mass")
        if not np.allclose(node["primary_mass"].to_numpy(float), mass):
            raise RuntimeError("node completeness curves are not on a common mass grid")
        blended += weight * node["recovery_isotonic"].to_numpy(float)
    return pd.DataFrame({"primary_mass": mass, "recovery_isotonic": blended})


def _observed_counts(
    masses: pd.DataFrame,
    subgroup: str,
    mass_col: str,
    bin_edges: np.ndarray,
    observed_mass_samples: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, int, str]:
    """The observed-count block of fit_one, verbatim in behaviour."""
    subgroup_mask = masses["subgroup"].eq(subgroup).to_numpy()
    point_window = masses[mass_col].between(bin_edges[0], w.CALIBRATION_HI).to_numpy()
    sample = masses[subgroup_mask & point_window].copy()
    raw_n = int(len(sample))
    if raw_n < 50:
        raise RuntimeError(f"{subgroup} has only {raw_n} calibration sources")
    if observed_mass_samples is None:
        weighted_counts = np.array(
            [
                sample.loc[
                    sample[mass_col].ge(lo)
                    & (
                        sample[mass_col].lt(hi)
                        if index < w.N_IMF_BINS - 1
                        else sample[mass_col].le(hi)
                    ),
                    "membership_probability",
                ].sum()
                for index, (lo, hi) in enumerate(
                    zip(bin_edges[:-1], bin_edges[1:], strict=True)
                )
            ],
            dtype=float,
        )
        count_method = "posterior_median_hard_bins"
    else:
        subgroup_samples = observed_mass_samples[subgroup_mask]
        subgroup_weights = masses.loc[
            subgroup_mask, "membership_probability"
        ].to_numpy(float)
        weighted_counts = np.array(
            [
                np.sum(
                    subgroup_weights
                    * np.mean(
                        (subgroup_samples >= lo)
                        & (
                            subgroup_samples
                            < (hi if index < w.N_IMF_BINS - 1 else hi + 1e-10)
                        ),
                        axis=1,
                    )
                )
                for index, (lo, hi) in enumerate(
                    zip(bin_edges[:-1], bin_edges[1:], strict=True)
                )
            ],
            dtype=float,
        )
        count_method = "membership_weighted_mass_posterior_soft_bins"
    raw_counts = np.array(
        [
            int(
                (
                    sample[mass_col].ge(lo)
                    & (
                        sample[mass_col].lt(hi)
                        if index < w.N_IMF_BINS - 1
                        else sample[mass_col].le(hi)
                    )
                ).sum()
            )
            for index, (lo, hi) in enumerate(
                zip(bin_edges[:-1], bin_edges[1:], strict=True)
            )
        ],
        dtype=int,
    )
    return weighted_counts, raw_counts, raw_n, count_method


def _category_counts(
    response: pd.DataFrame,
    bin_edges: np.ndarray,
    response_draw_columns: list[str] | None,
) -> tuple[np.ndarray, np.ndarray]:
    """The forward-response block of fit_one, verbatim in behaviour."""
    response = response[
        response["true_primary_mass"].between(
            float(w.MASS_GRID.min()), float(w.MASS_GRID.max())
        )
    ].copy()
    true_mass = np.sort(response["true_primary_mass"].unique())
    n_categories = w.N_IMF_BINS + 1
    category_counts = np.zeros((len(true_mass), n_categories), dtype=float)
    for mass_index, true_value in enumerate(true_mass):
        injected_rows = response[response["true_primary_mass"].eq(true_value)]
        if response_draw_columns:
            active = [
                column
                for column in response_draw_columns
                if injected_rows[column].notna().any()
            ]
            if not active:
                category_counts[mass_index, -1] = float(len(injected_rows))
                continue
            recovered_draws = injected_rows[active].to_numpy(float)
            for bin_index, (lo, hi) in enumerate(
                zip(bin_edges[:-1], bin_edges[1:], strict=True)
            ):
                in_bin = (recovered_draws >= lo) & (
                    recovered_draws
                    < (hi if bin_index < w.N_IMF_BINS - 1 else hi + 1e-10)
                )
                category_counts[mass_index, bin_index] = float(
                    np.sum(in_bin) / len(active)
                )
            category_counts[mass_index, -1] = (
                len(injected_rows) - category_counts[mass_index, :-1].sum()
            )
        else:
            recovered_mass = injected_rows["recovered_mass"].to_numpy(float)
            assigned = np.zeros(len(injected_rows), dtype=bool)
            for bin_index, (lo, hi) in enumerate(
                zip(bin_edges[:-1], bin_edges[1:], strict=True)
            ):
                in_bin = (recovered_mass >= lo) & (
                    recovered_mass
                    < (hi if bin_index < w.N_IMF_BINS - 1 else hi + 1e-10)
                )
                category_counts[mass_index, bin_index] = float(in_bin.sum())
                assigned |= in_bin
            category_counts[mass_index, -1] = float((~assigned).sum())
        if not np.isclose(category_counts[mass_index].sum(), len(injected_rows)):
            raise RuntimeError("injection response category accounting failed")
    return true_mass, category_counts


def _rate_bins(
    true_mass: np.ndarray, category_counts: np.ndarray, alpha: float
) -> tuple[np.ndarray, np.ndarray]:
    response_probability = (
        category_counts[:, : w.N_IMF_BINS] / category_counts.sum(axis=1)[:, None]
    )
    trapezoid_weight = np.empty(len(true_mass))
    trapezoid_weight[0] = 0.5 * (true_mass[1] - true_mass[0])
    trapezoid_weight[-1] = 0.5 * (true_mass[-1] - true_mass[-2])
    trapezoid_weight[1:-1] = 0.5 * (true_mass[2:] - true_mass[:-2])
    imf_weight = trapezoid_weight * true_mass ** (-alpha)
    return np.sum(response_probability * imf_weight[:, None], axis=0), imf_weight


def fit_joint(
    masses: pd.DataFrame,
    curves: dict[float, pd.DataFrame],
    responses: dict[float, pd.DataFrame],
    prior_weights: dict[float, float],
    subgroup: str,
    family: str,
    rv: float,
    alpha: float,
    rng: np.random.Generator,
    observed_mass_samples: np.ndarray | None = None,
    response_draw_columns: list[str] | None = None,
) -> tuple[dict, pd.DataFrame]:
    """Joint age--k fit for one subgroup on one branch."""
    ages = list(prior_weights)
    mass_col = w.mass_column(family, rv)
    blended = mixture_curve(curves, prior_weights)
    absolute_edge = F.absolute_95_edge(blended)
    plateau_edge, plateau = F.relative_plateau_edge(blended)
    if np.isfinite(absolute_edge):
        lower = max(w.CALIBRATION_NOMINAL_LO, absolute_edge)
        window_status = "absolute95_edge"
    else:
        lower = w.CALIBRATION_NOMINAL_LO
        window_status = "corrected_no_absolute95_edge"
    bin_edges = np.geomspace(lower, w.CALIBRATION_HI, w.N_IMF_BINS + 1)

    weighted_counts, raw_counts, raw_n, count_method = _observed_counts(
        masses, subgroup, mass_col, bin_edges, observed_mass_samples
    )
    weighted_n = float(weighted_counts.sum())

    node_true_mass = {}
    node_category_counts = {}
    node_rate_bins = {}
    log_marginal = {}
    for age in ages:
        true_mass, category_counts = _category_counts(
            responses[age], bin_edges, response_draw_columns
        )
        rate_bins, _ = _rate_bins(true_mass, category_counts, alpha)
        node_true_mass[age] = true_mass
        node_category_counts[age] = category_counts
        node_rate_bins[age] = rate_bins
        log_marginal[age] = float(
            np.sum(weighted_counts * np.log(rate_bins))
            - (weighted_n + 0.5) * np.log(rate_bins.sum())
        )

    log_posterior = {
        age: np.log(prior_weights[age]) + log_marginal[age] for age in ages
    }
    peak = max(log_posterior.values())
    unnormalized = {age: float(np.exp(value - peak)) for age, value in log_posterior.items()}
    total = sum(unnormalized.values())
    posterior_weights = {age: value / total for age, value in unnormalized.items()}

    rate_bins = np.zeros(w.N_IMF_BINS, dtype=float)
    for age in ages:
        rate_bins += posterior_weights[age] * node_rate_bins[age]
    rate_total = float(rate_bins.sum())
    k_mle = weighted_n / rate_total

    # Response uncertainty per node, then age selection, then k.  With one node
    # the selection draw is skipped, so the RNG stream matches fit_one exactly.
    node_rate_draws = {}
    for age in ages:
        draws = np.zeros((w.N_POSTERIOR_DRAWS, w.N_IMF_BINS))
        _, imf_weight = _rate_bins(node_true_mass[age], node_category_counts[age], alpha)
        for mass_index in range(len(node_true_mass[age])):
            probability_draw = rng.dirichlet(
                node_category_counts[age][mass_index].astype(float) + 0.5,
                size=w.N_POSTERIOR_DRAWS,
            )
            draws += probability_draw[:, : w.N_IMF_BINS] * imf_weight[mass_index]
        node_rate_draws[age] = draws
    if len(ages) == 1:
        rate_draws = node_rate_draws[ages[0]]
        age_draws = np.full(w.N_POSTERIOR_DRAWS, ages[0], dtype=float)
    else:
        selection = rng.choice(
            len(ages),
            size=w.N_POSTERIOR_DRAWS,
            p=np.array([posterior_weights[age] for age in ages]),
        )
        stacked = np.stack([node_rate_draws[age] for age in ages])
        rate_draws = stacked[selection, np.arange(w.N_POSTERIOR_DRAWS)]
        age_draws = np.array(ages, dtype=float)[selection]
    total_rate_draws = rate_draws.sum(axis=1)
    shape = weighted_n + 0.5
    k_draws = rng.gamma(shape=shape, scale=1.0 / total_rate_draws)
    k_summary = F.quantiles(k_draws)

    expected = k_summary["median"] * rate_bins
    residual = (weighted_counts - expected) / np.sqrt(
        np.maximum(expected, np.finfo(float).eps)
    )
    chi_square = float(np.sum(residual**2))
    chi_df = max(w.N_IMF_BINS - 1, 1)
    chi_p = float(stats.chi2.sf(chi_square, chi_df))
    spearman = stats.spearmanr(
        np.log10(np.sqrt(bin_edges[:-1] * bin_edges[1:])), residual
    )
    trend_p = float(spearman.pvalue) if np.isfinite(spearman.pvalue) else 1.0
    max_abs_residual = float(np.max(np.abs(residual)))
    residual_gate = bool(chi_p >= 0.01 and trend_p >= 0.05 and max_abs_residual <= 3.0)

    rows = []
    for index, (lo, hi) in enumerate(zip(bin_edges[:-1], bin_edges[1:], strict=True)):
        rows.append(
            {
                "subgroup": subgroup,
                "family": family,
                "R_V": rv,
                "alpha": alpha,
                "bin_index": index,
                "mass_lo": lo,
                "mass_hi": hi,
                "mass_geometric_center": np.sqrt(lo * hi),
                "raw_count": int(raw_counts[index]),
                "membership_weighted_count": weighted_counts[index],
                "completeness_weighted_imf_integral_per_k": rate_bins[index],
                "expected_count_at_k_median": expected[index],
                "pearson_residual": residual[index],
            }
        )
    summary = {
        "subgroup": subgroup,
        "family": family,
        "R_V": rv,
        "alpha": alpha,
        "calibration_lower_Msun": lower,
        "calibration_upper_Msun": w.CALIBRATION_HI,
        "window_status": window_status,
        "absolute_95_edge_Msun": absolute_edge,
        "relative_95_plateau_edge_Msun_diagnostic": plateau_edge,
        "bright_plateau_completeness": plateau,
        "raw_calibration_sources": raw_n,
        "membership_weighted_calibration_sources": weighted_n,
        "observed_count_method": count_method,
        "truth_age_nodes_Myr": [float(age) for age in ages],
        "truth_age_prior_weights": [float(prior_weights[age]) for age in ages],
        "truth_age_log_marginal_likelihood": [float(log_marginal[age]) for age in ages],
        "truth_age_posterior_weights": [float(posterior_weights[age]) for age in ages],
        "truth_age_posterior_mean_Myr": float(
            sum(age * posterior_weights[age] for age in ages)
        ),
        "truth_age_prior_mean_Myr": float(
            sum(age * prior_weights[age] for age in ages)
        ),
        "truth_age_posterior_map_Myr": float(
            max(posterior_weights, key=posterior_weights.get)
        ),
        "k_mle": k_mle,
        "k_median": k_summary["median"],
        "k_lo68": k_summary["lo68"],
        "k_hi68": k_summary["hi68"],
        "k_lo95": k_summary["lo95"],
        "k_hi95": k_summary["hi95"],
        "poisson_chi_square": chi_square,
        "poisson_chi_square_df": chi_df,
        "poisson_chi_square_p": chi_p,
        "residual_spearman_rho": float(spearman.statistic),
        "residual_trend_p": trend_p,
        "max_abs_pearson_residual": max_abs_residual,
        "residual_gate_pass": residual_gate,
    }
    return summary, pd.DataFrame(rows), k_draws, age_draws


def check_single_node_equivalence(
    masses: pd.DataFrame,
    curve: pd.DataFrame,
    response: pd.DataFrame,
    subgroup: str,
    family: str,
    rv: float,
    alpha: float,
    observed_mass_samples: np.ndarray | None,
    response_draw_columns: list[str] | None,
) -> dict:
    """A one-node joint fit must reproduce fit_one bit-for-bit."""
    reference_summary, reference_bins, reference_k = F.fit_one(
        masses,
        curve,
        response,
        subgroup,
        family,
        rv,
        alpha,
        np.random.default_rng(w.SEED),
        observed_mass_samples,
        response_draw_columns,
    )
    age = float(curve["age_fit_Myr"].iloc[0]) if "age_fit_Myr" in curve else 0.0
    joint_summary, joint_bins, joint_k, _ = fit_joint(
        masses,
        {age: curve},
        {age: response},
        {age: 1.0},
        subgroup,
        family,
        rv,
        alpha,
        np.random.default_rng(w.SEED),
        observed_mass_samples,
        response_draw_columns,
    )
    shared = [
        "calibration_lower_Msun",
        "membership_weighted_calibration_sources",
        "k_median",
        "poisson_chi_square",
        "poisson_chi_square_p",
        "residual_trend_p",
        "max_abs_pearson_residual",
        "residual_gate_pass",
    ]
    return {
        "subgroup": subgroup,
        "family": family,
        "R_V": rv,
        "alpha": alpha,
        "identical_summary_fields": {
            key: bool(reference_summary[key] == joint_summary[key]) for key in shared
        },
        "max_abs_residual_difference": float(
            np.max(
                np.abs(
                    reference_bins.sort_values("bin_index")["pearson_residual"].to_numpy()
                    - joint_bins.sort_values("bin_index")["pearson_residual"].to_numpy()
                )
            )
        ),
        "max_abs_k_draw_difference": float(np.max(np.abs(reference_k - joint_k))),
        "equivalent": bool(
            all(reference_summary[key] == joint_summary[key] for key in shared)
            and np.array_equal(reference_k, joint_k)
        ),
    }


def node_stem(subgroup: str, family: str, rv: float, age: float, version: str) -> str:
    label = subgroup.replace("CygOB2-", "")
    rv_tag = f"{rv:.1f}".replace(".", "p")
    age_tag = f"{age:.3f}".replace(".", "p")
    return f"wp5_agenode_{label}_{family}_rv{rv_tag}_age{age_tag}_{version}"


def node_response_path(
    subgroup: str, family: str, rv: float, age: float, version: str
) -> Path:
    return w.PROC / f"{node_stem(subgroup, family, rv, age, version)}_response.parquet"


def node_curve_path(
    subgroup: str, family: str, rv: float, age: float, version: str
) -> Path:
    return w.PROC / f"{node_stem(subgroup, family, rv, age, version)}_curve.parquet"


def reusable_scan_snapshot(
    subgroup: str, family: str, rv: float, age: float
) -> tuple[Path, Path] | None:
    """The gate-G2 scan already generated CygOB2-B's baseline nodes with the
    same fresh-seed recipe, so they are reused rather than recomputed."""
    # The stored snapshots were generated under the repair_v3 extinction.  Any
    # other upstream version has different photometry, so reuse would silently
    # mix extinction models.
    from wp3_repair_common import REPAIR_VERSION as _upstream

    if _upstream != "repair_v3":
        return None
    if subgroup != "CygOB2-B" or family != "PARSEC" or not np.isclose(rv, 3.1):
        return None
    tag = f"{age:.3f}".replace(".", "p")
    response = w.PROC / f"wp5_age_scan_B_response_age{tag}_repair_v3.parquet"
    curve = w.PROC / f"wp5_age_scan_B_curve_age{tag}_repair_v3.parquet"
    return (response, curve) if response.exists() and curve.exists() else None
