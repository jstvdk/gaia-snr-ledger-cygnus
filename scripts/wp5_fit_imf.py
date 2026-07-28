#!/usr/bin/env python3
"""WP5 steps 2--5: completeness-corrected Poisson IMF normalization.

The plan's absolute 95% completeness edge is measured first.  If no mass in
2--8 Msun reaches it, the fit is *not* silently relabelled complete.  The
normalization is still estimable because the injection curve enters the
Poisson intensity explicitly.  In that case the primary fit uses the full
nominal 2--8 Msun window, is labelled ``corrected_no_absolute95_edge``, and the
unmet edge is carried as a documented WP5 deviation.

Each observed source contributes its WP2 membership probability.  For mass bin
i, the model is

    lambda_i = k * integral_i C(m) m**(-alpha) dm.

A Jeffreys prior p(k) proportional to k**(-1/2) gives a Gamma conditional
posterior.  Injection uncertainty is propagated by drawing each completeness
grid point from a Jeffreys-binomial Beta posterior around the monotone
injection curve, then drawing k conditional on the resulting integral.
"""
from __future__ import annotations

import argparse
import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import scipy
from scipy import stats

import wp5_common as w


def absolute_95_edge(curve: pd.DataFrame) -> float:
    use = curve[
        curve["primary_mass"].between(
            w.CALIBRATION_NOMINAL_LO, w.CALIBRATION_HI
        )
    ].sort_values("primary_mass")
    values = use["recovery_isotonic"].to_numpy(float)
    masses = use["primary_mass"].to_numpy(float)
    # The edge must remain >=95% at every sampled mass above it, not cross once
    # because of Monte-Carlo noise.
    suffix_minimum = np.minimum.accumulate(values[::-1])[::-1]
    eligible = np.flatnonzero(suffix_minimum >= w.COMPLETENESS_TARGET)
    return float(masses[eligible[0]]) if len(eligible) else float("nan")


def relative_plateau_edge(curve: pd.DataFrame) -> tuple[float, float]:
    """Diagnostic-only edge at 95% of the 6--8 Msun bright plateau."""
    use = curve[
        curve["primary_mass"].between(
            w.CALIBRATION_NOMINAL_LO, w.CALIBRATION_HI
        )
    ].sort_values("primary_mass")
    plateau = float(
        use.loc[use["primary_mass"].ge(6.0), "recovery_isotonic"].median()
    )
    target = w.COMPLETENESS_TARGET * plateau
    eligible = use[use["recovery_isotonic"].ge(target)]
    edge = float(eligible["primary_mass"].min()) if len(eligible) else float("nan")
    return edge, plateau


def interpolation_basis_integrals(
    mass_points: np.ndarray,
    bin_edges: np.ndarray,
    alpha: float,
    n_fine: int = 2001,
) -> np.ndarray:
    """Matrix A[bin,point] so rate_bin = A @ completeness_points."""
    matrix = np.zeros((len(bin_edges) - 1, len(mass_points)), dtype=float)
    for bin_index, (lo, hi) in enumerate(zip(bin_edges[:-1], bin_edges[1:], strict=True)):
        fine = np.geomspace(lo, hi, n_fine)
        imf = fine ** (-alpha)
        for point_index in range(len(mass_points)):
            basis = np.zeros(len(mass_points))
            basis[point_index] = 1.0
            completeness = np.interp(fine, mass_points, basis)
            matrix[bin_index, point_index] = np.trapezoid(
                completeness * imf, fine
            )
    return matrix


def quantiles(values: np.ndarray) -> dict[str, float]:
    q = np.quantile(values, [0.025, 0.16, 0.5, 0.84, 0.975])
    return {
        "lo95": float(q[0]),
        "lo68": float(q[1]),
        "median": float(q[2]),
        "hi68": float(q[3]),
        "hi95": float(q[4]),
    }


def fit_one(
    masses: pd.DataFrame,
    curve: pd.DataFrame,
    response: pd.DataFrame,
    subgroup: str,
    family: str,
    rv: float,
    alpha: float,
    rng: np.random.Generator,
    observed_mass_samples: np.ndarray | None = None,
    response_draw_columns: list[str] | None = None,
) -> tuple[dict, pd.DataFrame, np.ndarray]:
    mass_col = w.mass_column(family, rv)
    absolute_edge = absolute_95_edge(curve)
    plateau_edge, plateau = relative_plateau_edge(curve)
    if np.isfinite(absolute_edge):
        lower = max(w.CALIBRATION_NOMINAL_LO, absolute_edge)
        window_status = "absolute95_edge"
    else:
        # Explicit, predeclared fallback: corrected likelihood over the nominal
        # immortal interval.  This is not called 95% complete.
        lower = w.CALIBRATION_NOMINAL_LO
        window_status = "corrected_no_absolute95_edge"

    subgroup_mask = masses["subgroup"].eq(subgroup).to_numpy()
    point_window = masses[mass_col].between(lower, w.CALIBRATION_HI).to_numpy()
    sample = masses[subgroup_mask & point_window].copy()
    raw_n = int(len(sample))
    if raw_n < 50:
        raise RuntimeError(
            f"{subgroup}/{family}/R_V={rv} has only {raw_n} calibration sources"
        )

    bin_edges = np.geomspace(lower, w.CALIBRATION_HI, w.N_IMF_BINS + 1)
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
    weighted_n = float(weighted_counts.sum())
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

    # Forward response: the parent interval is the *injected* range, not the
    # observed window.  True masses on both sides of 2--8 Msun are measured into
    # it -- from below by binary brightening and mass-estimation scatter, from
    # above by the same scatter acting on the living members heavier than
    # 8 Msun.  Clipping the parent at CALIBRATION_HI silently deleted the second
    # channel and under-predicted the top observed bin; see
    # reports/WP5_RESIDUAL_DIAGNOSIS_CORRECTION_repair_v1.md.  Integrate the
    # whole injected parent interval and predict counts in the *observed* bins.
    response = response[
        response["true_primary_mass"].between(
            float(w.MASS_GRID.min()), float(w.MASS_GRID.max())
        )
    ].copy()
    true_mass = np.sort(response["true_primary_mass"].unique())
    n_categories = w.N_IMF_BINS + 1  # six observed bins plus lost/outside
    category_counts = np.zeros((len(true_mass), n_categories), dtype=float)
    for mass_index, true_value in enumerate(true_mass):
        injected_rows = response[response["true_primary_mass"].eq(true_value)]
        if response_draw_columns:
            active_draw_columns = [
                column
                for column in response_draw_columns
                if injected_rows[column].notna().any()
            ]
            if not active_draw_columns:
                category_counts[mass_index, -1] = float(len(injected_rows))
                continue
            recovered_draws = injected_rows[active_draw_columns].to_numpy(float)
            for bin_index, (lo, hi) in enumerate(
                zip(bin_edges[:-1], bin_edges[1:], strict=True)
            ):
                in_bin = (recovered_draws >= lo) & (
                    recovered_draws
                    < (hi if bin_index < w.N_IMF_BINS - 1 else hi + 1e-10)
                )
                category_counts[mass_index, bin_index] = float(
                    np.sum(in_bin) / len(active_draw_columns)
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

    response_probability = (
        category_counts[:, : w.N_IMF_BINS]
        / category_counts.sum(axis=1)[:, None]
    )
    trapezoid_weight = np.empty(len(true_mass))
    trapezoid_weight[0] = 0.5 * (true_mass[1] - true_mass[0])
    trapezoid_weight[-1] = 0.5 * (true_mass[-1] - true_mass[-2])
    trapezoid_weight[1:-1] = 0.5 * (true_mass[2:] - true_mass[:-2])
    imf_weight = trapezoid_weight * true_mass ** (-alpha)
    rate_bins = np.sum(
        response_probability * imf_weight[:, None], axis=0
    )
    rate_total = float(rate_bins.sum())
    k_mle = weighted_n / rate_total

    # Response uncertainty: a Jeffreys-multinomial Dirichlet posterior at each
    # true-mass grid point, including the lost/outside category.  This carries
    # both completeness and mass-bin migration uncertainty into k.
    rate_draws = np.zeros((w.N_POSTERIOR_DRAWS, w.N_IMF_BINS))
    for mass_index in range(len(true_mass)):
        probability_draw = rng.dirichlet(
            category_counts[mass_index].astype(float) + 0.5,
            size=w.N_POSTERIOR_DRAWS,
        )
        rate_draws += (
            probability_draw[:, : w.N_IMF_BINS]
            * imf_weight[mass_index]
        )
    total_rate_draws = rate_draws.sum(axis=1)
    shape = weighted_n + 0.5
    k_draws = rng.gamma(shape=shape, scale=1.0 / total_rate_draws)
    k_summary = quantiles(k_draws)

    expected = k_summary["median"] * rate_bins
    residual = (weighted_counts - expected) / np.sqrt(
        np.maximum(expected, np.finfo(float).eps)
    )
    chi_square = float(np.sum(residual**2))
    chi_df = max(w.N_IMF_BINS - 1, 1)
    chi_p = float(stats.chi2.sf(chi_square, chi_df))
    spearman = stats.spearmanr(
        np.log10(np.sqrt(bin_edges[:-1] * bin_edges[1:])),
        residual,
    )
    trend_p = float(spearman.pvalue) if np.isfinite(spearman.pvalue) else 1.0
    max_abs_residual = float(np.max(np.abs(residual)))
    residual_gate = bool(
        chi_p >= 0.01 and trend_p >= 0.05 and max_abs_residual <= 3.0
    )

    rows = []
    for index, (lo, hi) in enumerate(
        zip(bin_edges[:-1], bin_edges[1:], strict=True)
    ):
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
    bins = pd.DataFrame(rows)
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
        "response_true_mass_lower_Msun": float(true_mass.min()),
        "response_true_mass_upper_Msun": float(true_mass.max()),
        "raw_calibration_sources": raw_n,
        "membership_weighted_calibration_sources": weighted_n,
        "observed_count_method": count_method,
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
    return summary, bins, k_draws


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--repair-version",
        default=None,
        help="consume and write version-suffixed repair products (for example repair_v1)",
    )
    parser.add_argument(
        "--wp5-version",
        default=None,
        help=(
            "override the suffix for the WP5 injection inputs and every WP5 "
            "output, leaving --repair-version to select the upstream WP3/WP4 "
            "products.  Use when a WP5-only re-run must not overwrite an "
            "earlier repair version."
        ),
    )
    args = parser.parse_args()
    version = args.repair_version
    suffix = f"_{version}" if version else ""
    wp5_version = args.wp5_version or version
    wp5_suffix = f"_{wp5_version}" if wp5_version else ""
    masses_input = (
        w.PROC / f"wp4_mass_posteriors_{version}.parquet"
        if version
        else w.PROC / "wp4_masses.parquet"
    )
    curves_input = w.PROC / f"wp5_completeness_curves{wp5_suffix}.parquet"
    response_input = w.PROC / f"wp5_injection_response{wp5_suffix}.parquet"
    masses = pd.read_parquet(masses_input)
    curves = pd.read_parquet(curves_input)
    responses = pd.read_parquet(response_input)
    observed_samples_by_branch = None
    response_draw_columns = None
    if version:
        sample_store = np.load(
            w.PROC / f"wp4_mass_posterior_samples_{version}.npz"
        )
        if not np.array_equal(
            sample_store["source_id"].astype("int64"),
            masses["source_id"].to_numpy("int64"),
        ):
            raise RuntimeError("repair mass posterior samples are not row-aligned")
        observed_samples_by_branch = sample_store["samples"]
        response_draw_columns = sorted(
            column
            for column in responses.columns
            if column.startswith("recovered_mass_draw_")
        )
        if not response_draw_columns:
            raise RuntimeError("repair response lacks mass-posterior draw columns")
    rng = np.random.default_rng(w.SEED)

    summaries = []
    bin_tables = []
    draw_store: dict[str, np.ndarray] = {}
    for family in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            for alpha in w.IMF_SLOPES:
                for subgroup in w.SUBGROUPS:
                    curve = curves[
                        curves["family"].eq(family)
                        & curves["R_V"].eq(rv)
                        & curves["subgroup"].eq(subgroup)
                    ]
                    response = responses[
                        responses["family"].eq(family)
                        & responses["R_V"].eq(rv)
                        & responses["subgroup"].eq(subgroup)
                    ]
                    summary, bins, k_draws = fit_one(
                        masses,
                        curve,
                        response,
                        subgroup,
                        family,
                        rv,
                        alpha,
                        rng,
                        (
                            observed_samples_by_branch[
                                :,
                                w.FAMILIES.index(family) * len(w.R_V_BRANCHES)
                                + w.R_V_BRANCHES.index(rv),
                                :,
                            ]
                            if observed_samples_by_branch is not None
                            else None
                        ),
                        response_draw_columns,
                    )
                    summaries.append(summary)
                    bin_tables.append(bins)
                    key = (
                        f"k__{subgroup.replace('-', '_')}__{family}"
                        f"__rv{rv:.1f}".replace(".", "p")
                        + f"__a{alpha:.1f}".replace(".", "p")
                    )
                    draw_store[key] = k_draws

    normalization = pd.DataFrame(summaries)
    mass_bins = pd.concat(bin_tables, ignore_index=True)

    association_rows = []
    for family in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            for alpha in w.IMF_SLOPES:
                keys = [
                    (
                        f"k__{subgroup.replace('-', '_')}__{family}"
                        f"__rv{rv:.1f}".replace(".", "p")
                        + f"__a{alpha:.1f}".replace(".", "p")
                    )
                    for subgroup in w.SUBGROUPS
                ]
                k_total = sum(draw_store[key] for key in keys)
                primary_factor = w.primary_system_mass_per_k(alpha)
                companion_factor = w.companion_mass_per_k(alpha)
                primary_mass = k_total * primary_factor
                stellar_mass = k_total * (primary_factor + companion_factor)
                kq = quantiles(k_total)
                pq = quantiles(primary_mass)
                sq = quantiles(stellar_mass)
                within = (
                    w.LITERATURE_MASS_MSUN / w.LITERATURE_FACTOR_GATE
                    <= sq["median"]
                    <= w.LITERATURE_MASS_MSUN * w.LITERATURE_FACTOR_GATE
                )
                association_rows.append(
                    {
                        "family": family,
                        "R_V": rv,
                        "alpha": alpha,
                        "k_total_median": kq["median"],
                        "k_total_lo68": kq["lo68"],
                        "k_total_hi68": kq["hi68"],
                        "primary_system_mass_median_Msun": pq["median"],
                        "primary_system_mass_lo68_Msun": pq["lo68"],
                        "primary_system_mass_hi68_Msun": pq["hi68"],
                        "multiplicity_adjusted_mass_median_Msun": sq["median"],
                        "multiplicity_adjusted_mass_lo68_Msun": sq["lo68"],
                        "multiplicity_adjusted_mass_hi68_Msun": sq["hi68"],
                        "companion_mass_fraction": companion_factor
                        / (primary_factor + companion_factor),
                        "literature_mass_Msun": w.LITERATURE_MASS_MSUN,
                        "within_factor_two_literature": bool(within),
                    }
                )
                draw_store[
                    (
                        f"association_mass__{family}"
                        f"__rv{rv:.1f}".replace(".", "p")
                        + f"__a{alpha:.1f}".replace(".", "p")
                    )
                ] = stellar_mass

    association = pd.DataFrame(association_rows)
    baseline = association[
        association["family"].eq("PARSEC")
        & association["R_V"].eq(3.1)
        & association["alpha"].eq(2.3)
    ].iloc[0]
    baseline_residual = normalization[
        normalization["family"].eq("PARSEC")
        & normalization["R_V"].eq(3.1)
        & normalization["alpha"].eq(2.3)
    ]
    gate = {
        "absolute_95_edge_exists_all_branches": bool(
            normalization["absolute_95_edge_Msun"].notna().all()
        ),
        "minimum_50_sources_per_subgroup_branch": bool(
            normalization["raw_calibration_sources"].ge(50).all()
        ),
        "baseline_all_subgroups_residual_gate": bool(
            baseline_residual["residual_gate_pass"].all()
        ),
        "all_branch_subgroups_residual_gate": bool(
            normalization["residual_gate_pass"].all()
        ),
        "baseline_mass_within_factor_two": bool(
            baseline["within_factor_two_literature"]
        ),
        "all_branch_masses_within_factor_two": bool(
            association["within_factor_two_literature"].all()
        ),
    }

    norm_path = w.PROC / f"wp5_imf_normalization{wp5_suffix}.parquet"
    bins_path = w.PROC / f"wp5_mass_function_bins{wp5_suffix}.parquet"
    assoc_path = w.PROC / f"wp5_association_mass{wp5_suffix}.parquet"
    draws_path = w.PROC / f"wp5_imf_posterior_draws{wp5_suffix}.npz"
    normalization.to_parquet(norm_path, index=False)
    mass_bins.to_parquet(bins_path, index=False)
    association.to_parquet(assoc_path, index=False)
    np.savez_compressed(draws_path, **draw_store)

    inputs = [
        masses_input,
        curves_input,
        response_input,
    ]
    outputs = [norm_path, bins_path, assoc_path, draws_path]
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_fit_imf.py",
        "repair_version": version,
        "wp5_version": wp5_version,
        "status": (
            "WP5_REPAIR_BASELINE_ACCEPTED"
            if gate["baseline_all_subgroups_residual_gate"]
            and gate["baseline_mass_within_factor_two"]
            else "WP5_BLOCKED_RESIDUAL_OR_MASS_GATE"
        ),
        "seed": w.SEED,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
        "inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path) for path in inputs
        },
        "model": {
            "imf": "dN_primary_system/dm = k m^-alpha above 0.5 Msun",
            "slopes": w.IMF_SLOPES,
            "subgroups": w.SUBGROUPS,
            "families": w.FAMILIES,
            "R_V": w.R_V_BRANCHES,
            "nominal_window_Msun": [
                w.CALIBRATION_NOMINAL_LO,
                w.CALIBRATION_HI,
            ],
            "parent_mass_range_Msun": [
                float(w.MASS_GRID.min()),
                float(w.MASS_GRID.max()),
            ],
            "parent_range_note": (
                "The injected parent range is deliberately wider than the "
                "observed window on BOTH sides; the mass kernel migrates true "
                "masses across each edge.  These are two separate numbers with "
                "separate justifications and must not be tied together."
            ),
            "bins": w.N_IMF_BINS,
            "count": "sum of WP2 membership probabilities",
            "mass_posterior_counting": (
                "soft bin probabilities for repair versions; posterior medians "
                "for the frozen unversioned path"
            ),
            "likelihood": (
                "independent Poisson observed-mass bins with lambda_i=k "
                "integral R_i(m_true)m_true^-alpha dm_true; R includes "
                "selection and WP4 recovered-mass migration"
            ),
            "k_prior": "Jeffreys p(k) proportional to k^-1/2",
            "posterior_draws": w.N_POSTERIOR_DRAWS,
            "response_uncertainty": (
                "Jeffreys-multinomial Dirichlet draws over six observed-mass "
                "bins plus lost/outside, at every injected true mass"
            ),
            "total_mass_imf": "Kroupa-like slope 1.3 from 0.08--0.5 Msun, branch alpha from 0.5--120 Msun, continuous at 0.5",
            "multiplicity": f"primary-system k plus companions with f_bin={w.F_BINARY}, q~U[{w.Q_MIN},1]",
        },
        "calibration_window_result": {
            "absolute_target": w.COMPLETENESS_TARGET,
            "absolute_edge_found_count": int(
                normalization["absolute_95_edge_Msun"].notna().sum()
            ),
            "branch_subgroup_fit_count": int(len(normalization)),
            "fallback": "full 2--8 Msun corrected likelihood; explicitly not labelled 95% complete",
            "reason": "magnitude-independent WP2 quality and membership losses keep the bright plateau below 95%",
            "relative_plateau_edges_are_diagnostic_only": True,
        },
        "unassigned_policy": {
            "excluded_from_subgroup_normalizations": True,
            "reason": "no frozen subgroup footprint/completeness curve",
            "raw_2_8_baseline_sources": int(
                (
                    masses["subgroup"].eq("unassigned")
                    & masses["mass_PARSEC_rv3.1"].between(2.0, 8.0)
                ).sum()
            ),
            "membership_weighted_2_8_baseline_sources": float(
                masses.loc[
                    masses["subgroup"].eq("unassigned")
                    & masses["mass_PARSEC_rv3.1"].between(2.0, 8.0),
                    "membership_probability",
                ].sum()
            ),
        },
        "gate": gate,
        "baseline": baseline.to_dict(),
        "outputs": {
            str(path.relative_to(w.ROOT)): {
                "sha256": w.sha256(path),
                "bytes": path.stat().st_size,
            }
            for path in outputs
        },
    }
    provenance_name = (
        f"wp5_imf_fit_execution_{wp5_version}.json"
        if wp5_version
        else "wp5_imf_fit_execution.json"
    )
    w.write_json(w.PROVENANCE / provenance_name, record)
    print(json.dumps({"gate": gate, "baseline": baseline.to_dict()}, indent=2))


if __name__ == "__main__":
    main()
