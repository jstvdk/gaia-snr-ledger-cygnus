#!/usr/bin/env python3
"""Step 3c, part 2: the repair_v4 WP5 fit with truth-side age marginalization.

Runs ``wp5_joint_age_fit.fit_joint`` over the full 54-branch grid (2 isochrone
families x 3 R_V x 3 IMF slopes x 3 subgroups), writes the repair_v4 WP5
products alongside the preserved earlier versions, and evaluates gate G3 of
the #1c plan.  The gate thresholds are the unmodified WP5 ones.

Everything except the forward response is identical to ``wp5_fit_imf.py``:
same observed counts, same Jeffreys k prior, same Dirichlet response
uncertainty, same association-mass and literature comparison.  The response is
now the WP4-posterior mixture over truth-age nodes, with node weights updated
by the Poisson likelihood of the observed counts (no gate statistic enters
them) and applied identically to all three subgroups.

Run:
  WP_REPAIR_VERSION=repair_v3 WP3_ANCHOR_PRIOR_MODE=variogram \
  PYTHONPATH=scripts python3 scripts/wp5_fit_imf_joint.py \
      --upstream-version repair_v3 --wp5-version repair_v4
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

import wp5_common as w
import wp5_fit_imf as F
import wp5_joint_age_fit as J


def load_nodes(
    subgroup: str,
    family: str,
    rv: float,
    prior: dict[float, float],
    version: str,
    stored_curves: pd.DataFrame,
    stored_responses: pd.DataFrame,
) -> tuple[dict[float, pd.DataFrame], dict[float, pd.DataFrame], dict[str, str]]:
    curves: dict[float, pd.DataFrame] = {}
    responses: dict[float, pd.DataFrame] = {}
    sources: dict[str, str] = {}
    for age in prior:
        node_response = J.node_response_path(subgroup, family, rv, age, version)
        node_curve = J.node_curve_path(subgroup, family, rv, age, version)
        reuse = (
            None
            if J.uses_age_interpolation(version)
            else J.reusable_scan_snapshot(subgroup, family, rv, age)
        )
        if node_response.exists() and node_curve.exists():
            responses[age] = pd.read_parquet(node_response)
            curves[age] = pd.read_parquet(node_curve)
            sources[f"{age:.3f}"] = str(node_response.relative_to(w.ROOT))
        elif reuse is not None:
            responses[age] = pd.read_parquet(reuse[0])
            curves[age] = pd.read_parquet(reuse[1])
            sources[f"{age:.3f}"] = str(reuse[0].relative_to(w.ROOT))
        else:
            raise RuntimeError(
                f"missing node response for {subgroup}/{family}/R_V={rv} at "
                f"{age:.3f} Myr — run scripts/wp5_injections_agenodes.py first"
            )
    return curves, responses, sources


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--upstream-version", default="repair_v3")
    parser.add_argument("--wp5-version", default="repair_v4")
    parser.add_argument(
        "--compare-version",
        default=None,
        help=(
            "previous WP5 normalization used for the A/C no-regression check; "
            "defaults to the upstream version when it has one on disk"
        ),
    )
    args = parser.parse_args()
    upstream = args.upstream_version
    version = args.wp5_version

    masses_input = w.PROC / f"wp4_mass_posteriors_{upstream}.parquet"
    masses = pd.read_parquet(masses_input)
    # Legacy single-age WP5 products exist only for upstream versions that had a
    # non-node injection run; a node-only version (repair_v5 onward) has none.
    curves_input = w.PROC / f"wp5_completeness_curves_{upstream}.parquet"
    response_input = w.PROC / f"wp5_injection_response_{upstream}.parquet"
    stored_curves = (
        pd.read_parquet(curves_input) if curves_input.exists() else None
    )
    stored_responses = (
        pd.read_parquet(response_input) if response_input.exists() else None
    )
    age_posterior = pd.read_parquet(w.PROC / f"wp4_age_posteriors_{upstream}.parquet")
    sample_store = np.load(w.PROC / f"wp4_mass_posterior_samples_{upstream}.npz")
    if not np.array_equal(
        sample_store["source_id"].astype("int64"),
        masses["source_id"].to_numpy("int64"),
    ):
        raise RuntimeError("mass posterior samples are not row-aligned")
    observed_samples = sample_store["samples"]
    if stored_responses is not None:
        probe_columns = stored_responses.columns
    else:
        probe = sorted(w.PROC.glob(f"wp5_agenode_*_{version}_response.parquet"))
        if not probe:
            raise RuntimeError(
                f"no node responses found for {version}; run "
                "scripts/wp5_injections_agenodes.py first"
            )
        probe_columns = pd.read_parquet(probe[0]).columns
    draw_columns = sorted(
        c for c in probe_columns if c.startswith("recovered_mass_draw_")
    )
    if not draw_columns:
        raise RuntimeError("response lacks mass-posterior draw columns")
    native = {family: J.native_isochrone_ages(family) for family in w.FAMILIES}
    interpolate = J.uses_age_interpolation(version)

    summaries = []
    bin_tables = []
    draw_store: dict[str, np.ndarray] = {}
    node_sources: dict[str, dict[str, str]] = {}
    age_posteriors_out = []
    mixture_curves = []
    for family in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            branch_samples = observed_samples[
                :,
                w.FAMILIES.index(family) * len(w.R_V_BRANCHES)
                + list(w.R_V_BRANCHES).index(rv),
                :,
            ]
            for subgroup in w.SUBGROUPS:
                prior = J.truth_age_nodes(
                    age_posterior, subgroup, family, rv, native[family],
                    snap=not interpolate,
                )
                curves, responses, sources = load_nodes(
                    subgroup, family, rv, prior, version,
                    stored_curves, stored_responses,
                )
                node_sources[f"{subgroup}|{family}|rv{rv:.1f}"] = sources
                # Canonical branch completeness curve: the PRIOR-weighted node
                # mixture.  It is deliberately not posterior-weighted -- the
                # calibration window must be defined before the fit, never from
                # the outcome (CUTS_AND_THRESHOLDS.md 6.4).  This is also the
                # curve fit_joint uses for the window edge.
                blended = J.mixture_curve(curves, prior)
                reference = curves[next(iter(prior))].sort_values("primary_mass")
                stage_columns = [
                    column
                    for column in [
                        "query_fraction",
                        "quality_fraction",
                        "cmd_ready_fraction",
                        "membership_fraction",
                        "recovery_fraction",
                    ]
                    if column in reference.columns
                ]
                stages = {
                    column: sum(
                        weight
                        * curves[age]
                        .sort_values("primary_mass")[column]
                        .to_numpy(float)
                        for age, weight in prior.items()
                    )
                    for column in stage_columns
                }
                mixture_curves.append(
                    pd.DataFrame(
                        {
                            "subgroup": subgroup,
                            "family": family,
                            "R_V": rv,
                            "primary_mass": blended["primary_mass"].to_numpy(float),
                            "n_truth_age_nodes": len(prior),
                            "recovery_isotonic": blended["recovery_isotonic"].to_numpy(float),
                            **stages,
                        }
                    )
                )
                for alpha in w.IMF_SLOPES:
                    summary, bins, k_draws, age_draws = J.fit_joint(
                        masses, curves, responses, prior, subgroup, family, rv,
                        alpha, np.random.default_rng(w.SEED), branch_samples,
                        draw_columns,
                    )
                    summaries.append(summary)
                    bin_tables.append(bins)
                    key = (
                        f"k__{subgroup.replace('-', '_')}__{family}"
                        f"__rv{rv:.1f}".replace(".", "p")
                        + f"__a{alpha:.1f}".replace(".", "p")
                    )
                    draw_store[key] = k_draws
                    draw_store[f"truth_age_draws__{key[3:]}"] = age_draws
                    age_posteriors_out.append(
                        {
                            "subgroup": subgroup,
                            "family": family,
                            "R_V": rv,
                            "alpha": alpha,
                            "nodes_Myr": summary["truth_age_nodes_Myr"],
                            "prior_weights": summary["truth_age_prior_weights"],
                            "posterior_weights": summary["truth_age_posterior_weights"],
                            "prior_mean_Myr": summary["truth_age_prior_mean_Myr"],
                            "posterior_mean_Myr": summary["truth_age_posterior_mean_Myr"],
                            "posterior_map_Myr": summary["truth_age_posterior_map_Myr"],
                        }
                    )

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
                kq = F.quantiles(k_total)
                pq = F.quantiles(k_total * primary_factor)
                sq = F.quantiles(k_total * (primary_factor + companion_factor))
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
                ] = k_total * (primary_factor + companion_factor)

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

    # Gate G3 of the #1c plan, on the unmodified WP5 gate thresholds.
    compare = args.compare_version or upstream
    previous_path = w.PROC / f"wp5_imf_normalization_{compare}.parquet"
    previous = pd.read_parquet(previous_path) if previous_path.exists() else None
    regressions = []
    for subgroup in ([] if previous is None else ["CygOB2-A", "CygOB2-C"]):
        for family in w.FAMILIES:
            for rv in w.R_V_BRANCHES:
                for alpha in w.IMF_SLOPES:
                    before = previous[
                        previous.subgroup.eq(subgroup) & previous.family.eq(family)
                        & previous.R_V.eq(rv) & previous.alpha.eq(alpha)
                    ].iloc[0]
                    after = normalization[
                        normalization.subgroup.eq(subgroup)
                        & normalization.family.eq(family)
                        & normalization.R_V.eq(rv) & normalization.alpha.eq(alpha)
                    ].iloc[0]
                    if bool(before.residual_gate_pass) and not bool(
                        after.residual_gate_pass
                    ):
                        regressions.append(
                            {
                                "subgroup": subgroup,
                                "family": family,
                                "R_V": float(rv),
                                "alpha": float(alpha),
                                "max_abs_residual_before": float(
                                    before.max_abs_pearson_residual
                                ),
                                "max_abs_residual_after": float(
                                    after.max_abs_pearson_residual
                                ),
                            }
                        )
    gate = {
        "baseline_all_subgroups_residual_gate": bool(
            baseline_residual["residual_gate_pass"].all()
        ),
        "no_A_or_C_branch_regression": not regressions,
        "regression_baseline_version": compare if previous is not None else None,
        "A_or_C_regressions": regressions,
        "baseline_mass_within_factor_two": bool(
            baseline["within_factor_two_literature"]
        ),
        "all_branch_subgroups_residual_gate": bool(
            normalization["residual_gate_pass"].all()
        ),
        "all_branch_masses_within_factor_two": bool(
            association["within_factor_two_literature"].all()
        ),
        "branches_passing": int(normalization["residual_gate_pass"].sum()),
        "branches_total": int(len(normalization)),
        "minimum_50_sources_per_subgroup_branch": bool(
            normalization["raw_calibration_sources"].ge(50).all()
        ),
        "G3_pass": bool(
            baseline_residual["residual_gate_pass"].all()
            and not regressions
            and baseline["within_factor_two_literature"]
        ),
    }

    norm_path = w.PROC / f"wp5_imf_normalization_{version}.parquet"
    bins_path = w.PROC / f"wp5_mass_function_bins_{version}.parquet"
    assoc_path = w.PROC / f"wp5_association_mass_{version}.parquet"
    draws_path = w.PROC / f"wp5_imf_posterior_draws_{version}.npz"
    curves_path = w.PROC / f"wp5_completeness_curves_{version}.parquet"
    normalization.to_parquet(norm_path, index=False)
    mass_bins.to_parquet(bins_path, index=False)
    association.to_parquet(assoc_path, index=False)
    pd.concat(mixture_curves, ignore_index=True).to_parquet(curves_path, index=False)
    np.savez_compressed(draws_path, **draw_store)

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_fit_imf_joint.py",
        "upstream_repair_version": upstream,
        "wp5_version": version,
        "status": (
            "WP5_REPAIR_BASELINE_ACCEPTED"
            if gate["G3_pass"]
            else "WP5_BLOCKED_RESIDUAL_OR_MASS_GATE"
        ),
        "issue": "#1c step 3 — adoption run with truth-side joint age–k fit",
        "seed": w.SEED,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
        "model_change_vs_repair_v3": (
            "the injection truth age is no longer a single upper-MS MAP: the "
            "forward response is the WP4-posterior mixture over truth-age "
            "nodes, with node weights updated by the Poisson likelihood of the "
            "observed counts (Jeffreys k integrated analytically).  Identical "
            "machinery for all three subgroups; no gate statistic enters the "
            "weights; no gate threshold moved.  Upstream WP3/WP4 unchanged."
        ),
        "truth_age_interpolation": interpolate,
        "node_rule": (
            (
                "wp4_repair_common.age_posterior_nodes unsnapped, weight 1/N "
                "each; the truth isochrone is interpolated between native ages "
                "(issue #13), matching the recovery side"
            )
            if interpolate
            else (
                "wp4_repair_common.age_posterior_nodes snapped to native "
                "isochrone ages, prior weight = summed node count"
            )
        ),
        "single_node_equivalence": (
            "with one node fit_joint reproduces wp5_fit_imf.fit_one bit-for-bit; "
            "verified in provenance/wp5_joint_fit_baseline_check_execution.json"
        ),
        "inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in [
                masses_input,
                w.PROC / f"wp4_age_posteriors_{upstream}.parquet",
                curves_input,
                response_input,
                previous_path,
            ]
            if path.exists()
        },
        "node_response_inputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in sorted(
                w.PROC.glob(f"wp5_agenode_*_{version}_response.parquet")
            )
        },
        "node_response_sources": node_sources,
        "truth_age_posteriors": age_posteriors_out,
        "gate": gate,
        "baseline": baseline.to_dict(),
        "completeness_curve_convention": (
            "wp5_completeness_curves_{version}.parquet holds the PRIOR-weighted "
            "truth-age node mixture per branch and subgroup.  It is not "
            "posterior-weighted: the calibration window must be fixed before "
            "the fit, never derived from the outcome (CUTS_AND_THRESHOLDS.md "
            "6.4).  Fitted node posteriors are in truth_age_posteriors above."
        ),
        "outputs": {
            str(path.relative_to(w.ROOT)): {
                "sha256": w.sha256(path),
                "bytes": path.stat().st_size,
            }
            for path in [norm_path, bins_path, assoc_path, curves_path, draws_path]
        },
    }
    w.write_json(w.PROVENANCE / f"wp5_imf_fit_execution_{version}.json", record)
    print(json.dumps({"gate": gate}, indent=2))
    for row in baseline_residual.itertuples():
        print(
            f"  {row.subgroup}: chi2p={row.poisson_chi_square_p:.4f} "
            f"trendp={row.residual_trend_p:.3f} "
            f"max|r|={row.max_abs_pearson_residual:.2f} "
            f"{'PASS' if row.residual_gate_pass else 'FAIL'}"
        )
    print(f"wrote provenance/wp5_imf_fit_execution_{version}.json")


if __name__ == "__main__":
    main()
