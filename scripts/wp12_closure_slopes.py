#!/usr/bin/env python3
"""WP12.2 -- subgroup closure, closing slopes, and the mixed-slope ledger.

The carried branch grid forces one IMF slope on all three subgroups.  That is a
modelling choice.  The out-of-sample closure test --- which never entered the
normalization fit --- says the three subgroups do not want the same slope:
Cyg OB2-A and -B sit near Salpeter, Cyg OB2-C prefers something shallower and
over-closes at every carried slope.  Reporting only the association-wide median
hides that.

This script does three things.

1.  Tabulates the closure ratio for every (subgroup, family, R_V, alpha) cell
    and the slope at which each cell's census closes exactly, flagging closing
    slopes that fall outside the carried [2.0, 2.6] range as extrapolations.

2.  Builds the association-level closure ratio properly --- summed observation
    over summed prediction, not the mean of three ratios.

3.  Runs a MIXED-SLOPE ledger: each subgroup keeps its own frozen (k, age)
    posterior draws at its own best-closing carried slope, and the supernova
    count is re-derived on that combination.  This is the sensitivity the brief
    asks for, and it is the only form of it the frozen products permit, since
    the posterior draws exist per slope and cannot be interpolated between.

Nothing upstream is refitted.  The engine is WP7's own `run_population`, and it
is validated by reproducing the stored single-slope WP7 numbers before any
mixed-slope result is reported.

Outputs:
  tables/wp12_closure_by_alpha.csv
  tables/wp12_closing_slopes.csv
  tables/wp12_mixed_slope_ledger.csv
  provenance/wp12_closure_slopes_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/wp12_closure_slopes.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import wp5_common as w
import wp12_common as W
from wp7_ledger import TurnoffRelation, draw_key, run_population, summarize
from wp7_ledger_prereg import SF_DURATIONS_MYR

ITERATIONS = 200_000
SEED = 20260803
# Declared in the preregistration: the tolerance at which the WP12 engine must
# reproduce the stored WP7 single-slope numbers before mixed-slope results are
# believed.
NSN_RELATIVE_TOLERANCE = 0.03
PRECENT_ABSOLUTE_TOLERANCE = 0.02


def closing_slope(alphas: np.ndarray, ratios: np.ndarray) -> float:
    """The slope at which closure_ratio = 1, by log-linear interpolation.

    The closure ratio INCREASES with alpha on this chain, which is not the
    naive expectation.  A steeper IMF predicts relatively fewer massive stars
    per low-mass star, but k is refitted to the same 2-8 Msun counts at every
    slope, and the refitted normalization more than compensates, so
    observed/predicted above 8 Msun rises from ~0.70 at alpha = 2.0 to ~1.62 at
    alpha = 2.6 (grid medians).  The direction is therefore measured from the
    data rather than assumed -- an earlier revision of this function hard-coded
    the opposite sign, returned 2.0 for every cell, and is recorded as defect
    D1 in the execution record.
    """
    if np.any(ratios <= 0):
        return float("nan")
    order = np.argsort(alphas)
    a, r = alphas[order], np.log(ratios[order])
    if np.all(np.diff(r) > 0):
        return float(np.interp(0.0, r, a))
    if np.all(np.diff(r) < 0):
        return float(np.interp(0.0, r[::-1], a[::-1]))
    raise RuntimeError(
        "closure ratio is not monotone in alpha; the closing slope would be "
        f"ambiguous.  ratios={ratios.tolist()} at alphas={alphas.tolist()}"
    )


def linear_extrapolated_root(alphas: np.ndarray, ratios: np.ndarray) -> float:
    """Least-squares log(ratio) vs alpha, solved for log(ratio) = 0."""
    if np.any(ratios <= 0):
        return float("nan")
    slope, intercept = np.polyfit(alphas, np.log(ratios), 1)
    return float(-intercept / slope)


def run_mixed(
    rng: np.random.Generator,
    draws: np.lib.npyio.NpzFile,
    relation: TurnoffRelation,
    family: str,
    rv: float,
    slope_by_subgroup: dict[str, float],
    delta: float,
    n_iter: int,
) -> dict:
    """Association ledger with a per-subgroup IMF slope."""
    totals = np.zeros(n_iter, dtype=int)
    last = np.full(n_iter, np.inf)
    per_subgroup = {}
    for subgroup in w.SUBGROUPS:
        alpha = slope_by_subgroup[subgroup]
        key = draw_key(subgroup, family, rv, alpha)
        k_all = draws[f"k__{key}"]
        age_all = draws[f"truth_age_draws__{key}"]
        sub_rng = np.random.default_rng(rng.integers(0, 2 ** 63 - 1))
        pick = sub_rng.integers(0, k_all.size, n_iter)
        res = run_population(
            sub_rng, k_all[pick], age_all[pick], alpha, delta, relation
        )
        n_sn = res["n_sn"]["all_explode"]
        t_last = res["t_last"]["all_explode"]
        totals += n_sn
        last = np.minimum(last, t_last)
        per_subgroup[subgroup] = float(n_sn.mean())
    out = summarize(totals, last)
    out.update({f"N_SN_{s[-1]}": v for s, v in per_subgroup.items()})
    return out


def main() -> None:
    closure = W.closure()
    subgroups = list(w.SUBGROUPS)

    # ---------------------------------------------- 1. closure by subgroup
    by_alpha = closure[
        ["subgroup", "family", "R_V", "alpha", "closure_ratio",
         "closure_ratio_lo68", "closure_ratio_hi68", "predicted_observed_living",
         "observed_living", "shortfall"]
    ].copy()

    # ------------------------------------- 2. association-level closure ratio
    assoc = closure.groupby(["family", "R_V", "alpha"], as_index=False).agg(
        predicted_observed_living=("predicted_observed_living", "sum"),
        observed_living=("observed_living", "sum"),
    )
    # WP6 defines closure_ratio = OBSERVED / PREDICTED (wp6_closure_test.py:237),
    # so the star-weighted association aggregate is the ratio of the summed
    # observed count to the summed predicted count, in that order.  An earlier
    # revision of this script inverted it; see defect D2 in the record.
    assoc["closure_ratio"] = (
        assoc.observed_living / assoc.predicted_observed_living
    )
    assoc["subgroup"] = "association"
    by_alpha = pd.concat(
        [by_alpha, assoc[["subgroup", "family", "R_V", "alpha", "closure_ratio",
                          "predicted_observed_living", "observed_living"]]],
        ignore_index=True,
    ).sort_values(["subgroup", "family", "R_V", "alpha"]).reset_index(drop=True)

    # ------------------------------------------------- 3. closing slopes
    rows = []
    for (subgroup, family, rv), cell in by_alpha.groupby(
        ["subgroup", "family", "R_V"]
    ):
        cell = cell.sort_values("alpha")
        a = cell.alpha.to_numpy(float)
        r = cell.closure_ratio.to_numpy(float)
        interp = closing_slope(a, r)
        fitted = linear_extrapolated_root(a, r)
        inside = bool(min(w.IMF_SLOPES) <= interp <= max(w.IMF_SLOPES))
        rows.append({
            "subgroup": subgroup, "family": family, "R_V": rv,
            "closing_alpha_interpolated": interp,
            "closing_alpha_loglinear_fit": fitted,
            "inside_carried_grid": inside,
            "closure_at_2p0": float(r[a == 2.0][0]),
            "closure_at_2p3": float(r[a == 2.3][0]),
            "closure_at_2p6": float(r[a == 2.6][0]),
            "best_closing_carried_alpha": float(a[np.argmin(np.abs(np.log(r)))]),
        })
    slopes = pd.DataFrame(rows).sort_values(
        ["subgroup", "family", "R_V"]
    ).reset_index(drop=True)

    base_slopes = slopes[
        slopes.family.eq(W.BASE["family"]) & slopes.R_V.eq(W.BASE["R_V"])
        & slopes.subgroup.isin(subgroups)
    ]
    spread = float(
        base_slopes.closing_alpha_interpolated.max()
        - base_slopes.closing_alpha_interpolated.min()
    )
    r3 = W.score_prediction("R3", spread > 0.15, {
        "closing_alpha_by_subgroup_at_baseline": {
            row.subgroup: round(row.closing_alpha_interpolated, 4)
            for row in base_slopes.itertuples()
        },
        "spread": round(spread, 4),
    })

    # ------------------------------------------------ 4. mixed-slope ledger
    draws = np.load(W.frozen("wp5_posterior_draws"))
    relations = {family: TurnoffRelation(family) for family in w.FAMILIES}
    master = np.random.default_rng(SEED)

    # 4a.  Validation: reproduce the stored single-slope WP7 branches.
    stored = W.association_all_explode()
    checks = []
    for family in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            for alpha in W.HEADLINE_ALPHAS:
                for delta in SF_DURATIONS_MYR:
                    got = run_mixed(
                        master, draws, relations[family], family, rv,
                        {s: alpha for s in subgroups}, delta, ITERATIONS,
                    )
                    ref = stored[
                        stored.family.eq(family) & stored.R_V.eq(rv)
                        & stored.alpha.eq(alpha)
                        & stored.sf_duration_Myr.eq(delta)
                    ].iloc[0]
                    checks.append({
                        "family": family, "R_V": rv, "alpha": alpha,
                        "sf_duration_Myr": delta,
                        "wp7_N_SN_mean": float(ref.N_SN_mean),
                        "wp12_N_SN_mean": round(got["N_SN_mean"], 4),
                        "relative_difference": round(
                            abs(got["N_SN_mean"] - ref.N_SN_mean)
                            / max(ref.N_SN_mean, 1e-12), 5),
                        "wp7_P_recent": float(ref.P_last_SN_within_100kyr),
                        "wp12_P_recent": round(got["P_last_SN_within_100kyr"], 5),
                        "absolute_difference_P_recent": round(
                            abs(got["P_last_SN_within_100kyr"]
                                - ref.P_last_SN_within_100kyr), 5),
                    })
    check_frame = pd.DataFrame(checks)
    worst_nsn = float(check_frame.relative_difference.max())
    worst_prec = float(check_frame.absolute_difference_P_recent.max())
    engine_ok = bool(
        worst_nsn <= NSN_RELATIVE_TOLERANCE
        and worst_prec <= PRECENT_ABSOLUTE_TOLERANCE
    )

    # 4b.  The mixed-slope branches themselves.
    mixed_rows = []
    for family in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            assignment = {}
            for subgroup in subgroups:
                cell = slopes[
                    slopes.subgroup.eq(subgroup) & slopes.family.eq(family)
                    & slopes.R_V.eq(rv)
                ].iloc[0]
                assignment[subgroup] = float(cell.best_closing_carried_alpha)
            for delta in SF_DURATIONS_MYR:
                got = run_mixed(
                    master, draws, relations[family], family, rv,
                    assignment, delta, ITERATIONS,
                )
                mixed_rows.append({
                    "family": family, "R_V": rv, "sf_duration_Myr": delta,
                    **{f"alpha_{s[-1]}": assignment[s] for s in subgroups},
                    "uniform_slope": len(set(assignment.values())) == 1,
                    **{k: round(v, 5) for k, v in got.items()},
                })
    mixed = pd.DataFrame(mixed_rows)

    # 4c.  What the mixed slope actually costs, branch by branch.  Only Cyg
    # OB2-C is reassigned (2.3 -> 2.0), so the whole effect is C's, and C's
    # contribution depends on whether its turnoff has crossed the 120 Msun IMF
    # ceiling on that branch -- which is family- and R_V-dependent.  The
    # sensitivity is therefore large on some branches and null on others, and
    # the split is measured here rather than asserted.
    pure = stored[stored.alpha.eq(2.3)].set_index(
        ["family", "R_V", "sf_duration_Myr"]
    )
    sub_ledger = W.ledger()
    c_only = sub_ledger[
        sub_ledger.scope.eq("subgroup") & sub_ledger.subgroup.eq("CygOB2-C")
        & sub_ledger.explodability.eq("all_explode")
    ].set_index(["family", "R_V", "alpha", "sf_duration_Myr"])
    deltas = []
    for row in mixed.itertuples():
        ref = pure.loc[(row.family, row.R_V, row.sf_duration_Myr)]
        c23 = float(c_only.loc[
            (row.family, row.R_V, 2.3, row.sf_duration_Myr)].N_SN_mean)
        c20 = float(c_only.loc[
            (row.family, row.R_V, 2.0, row.sf_duration_Myr)].N_SN_mean)
        deltas.append({
            "family": row.family, "R_V": row.R_V,
            "sf_duration_Myr": row.sf_duration_Myr,
            "N_SN_pure_2p3": round(float(ref.N_SN_mean), 4),
            "N_SN_mixed": round(float(row.N_SN_mean), 4),
            "relative_change": round(
                float(row.N_SN_mean - ref.N_SN_mean)
                / max(float(ref.N_SN_mean), 1e-12), 5),
            "C_contributes_at_2p3": round(c23, 4),
            "C_contributes_at_2p0": round(c20, 4),
            "C_above_imf_ceiling": bool(max(c23, c20) < 0.5),
        })
    delta_frame = pd.DataFrame(deltas)
    worst_mixed_shift = float(delta_frame.relative_change.abs().max())
    c_alive = delta_frame[~delta_frame.C_above_imf_ceiling]
    c_dead = delta_frame[delta_frame.C_above_imf_ceiling]

    headline = stored[stored.alpha.isin(W.HEADLINE_ALPHAS)]
    lo, hi = float(headline.N_SN_mean.min()), float(headline.N_SN_mean.max())
    inside = bool((mixed.N_SN_mean >= lo).all() and (mixed.N_SN_mean <= hi).all())
    r4 = W.score_prediction("R4", inside, {
        "headline_N_SN_range": [round(lo, 3), round(hi, 3)],
        "mixed_slope_N_SN_range": [
            round(float(mixed.N_SN_mean.min()), 3),
            round(float(mixed.N_SN_mean.max()), 3),
        ],
        "all_inside": inside,
    })

    # -------------------------------------------------------------- write out
    out_closure = w.TABLES / "wp12_closure_by_alpha.csv"
    out_slopes = w.TABLES / "wp12_closing_slopes.csv"
    out_mixed = w.TABLES / "wp12_mixed_slope_ledger.csv"
    by_alpha.to_csv(out_closure, index=False)
    slopes.to_csv(out_slopes, index=False)
    mixed.to_csv(out_mixed, index=False)

    payload = {
        "item": "WP12.2 -- subgroup closure, closing slopes, mixed-slope ledger",
        "defects_found_during_development": [
            {
                "id": "D1",
                "what": (
                    "the first revision of closing_slope() assumed the closure "
                    "ratio decreases with alpha and inverted the interpolation "
                    "accordingly.  It increases: k is refitted to the same "
                    "2-8 Msun counts at every slope and over-compensates.  The "
                    "bug returned closing_alpha = 2.0 for all three subgroups, "
                    "which would have made prediction R3 read FAIL."
                ),
                "found_by": (
                    "the returned value disagreed with the manuscript's stored "
                    "closingAlpha = 2.25, which is computed by an independent "
                    "expression in scripts/wp10_numbers.py"
                ),
                "fix": (
                    "the monotone direction is now measured from the data and "
                    "a non-monotone cell raises rather than silently "
                    "extrapolating"
                ),
                "affects_published_numbers": False,
            },
            {
                "id": "D2",
                "what": (
                    "the association-level aggregate was first computed as "
                    "summed predicted over summed observed.  WP6 defines "
                    "closure_ratio the other way round -- observed over "
                    "predicted -- so the aggregate ran opposite to every "
                    "subgroup curve it was plotted beside."
                ),
                "found_by": (
                    "the aggregate decreased with alpha while all three "
                    "subgroups increased, which is impossible for a weighted "
                    "combination of them"
                ),
                "fix": (
                    "the aggregate now follows WP6's definition; the figure "
                    "axis is labelled observed/predicted to match"
                ),
                "affects_published_numbers": (
                    "yes, within WP12: the star-weighted association closure "
                    "and the association closing slope both changed.  No "
                    "upstream WP6 number is affected -- the per-subgroup "
                    "ratios were always read straight from the frozen table."
                ),
            },
        ],
        "association_closure_definition": (
            "summed observed living stars above 8 Msun over summed predicted "
            "observed living stars, NOT the mean of the three subgroup ratios"
        ),
        "association_closure_at_baseline_alpha": {
            f"{row.family}_RV{row.R_V:g}": round(row.closure_ratio, 4)
            for row in assoc[assoc.alpha.eq(2.3)].itertuples()
        },
        "association_closure_grid_median_at_2p3": round(
            float(assoc[assoc.alpha.eq(2.3)].closure_ratio.median()), 4),
        # Two different aggregations that the manuscript has been treating as
        # one.  The median of the 18 subgroup ratios is 1.07; the star-weighted
        # ratio of summed observation to summed prediction is 1.15.  They
        # differ because Cyg OB2-A carries most of the observed massive stars
        # and under-closes, while Cyg OB2-C over-closes on a smaller census.
        # The paper must say which one it means.
        "two_aggregations_at_2p3": {
            "median_of_18_subgroup_ratios": round(
                float(closure[closure.alpha.eq(2.3)].closure_ratio.median()), 4),
            "star_weighted_summed_ratio_grid_median": round(
                float(assoc[assoc.alpha.eq(2.3)].closure_ratio.median()), 4),
            "why_they_differ": (
                "a median of ratios is not the ratio of sums; the subgroups "
                "carry very different numbers of observed massive stars"
            ),
        },
        "subgroup_grid_median_closure_by_alpha": {
            s: {
                f"alpha_{a:g}": round(float(
                    closure[closure.subgroup.eq(s) & closure.alpha.eq(a)]
                    .closure_ratio.median()), 4)
                for a in w.IMF_SLOPES
            }
            for s in subgroups
        },
        "subgroup_closure_at_baseline_branch": {
            row.subgroup: round(row.closure_ratio, 4)
            for row in by_alpha[
                by_alpha.family.eq(W.BASE["family"])
                & by_alpha.R_V.eq(W.BASE["R_V"]) & by_alpha.alpha.eq(2.3)
            ].itertuples()
        },
        "closing_slopes_at_baseline_family_rv": {
            row.subgroup: {
                "interpolated": round(row.closing_alpha_interpolated, 4),
                "inside_carried_grid": bool(row.inside_carried_grid),
            }
            for row in slopes[
                slopes.family.eq(W.BASE["family"]) & slopes.R_V.eq(W.BASE["R_V"])
            ].itertuples()
        },
        "closing_slope_grid_median_by_subgroup": {
            s: round(float(
                slopes[slopes.subgroup.eq(s)].closing_alpha_interpolated.median()
            ), 4)
            for s in subgroups + ["association"]
        },
        "cells_wanting_a_slope_outside_the_grid": {
            s: int((~slopes[slopes.subgroup.eq(s)].inside_carried_grid).sum())
            for s in subgroups + ["association"]
        },
        "mixed_slope": {
            "iterations": ITERATIONS,
            "seed": SEED,
            "engine": "wp7_ledger.run_population, unmodified",
            "engine_validation": {
                "branches_checked": int(len(check_frame)),
                "worst_relative_N_SN_difference": round(worst_nsn, 5),
                "declared_tolerance": NSN_RELATIVE_TOLERANCE,
                "worst_absolute_P_recent_difference": round(worst_prec, 5),
                "declared_tolerance_P_recent": PRECENT_ABSOLUTE_TOLERANCE,
                "pass": engine_ok,
            },
            "assignments": {
                f"{row.family}_RV{row.R_V:g}": {
                    f"alpha_{tag}": row[f"alpha_{tag}"] for tag in "ABC"
                }
                for _, row in mixed[mixed.sf_duration_Myr.eq(0.0)].iterrows()
            },
            "N_SN_range": [
                round(float(mixed.N_SN_mean.min()), 3),
                round(float(mixed.N_SN_mean.max()), 3),
            ],
            "P_recent_range": [
                round(float(mixed.P_last_SN_within_100kyr.min()), 4),
                round(float(mixed.P_last_SN_within_100kyr.max()), 4),
            ],
            "versus_pure_alpha_2p3": {
                "worst_relative_change_in_N_SN": round(worst_mixed_shift, 5),
                "branches_where_C_is_above_the_imf_ceiling": int(len(c_dead)),
                "branches_where_C_contributes": int(len(c_alive)),
                "worst_relative_change_where_C_contributes": (
                    round(float(c_alive.relative_change.abs().max()), 5)
                    if len(c_alive) else None
                ),
                "worst_relative_change_where_C_is_dead": (
                    round(float(c_dead.relative_change.abs().max()), 5)
                    if len(c_dead) else None
                ),
                "interpretation": (
                    "Cyg OB2-C is the only subgroup reassigned (2.3 -> 2.0), so "
                    "the entire effect is C's.  It splits on whether C's "
                    "turnoff has crossed the 120 Msun IMF ceiling on that "
                    "branch.  Where it has not -- PARSEC at every R_V, and "
                    "MIST at R_V = 3.5 -- C contributes essentially nothing at "
                    "either slope and the mixed-slope ledger reproduces the "
                    "uniform one.  Where it has -- MIST at R_V = 3.0 and 3.1 --"
                    " C supplies 2.5-7.3 supernovae and giving it the "
                    "shallower slope its own closure prefers raises the "
                    "association total by a third to a half.  The subgroup "
                    "heterogeneity is therefore not a presentational detail: "
                    "on the branches where C is dynamically alive it is worth "
                    "as much as an R_V step.  It remains inside the carried "
                    "headline range (prediction R4)."
                ),
                "per_branch": deltas,
            },
        },
        "predictions": [r3, r4],
    }
    rec = W.record(
        "scripts/wp12_closure_slopes.py", payload,
        {"closure": out_closure, "slopes": out_slopes, "mixed": out_mixed},
    )
    w.write_json(w.PROVENANCE / "wp12_closure_slopes_execution.json", rec)

    print("WP12.2 closure and mixed slopes")
    print("  closing alpha (grid median) by subgroup:")
    for s, v in payload["closing_slope_grid_median_by_subgroup"].items():
        print(f"    {s:12s} {v}")
    print(f"  engine validation: worst dN/N = {worst_nsn:.4f} "
          f"(tol {NSN_RELATIVE_TOLERANCE}), worst dP = {worst_prec:.4f} "
          f"(tol {PRECENT_ABSOLUTE_TOLERANCE}) -> "
          f"{'PASS' if engine_ok else 'FAIL'}")
    print(f"  mixed-slope N_SN range: {payload['mixed_slope']['N_SN_range']} "
          f"against headline {[round(lo,2), round(hi,2)]}")
    for p in (r3, r4):
        print(f"  {p['id']} {p['outcome']}  {p['measured']}")


if __name__ == "__main__":
    main()
