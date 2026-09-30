#!/usr/bin/env python3
"""WP12.3 + WP12.4 -- the cocoon quantity, relabelled and made a function of its
weakest terms.

WP9 reports P_verdict = C1 x C3 x C4.  The arithmetic is right and the label is
wrong.  Three of the four things a calibrated probability needs are missing:

  C4 = 0.854 is an UPPER BOUND, not a sampled probability.  It comes from a
       two-dimensional, footprint-limited runaway census of LIVING stars below
       the turnoff, and it is applied to already-dead, more massive progenitors
       whose escape fraction need not match.
  C3 = 1.000 is a DETERMINISTIC MAPPING, not a measurement.  It is the
       indicator "every simulated progenitor exceeds 30 Msun", which is true,
       and it is read as "every explosion would have been observed as type
       Ib/c", which does not follow.
  independence of C1, C3 and C4 is ASSERTED, not demonstrated.

So this script does not compute a better probability.  It computes the same
product under an honest name --- a conditional scenario-availability score ---
and reports it against the two terms that are not measurements:

  WP12.3  the score as a function of C4 over [0.30, 1.00], with the two
          thresholds the manuscript currently conflates stated separately;
  WP12.4  the score under alternative progenitor-mass-to-supernova-type
          mappings: a threshold scan, a smooth logistic stripping probability,
          and an explicitly pessimistic 60 Msun case.

The population engine is WP7's, unmodified, run on the frozen repair_v7
posterior draws.  It is validated by reproducing the stored WP9 C1 and C3
before any alternative is reported.

Outputs:
  tables/wp12_c4_scan.csv
  tables/wp12_c3_subtype.csv
  tables/wp12_scenario_score.csv
  provenance/wp12_scenario_score_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/wp12_scenario_score.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import wp5_common as w
import wp12_common as W
from wp7_ledger import TurnoffRelation, draw_key, run_population
from wp7_ledger_prereg import SF_DURATIONS_MYR, SN_THRESHOLD_MSUN
from wp9_verdict_prereg import (
    AGE_PERMISSIVE_KYR,
    AGE_WINDOW_KYR,
    STRIPPED_PROGENITOR_MSUN,
)

ITERATIONS = 200_000
SEED = 20260804

C4_GRID = np.round(np.arange(0.30, 1.00001, 0.005), 4)
C4_ADOPTED = 0.854          # WP9's bound, retained as the reference point
SUPPORT = W.SUPPORT_THRESHOLD

# WP12.4 -- alternative mass-to-type mappings.  None of these is offered as a
# better model of envelope stripping.  They are a bracket: the point is that
# C3 = 1 is a property of the 30 Msun step, not of the population.
THRESHOLD_SCAN = (20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0)
LOGISTIC_GRID = tuple(
    (m0, dm) for m0 in (30.0, 45.0, 60.0) for dm in (5.0, 10.0)
)
PESSIMISTIC_MSUN = 60.0

# Declared tolerances for reproducing the stored 2e6-iteration WP9 values at
# 2e5 iterations here.
C1_ABSOLUTE_TOLERANCE = 0.01
C3_ABSOLUTE_TOLERANCE = 1e-9    # an indicator over a population with a floor


def logistic(mass: np.ndarray, m0: float, dm: float) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-(mass - m0) / dm))


def main() -> None:
    draws = np.load(W.frozen("wp5_posterior_draws"))
    relations = {family: TurnoffRelation(family) for family in w.FAMILIES}
    master = np.random.default_rng(SEED)

    lo, hi = AGE_WINDOW_KYR[0] / 1000.0, AGE_WINDOW_KYR[1] / 1000.0
    plo, phi = AGE_PERMISSIVE_KYR[0] / 1000.0, AGE_PERMISSIVE_KYR[1] / 1000.0

    rows, subtype_rows = [], []
    for family in w.FAMILIES:
        relation = relations[family]
        for rv in w.R_V_BRANCHES:
            for alpha in W.HEADLINE_ALPHAS:
                for delta in SF_DURATIONS_MYR:
                    in_window = np.zeros(ITERATIONS, dtype=int)
                    in_permissive = np.zeros(ITERATIONS, dtype=int)
                    sn_masses = []
                    for subgroup in w.SUBGROUPS:
                        key = draw_key(subgroup, family, rv, alpha)
                        k_all = draws[f"k__{key}"]
                        age_all = draws[f"truth_age_draws__{key}"]
                        rng = np.random.default_rng(
                            master.integers(0, 2 ** 63 - 1)
                        )
                        pick = rng.integers(0, k_all.size, ITERATIONS)
                        res = run_population(
                            rng, k_all[pick], age_all[pick], alpha, delta,
                            relation,
                        )
                        mass = res["dead_masses"]
                        epoch = res["epochs"]
                        it = res["dead_iteration"]
                        explodes = mass >= SN_THRESHOLD_MSUN
                        if explodes.any():
                            sel = explodes & (epoch >= lo) & (epoch <= hi)
                            in_window += np.bincount(it[sel], minlength=ITERATIONS)
                            selp = explodes & (epoch >= plo) & (epoch <= phi)
                            in_permissive += np.bincount(
                                it[selp], minlength=ITERATIONS)
                            sn_masses.append(mass[explodes])

                    masses = (
                        np.concatenate(sn_masses) if sn_masses else np.array([])
                    )
                    c1 = float((in_window >= 1).mean())
                    c1p = float((in_permissive >= 1).mean())

                    # ------------------------------------------ WP12.4  C3 maps
                    maps: dict[str, float] = {}
                    if masses.size:
                        maps["step_30_baseline"] = float(
                            (masses > STRIPPED_PROGENITOR_MSUN).mean())
                        for t in THRESHOLD_SCAN:
                            maps[f"step_{t:g}"] = float((masses > t).mean())
                        for m0, dm in LOGISTIC_GRID:
                            maps[f"logistic_m0_{m0:g}_dm_{dm:g}"] = float(
                                logistic(masses, m0, dm).mean())
                        maps["pessimistic_60"] = float(
                            (masses > PESSIMISTIC_MSUN).mean())
                    else:
                        maps = {k: float("nan") for k in
                                ["step_30_baseline", "pessimistic_60"]}

                    branch = dict(family=family, R_V=rv, alpha=alpha,
                                  sf_duration_Myr=delta)
                    for name, value in maps.items():
                        subtype_rows.append({**branch, "mapping": name,
                                             "C3": round(value, 6)})

                    c3_base = maps["step_30_baseline"]
                    c3_pess = maps["pessimistic_60"]
                    alt = [v for k, v in maps.items() if k != "step_30_baseline"]
                    rows.append({
                        **branch,
                        "C1_age": round(c1, 5),
                        "C1_age_permissive": round(c1p, 5),
                        "C3_baseline_step30": round(c3_base, 5),
                        "C3_pessimistic_step60": round(c3_pess, 5),
                        "C3_alt_min": round(float(np.min(alt)), 5),
                        "C3_alt_max": round(float(np.max(alt)), 5),
                        "C4_adopted_bound": C4_ADOPTED,
                        "score_baseline": round(c1 * c3_base * C4_ADOPTED, 5),
                        "score_permissive": round(c1p * c3_base * C4_ADOPTED, 5),
                        "score_pessimistic_C3": round(
                            c1 * c3_pess * C4_ADOPTED, 5),
                        "n_supernovae_sampled": int(masses.size),
                        "min_progenitor_Msun": (
                            round(float(masses.min()), 3) if masses.size else None
                        ),
                        "median_progenitor_Msun": (
                            round(float(np.median(masses)), 3)
                            if masses.size else None
                        ),
                    })
            print(f"  {family} R_V={rv} done", flush=True)

    score = pd.DataFrame(rows)
    subtype = pd.DataFrame(subtype_rows)

    # ------------------------------------------- validation against stored WP9
    stored = W.verdict()
    stored = stored[stored.in_headline_set].set_index(
        ["family", "R_V", "alpha", "sf_duration_Myr"])
    checks = []
    for row in score.itertuples():
        ref = stored.loc[(row.family, row.R_V, row.alpha, row.sf_duration_Myr)]
        checks.append({
            "family": row.family, "R_V": row.R_V, "alpha": row.alpha,
            "sf_duration_Myr": row.sf_duration_Myr,
            "wp9_C1": float(ref.C1_age), "wp12_C1": row.C1_age,
            "dC1": round(abs(row.C1_age - float(ref.C1_age)), 5),
            "wp9_C3": float(ref.C3_stripped_fraction),
            "wp12_C3": row.C3_baseline_step30,
            "dC3": round(abs(row.C3_baseline_step30
                             - float(ref.C3_stripped_fraction)), 9),
        })
    check = pd.DataFrame(checks)
    worst_c1 = float(check.dC1.max())
    worst_c3 = float(check.dC3.max())
    engine_ok = bool(worst_c1 <= C1_ABSOLUTE_TOLERANCE
                     and worst_c3 <= C3_ABSOLUTE_TOLERANCE)

    # ------------------------------------------------------- WP12.3  C4 scan
    # The scan uses the STORED WP9 C1 and C3 so that the published scan is a
    # transformation of the published verdict table, not of this script's
    # re-run.  The re-run exists to validate, not to replace.
    head = W.verdict()
    head = head[head.in_headline_set].copy()
    head["C1C3"] = head.C1_age * head.C3_stripped_fraction
    gate = W.combination_pass().set_index(["family", "R_V", "alpha"])
    head["all_subgroup_pass"] = [
        bool(gate.loc[(r.family, r.R_V, r.alpha)].all_subgroup_pass)
        for r in head.itertuples()
    ]

    scan_rows = []
    for row in head.itertuples():
        for c4 in C4_GRID:
            scan_rows.append({
                "family": row.family, "R_V": row.R_V, "alpha": row.alpha,
                "sf_duration_Myr": row.sf_duration_Myr,
                "all_subgroup_pass": row.all_subgroup_pass,
                "C4": float(c4),
                "score": round(float(row.C1C3 * c4), 6),
                "above_support_threshold": bool(row.C1C3 * c4 > SUPPORT),
            })
    scan = pd.DataFrame(scan_rows)

    head["C4_critical"] = SUPPORT / head.C1C3
    a20 = head[head.alpha.eq(2.0)]
    a23 = head[head.alpha.eq(2.3)]
    c4_any = float(a20.C4_critical.max())    # below this, >=1 branch drops
    c4_all = float(a20.C4_critical.min())    # below this, ALL branches drop
    # And the other direction: what C4 would be needed to lift alpha=2.3.
    c4_lift_first = float(a23.C4_critical.min())
    c4_lift_all = float(a23.C4_critical.max())

    strict20 = a20[a20.all_subgroup_pass]
    strict23 = a23[a23.all_subgroup_pass]

    # --------------------------------------------------------- prediction R5
    r5_measured = {
        "min_C3_at_60_Msun_over_headline_branches": round(
            float(score.C3_pessimistic_step60.min()), 5),
        "max_C3_at_60_Msun_over_headline_branches": round(
            float(score.C3_pessimistic_step60.max()), 5),
        "C3_at_30_Msun_is_exactly_one_everywhere": bool(
            (score.C3_baseline_step30 >= 1.0 - 1e-12).all()),
    }
    r5 = W.score_prediction(
        "R5", bool(score.C3_pessimistic_step60.min() < 0.95), r5_measured)

    out_scan = w.TABLES / "wp12_c4_scan.csv"
    out_subtype = w.TABLES / "wp12_c3_subtype.csv"
    out_score = w.TABLES / "wp12_scenario_score.csv"
    scan.to_csv(out_scan, index=False)
    subtype.to_csv(out_subtype, index=False)
    score.to_csv(out_score, index=False)

    payload = {
        "item": "WP12.3 + WP12.4 -- conditional scenario-availability score",
        "rename": {
            "from": "P_verdict",
            "to": "conditional scenario-availability score S",
            "why": (
                "S = C1 x C3 x C4 multiplies a sampled Monte-Carlo probability "
                "(C1) by a deterministic model mapping (C3) and an upper bound "
                "(C4), assuming independence.  It is a scenario score, not a "
                "calibrated probability, and the manuscript must not call it "
                "one."
            ),
            "what_S_is": (
                "the probability, under the adopted branch and CONDITIONAL on "
                "the adopted type mapping and on the in-situ bound holding "
                "with equality, that Cyg OB2 made an explosion of the required "
                "age, type and location"
            ),
        },
        "engine_validation": {
            "iterations": ITERATIONS,
            "seed": SEED,
            "reference": "tables/wp9_verdict.csv at 2e6 iterations",
            "worst_absolute_C1_difference": round(worst_c1, 5),
            "declared_tolerance_C1": C1_ABSOLUTE_TOLERANCE,
            "worst_absolute_C3_difference": round(worst_c3, 12),
            "declared_tolerance_C3": C3_ABSOLUTE_TOLERANCE,
            "pass": engine_ok,
        },
        "wp12_3_c4_sensitivity": {
            "adopted_C4": C4_ADOPTED,
            "adopted_C4_status": (
                "an UPPER BOUND on the in-situ fraction from a 2-D, "
                "footprint-limited runaway census of living stars below the "
                "turnoff, applied to dead, more massive progenitors"
            ),
            "C4_below_which_at_least_one_alpha2p0_branch_falls_below_0p5":
                round(c4_any, 4),
            "C4_below_which_every_alpha2p0_branch_falls_below_0p5":
                round(c4_all, 4),
            "C4_above_which_the_first_alpha2p3_branch_would_rise_above_0p5":
                round(c4_lift_first, 4),
            "C4_above_which_every_alpha2p3_branch_would_rise_above_0p5":
                round(c4_lift_all, 4),
            "alpha2p3_cannot_be_lifted_by_C4_alone": bool(c4_lift_first > 1.0),
            "correction_to_the_manuscript": (
                "the current text says C4 would have to fall below 0.60 to "
                "move ANY alpha = 2.0 branch under 0.5.  0.60 is close to the "
                "threshold at which ALL of them do (0.5806); the threshold at "
                "which the FIRST one does is 0.7198.  The two must be stated "
                "separately."
            ),
            "strict_subset": {
                "alpha2p0_all_subgroup_pass_branches": int(len(strict20)),
                "C4_any_strict": (
                    round(float(strict20.C4_critical.max()), 4)
                    if len(strict20) else None
                ),
                "C4_all_strict": (
                    round(float(strict20.C4_critical.min()), 4)
                    if len(strict20) else None
                ),
                "alpha2p3_all_subgroup_pass_branches": int(len(strict23)),
            },
        },
        "wp12_4_c3_subtype": {
            "baseline_mapping": (
                f"deterministic step at {STRIPPED_PROGENITOR_MSUN:g} Msun"),
            "baseline_C3_is_exactly_one": bool(
                (score.C3_baseline_step30 >= 1.0 - 1e-12).all()),
            "why_it_is_one": (
                "every sampled progenitor exceeds the threshold.  That is a "
                "fact about the ledger's mass floor, not evidence that every "
                "explosion would have been observed as type Ib/c."
            ),
            "threshold_scan_C3_range": {
                f"step_{t:g}": [
                    round(float(subtype[subtype.mapping.eq(f"step_{t:g}")].C3.min()), 4),
                    round(float(subtype[subtype.mapping.eq(f"step_{t:g}")].C3.max()), 4),
                ]
                for t in THRESHOLD_SCAN
            },
            "logistic_C3_range": {
                f"logistic_m0_{m0:g}_dm_{dm:g}": [
                    round(float(subtype[subtype.mapping.eq(
                        f"logistic_m0_{m0:g}_dm_{dm:g}")].C3.min()), 4),
                    round(float(subtype[subtype.mapping.eq(
                        f"logistic_m0_{m0:g}_dm_{dm:g}")].C3.max()), 4),
                ]
                for m0, dm in LOGISTIC_GRID
            },
            "pessimistic_60_Msun_C3_range": [
                round(float(score.C3_pessimistic_step60.min()), 4),
                round(float(score.C3_pessimistic_step60.max()), 4),
            ],
            "score_under_pessimistic_C3": {
                "range": [round(float(score.score_pessimistic_C3.min()), 4),
                          round(float(score.score_pessimistic_C3.max()), 4)],
                "alpha2p0_above_0p5": int(
                    (score[score.alpha.eq(2.0)].score_pessimistic_C3
                     > SUPPORT).sum()),
                "alpha2p0_branches": int((score.alpha.eq(2.0)).sum()),
                "alpha2p3_above_0p5": int(
                    (score[score.alpha.eq(2.3)].score_pessimistic_C3
                     > SUPPORT).sum()),
            },
            "progenitor_mass_floor": {
                "min_over_headline_branches_Msun": round(
                    float(score.min_progenitor_Msun.min()), 3),
                "max_over_headline_branches_Msun": round(
                    float(score.min_progenitor_Msun.max()), 3),
            },
            "sourcing_note": (
                "The dependence of envelope stripping on wind mass loss, "
                "metallicity, rotation and binarity is why a single ZAMS-mass "
                "threshold cannot be a measurement.  The manuscript cites "
                "primary sources for that dependence where it quotes C3."
            ),
        },
        "predictions": [r5],
    }
    rec = W.record(
        "scripts/wp12_scenario_score.py", payload,
        {"scan": out_scan, "subtype": out_subtype, "score": out_score},
    )
    w.write_json(w.PROVENANCE / "wp12_scenario_score_execution.json", rec)

    print("WP12.3/12.4 conditional scenario-availability score")
    print(f"  engine validation: worst dC1 = {worst_c1:.5f} "
          f"(tol {C1_ABSOLUTE_TOLERANCE}), worst dC3 = {worst_c3:.2e} -> "
          f"{'PASS' if engine_ok else 'FAIL'}")
    print(f"  C4 thresholds: first branch drops at {c4_any:.4f}, "
          f"all drop at {c4_all:.4f}")
    print(f"  C4 needed to lift the first alpha=2.3 branch: {c4_lift_first:.4f}"
          f"  (>1 means impossible: {c4_lift_first > 1.0})")
    print(f"  C3 at 30 Msun = 1.000 everywhere: "
          f"{payload['wp12_4_c3_subtype']['baseline_C3_is_exactly_one']}")
    print(f"  C3 at 60 Msun range: "
          f"{payload['wp12_4_c3_subtype']['pessimistic_60_Msun_C3_range']}")
    print(f"  score under pessimistic C3: "
          f"{payload['wp12_4_c3_subtype']['score_under_pessimistic_C3']}")
    print(f"  {r5['id']} {r5['outcome']}  {r5['measured']}")


if __name__ == "__main__":
    main()
