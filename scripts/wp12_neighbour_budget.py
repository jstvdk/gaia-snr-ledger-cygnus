#!/usr/bin/env python3
"""WP12.5 -- the neighbouring-association supernova budget (the plan's WP8.5).

The original execution plan called for a coarse, literature-based supernova
budget for the other Cygnus populations, so that the paper could bound the
probability that an explosion somewhere in the Cygnus X cavity came from
Cyg OB2 rather than from a neighbour.  WP8 explicitly declined to claim one
("No quantitative budget is claimed").  The manuscript nevertheless reads as
though the question had been addressed.  This script closes that gap.

SOURCE.  Martin, Knodlseder, Meynet & Diehl 2010, A&A 511, A86, Table 1,
"Characteristics of the Cygnus OB associations and stellar clusters
collectively referred to as the Cygnus complex", reproduced there from
Knodlseder, Cervino, Le Duigou, Meynet, Schaerer & von Ballmoos 2002, A&A 390,
945.  Each row gives an OBSERVED massive-star count inside a stated ZAMS mass
interval, together with a distance and an age.  Those three numbers are all
this calculation uses.

METHOD.  Normalize a single power law to the observed count over the quoted
interval (truncated at the turnoff, since stars above it are already gone),
integrate above the turnoff to get the cumulative number of deaths, and
difference it over a 100 kyr window to get the recent rate.  The turnoff
relation is THIS PROJECT'S OWN --- the same `turnoff_mass` the ledger uses ---
so the coarse estimator and the measured ledger differ in their inputs, not in
their stellar physics.

WHAT THIS IS NOT.  It is not a census, it is not completeness-corrected, and
its ages are heterogeneous literature values from 2002.  Every number it
produces is a lower bound with an age systematic that dominates it.  The
declared coarseness list in the preregistration is reproduced in the output.

Outputs:
  tables/wp12_neighbour_budget.csv
  tables/wp12_cavity_share.csv
  provenance/wp12_neighbour_budget_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/wp12_neighbour_budget.py
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import wp5_common as w
import wp12_common as W
from wp6_mass_extension_decision import IMF_UPPER_LIMIT, turnoff_mass
from wp7_ledger_prereg import RECENT_WINDOW_MYR

# --------------------------------------------------------------------- source
SOURCE = (
    "Martin et al. 2010, A&A 511, A86, Table 1 (from Knodlseder et al. 2002, "
    "A&A 390, 945)"
)

# name, observed count, mass interval lo/hi (Msun), distance (pc), age (Myr)
CYGNUS_COMPLEX = [
    ("Cyg OB1",   23, 15.0,  40.0, 1905,  4.0),
    ("Cyg OB2",  120, 20.0, 120.0, 1584,  2.5),
    ("Cyg OB3",   14, 25.0,  60.0, 2187,  3.5),
    ("Cyg OB7",   10,  7.0,  25.0,  832,  3.5),
    ("Cyg OB8",    6, 20.0,  40.0, 2399,  7.5),
    ("Cyg OB9",    8, 20.0,  40.0, 1259,  3.5),
    ("Ber 86",    11,  7.0,  25.0, 1660,  4.0),
    ("Ber 87",    24,  7.0,  25.0, 1905,  4.5),
    ("NGC 6871",  13,  7.0,  25.0, 2399,  5.5),
    ("NGC 6913",  13,  7.0,  40.0, 1820,  3.5),
    ("NGC 6910",   7,  7.0,  25.0, 1820,  4.5),
    ("NGC 6883",   2, 12.0,  15.0, 1820, 15.0),
    ("IC 4996",    6,  6.0,  25.0, 1660,  5.5),
]

# Which populations plausibly share the Cygnus X cavity with Cyg OB2.  This is
# a projection judgement, not a measurement, so two nested sets are carried and
# the answer is reported for both.
STRICT_SET = ["Cyg OB1", "Cyg OB9"]
WIDE_SET = STRICT_SET + [
    "Cyg OB3", "Cyg OB8", "Ber 86", "Ber 87", "NGC 6871", "NGC 6913",
    "NGC 6910", "IC 4996",
]
# Cyg OB7 is excluded from both: at 832 pc it is a foreground object, not part
# of the complex.  NGC 6883 is excluded from STRICT and kept in WIDE only
# through its 15 Myr age, where its turnoff is far below the others'.
EXCLUDED = ["Cyg OB7"]

# An illustrative age bracket.  NOT a literature uncertainty -- Knodlseder's
# ages are quoted without one.  Declared as an adopted bracket so the reader
# can see how much of the answer is the assumed age.
AGE_SCALES = (0.75, 1.0, 1.35)


def imf_count(alpha: float, lo: float, hi: float) -> float:
    """integral[lo, hi] M^-alpha dM, zero if the interval is empty."""
    if hi <= lo:
        return 0.0
    return (lo ** (1.0 - alpha) - hi ** (1.0 - alpha)) / (alpha - 1.0)


def population_budget(
    n_obs: int, m_lo: float, m_hi: float, age: float,
    family: str, alpha: float, window: float = RECENT_WINDOW_MYR,
) -> dict:
    """Cumulative and recent supernova counts for one literature population."""
    turnoff_now = min(turnoff_mass(family, age), IMF_UPPER_LIMIT)
    # Stars above the turnoff are gone, so the observed interval is truncated
    # there before the normalization is taken.
    obs_hi = min(m_hi, turnoff_now)
    denom = imf_count(alpha, m_lo, obs_hi)
    if denom <= 0:
        return {
            "k": float("nan"), "turnoff_Msun": turnoff_now,
            "N_dead_cumulative": float("nan"),
            "N_dead_recent_window": float("nan"),
            "usable": False,
            "why_unusable": (
                "the quoted mass interval lies entirely above the turnoff at "
                "the quoted age, so the count cannot normalize an IMF"
            ),
        }
    k = n_obs / denom
    cumulative = k * imf_count(alpha, turnoff_now, IMF_UPPER_LIMIT)
    earlier = max(age - window, 1e-6)
    turnoff_then = min(turnoff_mass(family, earlier), IMF_UPPER_LIMIT)
    recent = k * imf_count(alpha, turnoff_now, turnoff_then)
    return {
        "k": k, "turnoff_Msun": turnoff_now,
        "turnoff_one_window_ago_Msun": turnoff_then,
        "N_dead_cumulative": cumulative,
        "N_dead_recent_window": recent,
        "usable": True,
        "why_unusable": "",
    }


def measured_recent_expectation() -> tuple[float, float, float]:
    """Expected number of Cyg OB2 supernovae in the last 100 kyr, measured.

    From the ledger's own R_SN(t) curves on the baseline branch, so the
    comparison against the coarse neighbour numbers is like for like: both are
    expected counts in the same window.
    """
    rsn = pd.read_csv(w.ROOT / "tables" / "wp7_rsn_curves.csv")
    base = rsn[
        rsn.family.eq(W.BASE["family"]) & rsn.R_V.eq(W.BASE["R_V"])
        & rsn.alpha.eq(W.BASE["alpha"])
        & rsn.sf_duration_Myr.eq(W.BASE["sf_duration_Myr"])
    ] if "family" in rsn.columns else rsn
    width = float(np.median(np.diff(sorted(base.lookback_lo_Myr.unique()))))
    recent = base[base.lookback_lo_Myr < RECENT_WINDOW_MYR - 1e-9]
    expected = float(recent.rate_per_Myr.sum() * width)
    total = float(base.rate_per_Myr.sum() * width)
    return expected, total, width


def main() -> None:
    rows = []
    for name, n_obs, m_lo, m_hi, dist, age in CYGNUS_COMPLEX:
        for scale in AGE_SCALES:
            t = age * scale
            for family in w.FAMILIES:
                for alpha in w.IMF_SLOPES:
                    got = population_budget(n_obs, m_lo, m_hi, t, family, alpha)
                    rows.append({
                        "population": name, "observed_count": n_obs,
                        "mass_lo_Msun": m_lo, "mass_hi_Msun": m_hi,
                        "distance_pc": dist, "quoted_age_Myr": age,
                        "age_scale": scale, "age_used_Myr": round(t, 4),
                        "family": family, "alpha": alpha,
                        **{k: (round(v, 6) if isinstance(v, float) else v)
                           for k, v in got.items()},
                        "in_strict_cavity_set": name in STRICT_SET,
                        "in_wide_cavity_set": name in WIDE_SET,
                        "excluded_from_cavity": name in EXCLUDED,
                    })
    budget = pd.DataFrame(rows)
    budget.to_csv(w.TABLES / "wp12_neighbour_budget.csv", index=False)

    # ------------------------------------------------ R6: validate the method
    measured_recent, measured_total, width = measured_recent_expectation()
    stored = W.association_all_explode()
    base_nsn = float(stored[
        stored.family.eq(W.BASE["family"]) & stored.R_V.eq(W.BASE["R_V"])
        & stored.alpha.eq(W.BASE["alpha"])
        & stored.sf_duration_Myr.eq(W.BASE["sf_duration_Myr"])
    ].iloc[0].N_SN_mean)

    own = budget[
        budget.population.eq("Cyg OB2") & budget.age_scale.eq(1.0)
    ]
    own_usable = own[own.usable]
    ratios = (own_usable.N_dead_cumulative / base_nsn).to_numpy()
    r6_pass = bool(
        len(ratios) > 0 and np.all((ratios >= 0.1) & (ratios <= 10.0))
    )
    # A non-preregistered diagnostic, labelled as such: the same estimator run
    # at OUR measured age instead of Knodlseder's 2.5 Myr.  This separates "the
    # estimator is wrong" from "the 2002 age is wrong".
    diagnostic = []
    for family in w.FAMILIES:
        for alpha in w.IMF_SLOPES:
            got = population_budget(120, 20.0, 120.0, 4.0, family, alpha)
            diagnostic.append({
                "family": family, "alpha": alpha,
                "N_dead_cumulative": round(got["N_dead_cumulative"], 4)
                if got["usable"] else None,
                "ratio_to_measured_baseline": round(
                    got["N_dead_cumulative"] / base_nsn, 4)
                if got["usable"] else None,
            })
    diag = pd.DataFrame(diagnostic)
    diag_ok = bool(
        diag.ratio_to_measured_baseline.notna().all()
        and ((diag.ratio_to_measured_baseline >= 0.1)
             & (diag.ratio_to_measured_baseline <= 10.0)).all()
    )

    r6 = W.score_prediction("R6", r6_pass, {
        "measured_baseline_N_SN": round(base_nsn, 3),
        "coarse_at_Knodlseder_age_2p5_Myr": {
            f"{r.family}_alpha{r.alpha:g}": (
                round(r.N_dead_cumulative, 4) if r.usable else "UNUSABLE"
            )
            for r in own.itertuples()
        },
        "usable_cells": int(own.usable.sum()),
        "total_cells": int(len(own)),
        "ratio_range": (
            [round(float(ratios.min()), 4), round(float(ratios.max()), 4)]
            if len(ratios) else None
        ),
        "non_preregistered_diagnostic_at_our_measured_4p0_Myr": {
            "rows": diagnostic,
            "all_within_factor_10": diag_ok,
            "label": (
                "NOT part of R6.  Recorded to separate 'the coarse estimator "
                "is wrong' from 'the 2002 age for Cyg OB2 is wrong'."
            ),
        },
    })

    # ---------------------------------------------------- the cavity share
    share_rows = []
    for scale in AGE_SCALES:
        for family in w.FAMILIES:
            for alpha in w.IMF_SLOPES:
                cell = budget[
                    budget.age_scale.eq(scale) & budget.family.eq(family)
                    & budget.alpha.eq(alpha)
                ]
                for label, members in (("strict", STRICT_SET), ("wide", WIDE_SET)):
                    nb = cell[cell.population.isin(members)]
                    nb_recent = float(
                        nb.N_dead_recent_window.fillna(0.0).sum())
                    nb_cumulative = float(
                        nb.N_dead_cumulative.fillna(0.0).sum())
                    unusable = int((~nb.usable).sum())
                    denom = measured_recent + nb_recent
                    share_rows.append({
                        "cavity_set": label, "age_scale": scale,
                        "family": family, "alpha": alpha,
                        "n_neighbours": int(len(nb)),
                        "n_neighbours_unusable": unusable,
                        "cygob2_measured_recent_expected": round(
                            measured_recent, 5),
                        "neighbours_coarse_recent_expected": round(nb_recent, 5),
                        "neighbours_coarse_cumulative": round(nb_cumulative, 4),
                        "f_cygob2_recent": (
                            round(measured_recent / denom, 5)
                            if denom > 0 else float("nan")
                        ),
                    })
    share = pd.DataFrame(share_rows)
    share.to_csv(w.TABLES / "wp12_cavity_share.csv", index=False)

    baseline_share = share[
        share.cavity_set.eq("wide") & share.age_scale.eq(1.0)
        & share.family.eq(W.BASE["family"]) & share.alpha.eq(W.BASE["alpha"])
    ]
    f_wide = float(baseline_share.f_cygob2_recent.iloc[0])
    r7 = W.score_prediction("R7", bool(f_wide < 0.9), {
        "f_cygob2_recent_wide_baseline": round(f_wide, 5),
        "f_range_over_all_cells_wide": [
            round(float(share[share.cavity_set.eq("wide")].f_cygob2_recent.min()), 4),
            round(float(share[share.cavity_set.eq("wide")].f_cygob2_recent.max()), 4),
        ],
        "f_range_over_all_cells_strict": [
            round(float(share[share.cavity_set.eq("strict")].f_cygob2_recent.min()), 4),
            round(float(share[share.cavity_set.eq("strict")].f_cygob2_recent.max()), 4),
        ],
    })

    reportable = r6["outcome"] == "PASS"
    payload = {
        "item": "WP12.5 -- neighbouring-association supernova budget (plan WP8.5)",
        "source": SOURCE,
        "source_table_rows": [
            {"population": n, "observed_count": c, "mass_interval_Msun": [lo, hi],
             "distance_pc": d, "age_Myr": a}
            for n, c, lo, hi, d, a in CYGNUS_COMPLEX
        ],
        "method": (
            "normalize a single power law to the observed count over the "
            "quoted interval truncated at the turnoff; integrate above the "
            "turnoff for the cumulative deaths; difference over 100 kyr for "
            "the recent rate.  Turnoff from this project's own relation."
        ),
        "declared_coarseness": W.prereg()["wp12_5_neighbour_budget"][
            "declared_coarseness"],
        "age_bracket": {
            "scales": list(AGE_SCALES),
            "status": (
                "ADOPTED ILLUSTRATIVE BRACKET, not a literature uncertainty.  "
                "Knodlseder et al. quote ages without errors."
            ),
        },
        "cavity_sets": {
            "strict": STRICT_SET, "wide": WIDE_SET, "excluded": EXCLUDED,
            "status": "a projection judgement, not a measurement",
        },
        "measured_cygob2_reference": {
            "expected_SNe_in_last_100kyr_baseline": round(measured_recent, 5),
            "expected_SNe_total_over_rsn_grid": round(measured_total, 4),
            "baseline_N_SN_mean": round(base_nsn, 4),
            "rsn_bin_width_Myr": width,
        },
        "cavity_share_baseline_wide": round(f_wide, 5),
        "cavity_share_ranges": {
            "wide": r7["measured"]["f_range_over_all_cells_wide"],
            "strict": r7["measured"]["f_range_over_all_cells_strict"],
        },
        "gate": {
            "R6_is_the_gate": True,
            "reportable_as_a_number": reportable,
            "consequence_if_failed": (
                "per the preregistration, if R6 fails the neighbour budget is "
                "NOT reported as a number and the manuscript scopes the claim "
                "out instead, stating that it estimates only the availability "
                "of an event from Cyg OB2 and cannot quantify whether Cyg OB2 "
                "was the source of an event elsewhere in the cavity."
            ),
        },
        "predictions": [r6, r7],
    }
    rec = W.record(
        "scripts/wp12_neighbour_budget.py", payload,
        {"budget": w.TABLES / "wp12_neighbour_budget.csv",
         "share": w.TABLES / "wp12_cavity_share.csv"},
    )
    w.write_json(w.PROVENANCE / "wp12_neighbour_budget_execution.json", rec)

    print("WP12.5 neighbouring-association budget")
    print(f"  source: {SOURCE}")
    print(f"  measured Cyg OB2 expectation in last 100 kyr: "
          f"{measured_recent:.4f} supernovae")
    print(f"  R6 (gate) {r6['outcome']}")
    print(f"    coarse Cyg OB2 at Knodlseder's 2.5 Myr: "
          f"{r6['measured']['coarse_at_Knodlseder_age_2p5_Myr']}")
    print(f"    diagnostic at our measured 4.0 Myr, within factor 10: "
          f"{diag_ok}")
    print(f"  R7 {r7['outcome']}  f_CygOB2(wide, baseline) = {f_wide:.4f}")
    print(f"  reportable as a number: {reportable}")


if __name__ == "__main__":
    main()
