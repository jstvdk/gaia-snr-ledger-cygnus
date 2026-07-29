#!/usr/bin/env python3
"""WP9 — the verdict.

Pre-registered in provenance/wp9_verdict_prereg.json.  Predictions V1-V3 and the
framing rule are fixed there and are applied, not amended.

The population engine is IMPORTED from wp7_ledger rather than reimplemented, so
the verdict cannot drift from the ledger it is derived from.

Outputs:
  tables/wp9_verdict.csv            per-branch verdict terms
  tables/wp9_sensitivity.csv        verdict vs each branch axis
  provenance/wp9_verdict_execution.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp9_verdict.py
"""
from __future__ import annotations

import argparse
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
from wp7_ledger import TurnoffRelation, draw_key, run_population
from wp7_ledger_prereg import SF_DURATIONS_MYR, SN_THRESHOLD_MSUN
from wp9_verdict_prereg import (
    AGE_PERMISSIVE_KYR,
    AGE_WINDOW_KYR,
    EXCLUDED_ALPHA,
    HEADLINE_ALPHAS,
    IMPLAUSIBLE_THRESHOLD,
    PLAUSIBLE_THRESHOLD,
    PRIMARY_EXPLODABILITY,
    STRIPPED_PROGENITOR_MSUN,
)

WP5_VERSION = "repair_v7"
ITERATIONS = 2_000_000
# WP7's in-situ bound: at most 14.6% of supernovae occurred outside the
# association, from the runaway fraction.  A bound, applied as a constant.
IN_SITU_FRACTION = 1.0 - 0.146


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--iterations", type=int, default=ITERATIONS)
    args = parser.parse_args()
    n_iter = int(args.iterations)

    draws = np.load(w.PROC / f"wp5_imf_posterior_draws_{WP5_VERSION}.npz")
    relations = {family: TurnoffRelation(family) for family in w.FAMILIES}
    master = np.random.default_rng(w.SEED)

    lo, hi = AGE_WINDOW_KYR[0] / 1000.0, AGE_WINDOW_KYR[1] / 1000.0
    plo, phi = AGE_PERMISSIVE_KYR[0] / 1000.0, AGE_PERMISSIVE_KYR[1] / 1000.0

    rows = []
    for family in w.FAMILIES:
        relation = relations[family]
        for rv in w.R_V_BRANCHES:
            for alpha in w.IMF_SLOPES:
                for delta in SF_DURATIONS_MYR:
                    in_window = np.zeros(n_iter, dtype=int)
                    in_permissive = np.zeros(n_iter, dtype=int)
                    stripped = total_sn = 0
                    for subgroup in w.SUBGROUPS:
                        key = draw_key(subgroup, family, rv, alpha)
                        k_all = draws[f"k__{key}"]
                        age_all = draws[f"truth_age_draws__{key}"]
                        rng = np.random.default_rng(
                            master.integers(0, 2 ** 63 - 1)
                        )
                        pick = rng.integers(0, k_all.size, n_iter)
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
                            in_window += np.bincount(it[sel], minlength=n_iter)
                            selp = explodes & (epoch >= plo) & (epoch <= phi)
                            in_permissive += np.bincount(
                                it[selp], minlength=n_iter
                            )
                            total_sn += int(explodes.sum())
                            stripped += int(
                                (explodes
                                 & (mass > STRIPPED_PROGENITOR_MSUN)).sum()
                            )

                    c1 = float((in_window >= 1).mean())
                    c1p = float((in_permissive >= 1).mean())
                    c3 = float(stripped / total_sn) if total_sn else float("nan")
                    c4 = IN_SITU_FRACTION
                    rows.append(
                        {
                            "family": family, "R_V": rv, "alpha": alpha,
                            "sf_duration_Myr": delta,
                            "explodability": PRIMARY_EXPLODABILITY,
                            "C1_age": round(c1, 5),
                            "C1_age_permissive": round(c1p, 5),
                            "C3_stripped_fraction": round(c3, 5),
                            "C4_in_situ": round(c4, 5),
                            "P_verdict": round(c1 * c3 * c4, 5),
                            "P_verdict_permissive": round(c1p * c3 * c4, 5),
                            "in_headline_set": alpha in HEADLINE_ALPHAS,
                        }
                    )
                print(f"  {family} R_V={rv} alpha={alpha} done", flush=True)

    table = pd.DataFrame(rows)
    table.to_csv(w.TABLES / "wp9_verdict.csv", index=False)

    head = table[table.in_headline_set]
    excluded = table[~table.in_headline_set]

    # ---- framing rule, applied mechanically -------------------------------
    if bool((head.P_verdict >= PLAUSIBLE_THRESHOLD).all()):
        outcome = "SUPPORTED"
    elif bool((head.P_verdict <= IMPLAUSIBLE_THRESHOLD).all()):
        outcome = "DISFAVOURED"
    else:
        outcome = "INCONCLUSIVE"
    framing = "Letter" if outcome in ("SUPPORTED", "DISFAVOURED") else "regular article"

    # ---- sensitivity: spread contributed by each axis ---------------------
    sensitivity = []
    for axis in ("alpha", "family", "R_V", "sf_duration_Myr"):
        grouped = head.groupby(axis).P_verdict.median()
        sensitivity.append(
            {
                "axis": axis,
                "n_levels": int(len(grouped)),
                "P_verdict_min": round(float(grouped.min()), 4),
                "P_verdict_max": round(float(grouped.max()), 4),
                "spread": round(float(grouped.max() - grouped.min()), 4),
                "levels": {f"{k}": round(float(v), 4)
                           for k, v in grouped.items()},
            }
        )
    sensitivity.sort(key=lambda r: -r["spread"])
    pd.DataFrame(
        [{k: (v if not isinstance(v, dict) else str(v)) for k, v in r.items()}
         for r in sensitivity]
    ).to_csv(w.TABLES / "wp9_sensitivity.csv", index=False)

    # ---- score the predictions --------------------------------------------
    v1 = bool((head.C3_stripped_fraction > 0.95).all())
    c1_spread = float(head.C1_age.max() - head.C1_age.min())
    c3_spread = float(head.C3_stripped_fraction.max()
                      - head.C3_stripped_fraction.min())
    v2 = bool(c1_spread > c3_spread)
    v3 = bool(outcome == "INCONCLUSIVE")

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp9_verdict.py",
        "status": "SUCCESS",
        "work_package": "WP9",
        "prereg": "provenance/wp9_verdict_prereg.json",
        "iterations": n_iter,
        "wp5_version": WP5_VERSION,
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "headline_branches": int(len(head)),
        "verdict": {
            "P_verdict_min": round(float(head.P_verdict.min()), 4),
            "P_verdict_median": round(float(head.P_verdict.median()), 4),
            "P_verdict_max": round(float(head.P_verdict.max()), 4),
            "P_verdict_permissive_median": round(
                float(head.P_verdict_permissive.median()), 4
            ),
            "outcome": outcome,
            "framing": framing,
            "rule_applied": (
                f"SUPPORTED if P >= {PLAUSIBLE_THRESHOLD} on every headline "
                f"branch, DISFAVOURED if P <= {IMPLAUSIBLE_THRESHOLD} on every "
                f"headline branch, INCONCLUSIVE otherwise.  Applied "
                f"mechanically to the numbers as computed."
            ),
        },
        "terms": {
            "C1_age": {
                "window_kyr": list(AGE_WINDOW_KYR),
                "min": round(float(head.C1_age.min()), 4),
                "median": round(float(head.C1_age.median()), 4),
                "max": round(float(head.C1_age.max()), 4),
                "permissive_median": round(
                    float(head.C1_age_permissive.median()), 4
                ),
            },
            "C3_progenitor_type": {
                "threshold_Msun": STRIPPED_PROGENITOR_MSUN,
                "min": round(float(head.C3_stripped_fraction.min()), 4),
                "median": round(float(head.C3_stripped_fraction.median()), 4),
            },
            "C4_in_situ": {
                "value": IN_SITU_FRACTION,
                "nature": "a bound from the WP7 runaway fraction, not a distribution",
            },
            "C2_energy": (
                "NOT computed and NOT multiplied in.  Haerer requires "
                "3-5e51 erg, which is above the canonical 1e51 erg and needs an "
                "energetic stripped-envelope event.  This project measures "
                "progenitor masses, not explosion energies.  C2 is reported as "
                "a conditional: the verdict below is the probability that a "
                "supernova of the right AGE, TYPE and LOCATION was available, "
                "GIVEN that such an event can reach the required energy."
            ),
        },
        "sensitivity": sensitivity,
        "dominant_driver": sensitivity[0]["axis"],
        "excluded_branches_reported_not_deleted": {
            "alpha_2p6": {
                "P_verdict_median": round(
                    float(excluded.P_verdict.median()), 4
                ),
                "why_excluded": (
                    "the WP5 calibration-window evidence alone rejects it -- "
                    "best in 1 of 18 cells, worst median chi-square.  Excluding "
                    "it does not spend the census closure."
                ),
            },
            "islands_explodability": (
                "excluded by WP8: PSR J2032+4127 is a neutron star inside the "
                "association and the islands branch predicts exactly zero "
                "supernovae.  Not run here; its verdict would be identically "
                "zero."
            ),
        },
        "predictions": [
            {
                "id": "V1", "outcome": "PASS" if v1 else "FAIL",
                "statement": "C3 > 0.95 on every headline branch",
                "min_C3": round(float(head.C3_stripped_fraction.min()), 4),
            },
            {
                "id": "V2", "outcome": "PASS" if v2 else "FAIL",
                "statement": "the verdict is driven mainly by C1, the age term",
                "C1_spread": round(c1_spread, 4),
                "C3_spread": round(c3_spread, 4),
            },
            {
                "id": "V3", "outcome": "PASS" if v3 else "FAIL",
                "statement": "the outcome is INCONCLUSIVE, so a regular article",
                "actual_outcome": outcome,
            },
        ],
        "outputs": {
            str(path.relative_to(w.ROOT)): w.sha256(path)
            for path in [
                w.TABLES / "wp9_verdict.csv",
                w.TABLES / "wp9_sensitivity.csv",
            ]
        },
    }
    w.write_json(w.PROVENANCE / "wp9_verdict_execution.json", record)

    print("\nWP9 — the verdict\n")
    print(f"  P_verdict = P(C1 age) x P(C3 type) x P(C4 in-situ)")
    print(f"  over {len(head)} headline branches "
          f"(alpha in {list(HEADLINE_ALPHAS)}, {PRIMARY_EXPLODABILITY})\n")
    t = record["terms"]
    print(f"  C1 age {AGE_WINDOW_KYR[0]:g}-{AGE_WINDOW_KYR[1]:g} kyr : "
          f"{t['C1_age']['min']:.3f} - {t['C1_age']['max']:.3f}"
          f"   (median {t['C1_age']['median']:.3f})")
    print(f"  C3 stripped progenitor     : "
          f"median {t['C3_progenitor_type']['median']:.3f}")
    print(f"  C4 in-situ bound           : {IN_SITU_FRACTION:.3f}")
    v = record["verdict"]
    print(f"\n  P_verdict: {v['P_verdict_min']:.3f} - {v['P_verdict_max']:.3f}"
          f"   median {v['P_verdict_median']:.3f}")
    print(f"  permissive window (to {AGE_PERMISSIVE_KYR[1]:g} kyr): "
          f"median {v['P_verdict_permissive_median']:.3f}")
    print(f"\n  OUTCOME: {outcome}   ->   framing: {framing}")
    print(f"\n  dominant driver: {record['dominant_driver']}")
    for entry in sensitivity:
        print(f"    {entry['axis']:16s} spread {entry['spread']:.4f}")
    print("\n  predictions:")
    for entry in record["predictions"]:
        print(f"    {entry['id']}  {entry['outcome']:4s}  {entry['statement']}")
    print("\nwrote provenance/wp9_verdict_execution.json")


if __name__ == "__main__":
    main()
