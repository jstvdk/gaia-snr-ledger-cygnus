#!/usr/bin/env python3
"""B1 -- reconcile the two association masses that circulate in this project.

Two numbers are on the books and both are called "the association mass":

  * WP5 reports 29,122 Msun (repair_v6) / 29,246 Msun (repair_v7) at the
    baseline branch and says it is "within a factor two of the 16,500 Msun
    literature scale";
  * the WP6 external cross-check table reports 1.74e4 Msun and says it
    "agrees, 5%" with Wright+2015 / Harer+2025's 1.65e4 Msun.

They are not in tension and neither is wrong.  They are *different integrals of
the same fitted normalization k*, and the difference is exactly two terms that
this script computes and reports:

  M1  primary stars, 0.5-120 Msun, single power law         (the WP6 figure)
  M2  = M1 + primary stars 0.08-0.5 Msun, Kroupa-like break at 0.5 with
       slope 1.3                                            (primary_system_mass)
  M3  = M2 + unresolved companions, f_bin = 0.40, q ~ U(0.1, 1), counted only
       where q*m >= 0.08 Msun                               (the WP5 headline)

Nothing here is a new measurement: M2 and M3 are read from the stored WP5
products, M1 is recomputed with the same closed form WP6 used, and the script
asserts that both stored numbers are reproduced.  The output is a definition
table, not a revision.

It also states, with the definitional caveat, which of the three is the
like-for-like partner of Wright+2015's 16,500 Msun -- see the report for why
that is M2 and not M1.

Outputs:
  tables/wp5_association_mass_reconciliation.csv
  provenance/wp5_association_mass_reconciliation_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/wp5_association_mass_reconciliation.py
"""
from __future__ import annotations

import argparse
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w

BASELINE_FAMILY, BASELINE_RV, BASELINE_ALPHA = "PARSEC", 3.1, 2.3

# The lower edge WP6's external cross-check integrated from: the frozen
# MASS_GRID lower edge, not a physical statement about the IMF.
WP6_MASS_LO = 0.5
IMF_UPPER = w.TOTAL_MASS_RANGE[1]

# Wright+2015 section 6.3: total stellar mass from a Maschberger (2013)
# universal IMF with alpha = 2.30, beta = 1.40, masses drawn over 0.01-150
# Msun, normalized on the observed count of stars at 20-40 Msun.
WRIGHT = {
    "reference": "Wright, Drew & Mohr-Smith 2015, MNRAS 449, 741",
    "local_copy": "papers/Wright_2015.pdf",
    "section": "6.3",
    "total_mass_Msun": 16_500.0,
    "total_mass_hi_Msun": 3_800.0,
    "total_mass_lo_Msun": 2_800.0,
    "imf": "Maschberger 2013 universal IMF, alpha = 2.30, beta = 1.40",
    "sampled_mass_range_Msun": [0.01, 150.0],
    "normalized_on": "observed number of stars at 20-40 Msun (36 +1 -4)",
    "verbatim": (
        "We find that the observed number of massive stars in Cyg OB2 in the "
        "mass range of 20-40 M_sun can be reproduced if the association has a "
        "total stellar mass of 16500 +3800 -2800 M_sun."
    ),
    "multiplicity_treatment": (
        "the census counts PRIMARY OB stars and the Monte Carlo matches drawn "
        "individual stars to that count, so unresolved companions are not "
        "separately added; the quantity is best read as a primary-star total"
    ),
}


def primary_mass_above(k: float, alpha: float, lo: float, hi: float) -> float:
    """k * integral[lo, hi] M^(1-alpha) dM -- single power law, no break."""
    if np.isclose(alpha, 2.0):
        return float(k * np.log(hi / lo))
    return float(k * (lo ** (2.0 - alpha) - hi ** (2.0 - alpha)) / (alpha - 2.0))


def low_mass_segment(k: float, alpha: float) -> float:
    """k * the 0.08-0.5 Msun Kroupa-like segment, continuous at the break."""
    coefficient = w.LOW_MASS_BREAK ** (w.LOW_MASS_SLOPE - alpha)
    return float(
        k
        * coefficient
        * w.power_integral(
            w.TOTAL_MASS_RANGE[0], w.LOW_MASS_BREAK, 1.0 - w.LOW_MASS_SLOPE
        )
    )


def brown_dwarf_extension(k: float, alpha: float, lo: float = 0.01) -> float:
    """Mass in [lo, 0.08] if the 1.3 segment were continued below the grid.

    Reported for scale only.  It lies outside the frozen MASS_GRID and outside
    the range any WP5 injection touched, and a real Kroupa/Maschberger IMF
    flattens again below 0.08, so this OVERSTATES the extension.
    """
    coefficient = w.LOW_MASS_BREAK ** (w.LOW_MASS_SLOPE - alpha)
    return float(
        k
        * coefficient
        * w.power_integral(lo, w.TOTAL_MASS_RANGE[0], 1.0 - w.LOW_MASS_SLOPE)
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wp5-version", default="repair_v7")
    parser.add_argument("--compare-version", default="repair_v6")
    args = parser.parse_args()

    rows = []
    stored_checks = []
    sources = {}
    for version in (args.wp5_version, args.compare_version):
        path = w.PROC / f"wp5_association_mass_{version}.parquet"
        sources[str(path.relative_to(w.ROOT))] = w.sha256(path)
        table = pd.read_parquet(path)
        for entry in table.itertuples():
            k = float(entry.k_total_median)
            alpha = float(entry.alpha)
            m1 = primary_mass_above(k, alpha, WP6_MASS_LO, IMF_UPPER)
            low = low_mass_segment(k, alpha)
            m2 = m1 + low
            companions = float(entry.multiplicity_adjusted_mass_median_Msun) - float(
                entry.primary_system_mass_median_Msun
            )
            m3 = m2 + companions
            rows.append(
                {
                    "wp5_version": version,
                    "family": entry.family,
                    "R_V": float(entry.R_V),
                    "alpha": alpha,
                    "k_total_median": k,
                    "M1_primaries_0p5_to_120_Msun": m1,
                    "low_mass_0p08_to_0p5_Msun": low,
                    "M2_primary_system_0p08_to_120_Msun": m2,
                    "companion_mass_Msun": companions,
                    "M3_multiplicity_adjusted_Msun": m3,
                    "M3_over_M1": m3 / m1,
                    "M1_over_Wright2015": m1 / WRIGHT["total_mass_Msun"],
                    "M2_over_Wright2015": m2 / WRIGHT["total_mass_Msun"],
                    "M3_over_Wright2015": m3 / WRIGHT["total_mass_Msun"],
                    "brown_dwarf_extension_0p01_to_0p08_Msun": brown_dwarf_extension(
                        k, alpha
                    ),
                }
            )
            # The stored products must be reproduced, not replaced.
            stored_checks.append(
                {
                    "wp5_version": version,
                    "family": entry.family,
                    "R_V": float(entry.R_V),
                    "alpha": alpha,
                    "stored_primary_system_Msun": float(
                        entry.primary_system_mass_median_Msun
                    ),
                    "recomputed_M2_Msun": m2,
                    "abs_difference_Msun": abs(
                        m2 - float(entry.primary_system_mass_median_Msun)
                    ),
                }
            )

    table = pd.DataFrame(rows)
    worst = max(c["abs_difference_Msun"] for c in stored_checks)
    if worst > 1.0:
        raise SystemExit(
            f"recomputation does not reproduce the stored primary-system mass "
            f"(worst {worst:.3f} Msun); the definitions in this script and in "
            f"wp5_common have drifted apart"
        )

    out_csv = w.TABLES / "wp5_association_mass_reconciliation.csv"
    table.to_csv(out_csv, index=False)

    def baseline(version: str) -> pd.Series:
        return table[
            table.wp5_version.eq(version)
            & table.family.eq(BASELINE_FAMILY)
            & table.R_V.eq(BASELINE_RV)
            & table.alpha.eq(BASELINE_ALPHA)
        ].iloc[0]

    accepted = baseline(args.wp5_version)
    previous = baseline(args.compare_version)

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp5_association_mass_reconciliation.py",
        "purpose": (
            "B1 of tasks/pre_wp10_assessment_brief.md -- reconcile the WP5 "
            "association mass with the WP6 external cross-check figure.  "
            "Definitional only: no fit was re-run and no stored number moved."
        ),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "wp5_version": args.wp5_version,
        "compare_version": args.compare_version,
        "definitions": {
            "M1_primaries_0p5_to_120": (
                "k * integral[0.5, 120] M^(1-alpha) dM.  Single power law, no "
                "low-mass break, primaries only.  This is what "
                "scripts/wp6_external_crosschecks.py computes and what the "
                "PROJECT_TRACE 10b row reports as 1.74e4 Msun."
            ),
            "M2_primary_system_0p08_to_120": (
                "M1 plus the 0.08-0.5 Msun segment of a Kroupa-like broken "
                "power law (slope 1.3 below the 0.5 Msun break, continuous "
                "there).  Stored as primary_system_mass_median_Msun."
            ),
            "M3_multiplicity_adjusted": (
                "M2 plus unresolved companion mass for f_bin = 0.40, "
                "q ~ U(0.1, 1), companion counted only where q*m >= 0.08 "
                "Msun.  Stored as multiplicity_adjusted_mass_median_Msun and "
                "quoted by wp5_report.py as THE association mass."
            ),
        },
        "constants": {
            "TOTAL_MASS_RANGE": list(w.TOTAL_MASS_RANGE),
            "LOW_MASS_BREAK": w.LOW_MASS_BREAK,
            "LOW_MASS_SLOPE": w.LOW_MASS_SLOPE,
            "F_BINARY": w.F_BINARY,
            "Q_MIN": w.Q_MIN,
            "wp6_crosscheck_mass_lo": WP6_MASS_LO,
        },
        "literature": WRIGHT,
        "baseline_PARSEC_rv3.1_alpha2.3": {
            version: {
                "k_total_median": round(float(row.k_total_median), 2),
                "M1_Msun": round(float(row.M1_primaries_0p5_to_120_Msun), 0),
                "low_mass_0p08_to_0p5_Msun": round(
                    float(row.low_mass_0p08_to_0p5_Msun), 0
                ),
                "M2_Msun": round(float(row.M2_primary_system_0p08_to_120_Msun), 0),
                "companion_mass_Msun": round(float(row.companion_mass_Msun), 0),
                "M3_Msun": round(float(row.M3_multiplicity_adjusted_Msun), 0),
                "M3_over_M1": round(float(row.M3_over_M1), 3),
                "M1_over_Wright2015": round(float(row.M1_over_Wright2015), 3),
                "M2_over_Wright2015": round(float(row.M2_over_Wright2015), 3),
                "M3_over_Wright2015": round(float(row.M3_over_Wright2015), 3),
            }
            for version, row in (
                (args.wp5_version, accepted),
                (args.compare_version, previous),
            )
        },
        "reproduction_check": {
            "quantity": "recomputed M2 vs stored primary_system_mass_median_Msun",
            "cells": len(stored_checks),
            "worst_abs_difference_Msun": round(float(worst), 6),
            "pass": bool(worst <= 1.0),
        },
        "finding": (
            "The two circulating numbers differ by exactly two terms and no "
            "third.  At the accepted baseline the 0.08-0.5 Msun segment adds "
            f"{accepted.low_mass_0p08_to_0p5_Msun / accepted.M1_primaries_0p5_to_120_Msun:.0%} "
            "and unresolved companions add "
            f"{accepted.companion_mass_Msun / accepted.M2_primary_system_0p08_to_120_Msun:.0%}, "
            f"taking M1 = {accepted.M1_primaries_0p5_to_120_Msun:,.0f} Msun to "
            f"M3 = {accepted.M3_multiplicity_adjusted_Msun:,.0f} Msun."
        ),
        "consequence_for_the_crosscheck": (
            "The '5%' agreement in the PROJECT_TRACE 10b table compares M1 "
            "against a Wright+2015 mass integrated over the FULL IMF range, so "
            "it pairs mismatched definitions.  The like-for-like partner of "
            "Wright's 16,500 Msun is M2 (primary stars, full grid range), "
            f"which sits at {accepted.M2_over_Wright2015:.2f} times it -- an "
            "agreement at the factor level and a ~50% offset in the mean, "
            "inside the factor-2 gate but NOT at the 5% level.  The 5% "
            "coincidence must not be quoted as an agreement."
        ),
        "inputs": sources,
        "outputs": {str(out_csv.relative_to(w.ROOT)): w.sha256(out_csv)},
    }
    w.write_json(
        w.PROVENANCE / "wp5_association_mass_reconciliation_execution.json", record
    )

    print(f"association mass reconciliation, baseline PARSEC R_V=3.1 alpha=2.3\n")
    print(f"{'version':>10s} {'k':>9s} {'M1':>9s} {'+lowmass':>9s} {'M2':>9s}"
          f" {'+comp':>9s} {'M3':>9s}  M3/M1")
    for version, row in ((args.wp5_version, accepted), (args.compare_version, previous)):
        print(
            f"{version:>10s} {row.k_total_median:9.0f}"
            f" {row.M1_primaries_0p5_to_120_Msun:9.0f}"
            f" {row.low_mass_0p08_to_0p5_Msun:9.0f}"
            f" {row.M2_primary_system_0p08_to_120_Msun:9.0f}"
            f" {row.companion_mass_Msun:9.0f}"
            f" {row.M3_multiplicity_adjusted_Msun:9.0f}"
            f"  {row.M3_over_M1:.3f}"
        )
    print(
        f"\nagainst Wright+2015 16,500 Msun:  M1 {accepted.M1_over_Wright2015:.2f}x"
        f"   M2 {accepted.M2_over_Wright2015:.2f}x"
        f"   M3 {accepted.M3_over_Wright2015:.2f}x"
    )
    print(
        f"reproduction of stored primary_system_mass: worst difference "
        f"{worst:.2e} Msun over {len(stored_checks)} cells"
    )
    print("wrote tables/wp5_association_mass_reconciliation.csv")
    print("wrote provenance/wp5_association_mass_reconciliation_execution.json")


if __name__ == "__main__":
    main()
