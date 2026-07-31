#!/usr/bin/env python3
"""B2 -- reconcile the two subgroup ages that circulate in this project.

Two age sets are on the books:

  WP4 upper-MS MAP (repair_v5)        A 3.981   B 3.548   C 2.512 Myr
  obligation O1 in PROJECT_TRACE      A 4.00    B 4.07    C 2.52  Myr

The second set is not stale and is not a typo: it is the WP5 **fitted truth-age
posterior mean** on the baseline branch, the age at which the injection truth
model reproduces the observed 2-8 Msun mass function, with the WP4 posterior as
prior and the node weights updated by the Poisson likelihood of the counts.
This script demonstrates that identification, measures how far the two differ
per subgroup, and -- because N_SN is essentially the IMF integral above the
turnoff -- converts the difference into the quantity the paper is about.

It also measures something the age tables do not show: how much of CygOB2-B's
fitted posterior sits on the TOP node of its own WP4 support.  A posterior that
piles up against the edge of its prior grid is reporting a bound, not a
measurement, and the direction of that bound matters (older -> lower turnoff ->
more supernovae).

Nothing is refitted.  Every number is read from stored products.

Outputs:
  tables/wp4_wp5_age_reconciliation.csv
  provenance/wp4_wp5_age_reconciliation_execution.json

Run:
  PYTHONPATH=scripts python3 scripts/wp4_wp5_age_reconciliation.py
"""
from __future__ import annotations

import argparse
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
from wp6_mass_extension_decision import IMF_UPPER_LIMIT, turnoff_mass

BASELINE_FAMILY, BASELINE_RV, BASELINE_ALPHA = "PARSEC", 3.1, 2.3


def dead_fraction_per_k(turnoff: float, alpha: float) -> float:
    """integral[turnoff, 120] M^-alpha dM -- the SN count per unit k.

    Capped at the IMF ceiling: a turnoff above 120 Msun means nothing has died.
    """
    if turnoff >= IMF_UPPER_LIMIT:
        return 0.0
    return float(
        (turnoff ** (1.0 - alpha) - IMF_UPPER_LIMIT ** (1.0 - alpha)) / (alpha - 1.0)
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--wp5-version", default="repair_v7")
    parser.add_argument("--wp4-version", default="repair_v5")
    args = parser.parse_args()

    wp4_path = w.PROC / f"wp4_age_posteriors_{args.wp4_version}.parquet"
    wp5_path = w.PROC / f"wp5_imf_normalization_{args.wp5_version}.parquet"
    draws_path = w.PROC / f"wp5_imf_posterior_draws_{args.wp5_version}.npz"
    wp4 = pd.read_parquet(wp4_path)
    wp5 = pd.read_parquet(wp5_path)
    draws = np.load(draws_path)

    rows = []
    for entry in wp5.itertuples():
        subgroup, family, rv, alpha = (
            entry.subgroup, entry.family, float(entry.R_V), float(entry.alpha)
        )
        # The WP4 branch that fed this cell: same family and R_V, the baseline
        # binary fraction, the upper-MS indicator, no distance-modulus offset.
        prior_row = wp4[
            wp4.subgroup.eq(subgroup)
            & wp4.family.eq(family)
            & wp4.R_V.eq(rv)
            & wp4.f_bin.eq(w.F_BINARY)
            & wp4.indicator.eq("ums")
            & wp4.dmu.eq(0.0)
        ]
        if prior_row.empty:
            continue
        ums_map = float(prior_row.age_map.iloc[0])
        ums_lo, ums_hi = (
            float(prior_row.age_lo68.iloc[0]),
            float(prior_row.age_hi68.iloc[0]),
        )

        nodes = np.asarray(entry.truth_age_nodes_Myr, dtype=float)
        posterior = np.asarray(entry.truth_age_posterior_weights, dtype=float)
        fitted = float(entry.truth_age_posterior_mean_Myr)

        to_ums = turnoff_mass(family, ums_map)
        to_fit = turnoff_mass(family, fitted)
        dead_ums = dead_fraction_per_k(to_ums, alpha)
        dead_fit = dead_fraction_per_k(to_fit, alpha)

        rows.append(
            {
                "subgroup": subgroup,
                "family": family,
                "R_V": rv,
                "alpha": alpha,
                "wp4_ums_map_Myr": ums_map,
                "wp4_ums_lo68_Myr": ums_lo,
                "wp4_ums_hi68_Myr": ums_hi,
                "wp5_fitted_posterior_mean_Myr": fitted,
                "wp5_fitted_posterior_map_Myr": float(
                    entry.truth_age_posterior_map_Myr
                ),
                "age_shift_Myr": fitted - ums_map,
                "top_node_Myr": float(nodes[-1]),
                "top_node_posterior_weight": float(posterior[-1]),
                "top_two_nodes_posterior_weight": float(posterior[-2:].sum()),
                "turnoff_at_ums_map_Msun": to_ums,
                "turnoff_at_fitted_Msun": to_fit,
                "sn_per_k_at_ums_map": dead_ums,
                "sn_per_k_at_fitted": dead_fit,
                "sn_ratio_fitted_over_ums": (
                    dead_fit / dead_ums if dead_ums > 0 else np.nan
                ),
            }
        )

    table = pd.DataFrame(rows)
    out_csv = w.TABLES / "wp4_wp5_age_reconciliation.csv"
    table.to_csv(out_csv, index=False)

    baseline = table[
        table.family.eq(BASELINE_FAMILY)
        & table.R_V.eq(BASELINE_RV)
        & table.alpha.eq(BASELINE_ALPHA)
    ].set_index("subgroup")

    # The railing bound: what B's supernova count would be if its true age sat
    # at the top of the node support rather than at the fitted mean.
    b = baseline.loc["CygOB2-B"]
    dead_top = dead_fraction_per_k(
        turnoff_mass(BASELINE_FAMILY, b.top_node_Myr), BASELINE_ALPHA
    )
    railing_headroom = dead_top / b.sn_per_k_at_fitted - 1.0

    # Draw-level check: the ledger consumes these draws, not the summary means.
    draw_means = {}
    for subgroup in ("CygOB2-A", "CygOB2-B", "CygOB2-C"):
        key = (
            f"truth_age_draws__{subgroup.replace('-', '_')}__{BASELINE_FAMILY}"
            f"__rv{BASELINE_RV:.1f}".replace(".", "p")
            + f"__a{BASELINE_ALPHA:.1f}".replace(".", "p")
        )
        draw_means[subgroup] = float(draws[key].mean())

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp4_wp5_age_reconciliation.py",
        "purpose": (
            "B2 of tasks/pre_wp10_assessment_brief.md -- identify which age is "
            "which, and price the difference in supernovae.  Definitional and "
            "diagnostic only: nothing was refitted, no stored number moved."
        ),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "wp4_version": args.wp4_version,
        "wp5_version": args.wp5_version,
        "definitions": {
            "wp4_ums_map": (
                "MAP of the WP4 upper-main-sequence CMD age posterior at "
                "f_bin = 0.40, dmu = 0.  A photometric age from isochrone "
                "fitting to the observed CMD."
            ),
            "wp5_fitted_posterior_mean": (
                "mean of the WP5 joint age-k truth-age posterior: the WP4 "
                "posterior taken as prior over nine unsnapped quantile nodes, "
                "reweighted by the Poisson likelihood of the observed 2-8 Msun "
                "counts with k integrated out under a Jeffreys prior.  A "
                "counts-based age."
            ),
        },
        "identification": {
            "claim": (
                "obligation O1's 'A 4.00, B 4.07, C 2.52' is the WP5 fitted "
                "truth-age posterior mean on the baseline branch, not a stale "
                "or alternative WP4 number"
            ),
            "baseline_fitted_posterior_mean_Myr": {
                s: round(float(baseline.loc[s].wp5_fitted_posterior_mean_Myr), 3)
                for s in baseline.index
            },
            "baseline_wp4_ums_map_Myr": {
                s: round(float(baseline.loc[s].wp4_ums_map_Myr), 3)
                for s in baseline.index
            },
            "ledger_draw_means_Myr": {
                s: round(v, 3) for s, v in draw_means.items()
            },
            "consumed_downstream_by": (
                "wp7_ledger.py and wp9_verdict.py both read "
                "truth_age_draws__* from the WP5 posterior-draw archive, so "
                "the counts-based age is the one every downstream number "
                "rests on; the upper-MS MAP enters only as its prior"
            ),
        },
        "baseline_PARSEC_rv3.1_alpha2.3": [
            {
                "subgroup": s,
                "wp4_ums_map_Myr": round(float(r.wp4_ums_map_Myr), 3),
                "wp4_ums_68_Myr": [
                    round(float(r.wp4_ums_lo68_Myr), 3),
                    round(float(r.wp4_ums_hi68_Myr), 3),
                ],
                "wp5_fitted_mean_Myr": round(
                    float(r.wp5_fitted_posterior_mean_Myr), 3
                ),
                "shift_Myr": round(float(r.age_shift_Myr), 3),
                "turnoff_ums_Msun": round(float(r.turnoff_at_ums_map_Msun), 1),
                "turnoff_fitted_Msun": round(float(r.turnoff_at_fitted_Msun), 1),
                "sn_ratio_fitted_over_ums": (
                    None
                    if not np.isfinite(r.sn_ratio_fitted_over_ums)
                    else round(float(r.sn_ratio_fitted_over_ums), 3)
                ),
                "top_node_Myr": round(float(r.top_node_Myr), 3),
                "top_node_posterior_weight": round(
                    float(r.top_node_posterior_weight), 3
                ),
            }
            for s, r in baseline.iterrows()
        ],
        "railing_finding": {
            "subgroup": "CygOB2-B",
            "top_node_Myr": round(float(b.top_node_Myr), 3),
            "top_node_posterior_weight": round(float(b.top_node_posterior_weight), 3),
            "top_two_nodes_posterior_weight": round(
                float(b.top_two_nodes_posterior_weight), 3
            ),
            "grid_median_top_node_weight_B": round(
                float(
                    table[table.subgroup.eq("CygOB2-B")].top_node_posterior_weight.median()
                ),
                3,
            ),
            "grid_median_top_node_weight_A": round(
                float(
                    table[table.subgroup.eq("CygOB2-A")].top_node_posterior_weight.median()
                ),
                3,
            ),
            "grid_median_top_node_weight_C": round(
                float(
                    table[table.subgroup.eq("CygOB2-C")].top_node_posterior_weight.median()
                ),
                3,
            ),
            "sn_headroom_if_true_age_at_top_node": round(float(railing_headroom), 3),
            "reading": (
                "B's counts-based age piles up against the top of its own WP4 "
                "node support, so the fitted age is truncated by the prior "
                "grid, not by the data.  The bound is one-sided and its "
                "direction is known: older means a lower turnoff means MORE "
                "supernovae, so B's contribution to N_SN is a lower bound at "
                "the level quoted here."
            ),
        },
        "inputs": {
            str(p.relative_to(w.ROOT)): w.sha256(p)
            for p in (wp4_path, wp5_path, draws_path)
        },
        "outputs": {str(out_csv.relative_to(w.ROOT)): w.sha256(out_csv)},
    }
    w.write_json(w.PROVENANCE / "wp4_wp5_age_reconciliation_execution.json", record)

    print("B2 -- subgroup ages, baseline PARSEC R_V = 3.1 alpha = 2.3\n")
    print(
        f"{'subgroup':<12s} {'UMS MAP':>8s} {'fitted':>8s} {'shift':>7s}"
        f" {'M_to(UMS)':>10s} {'M_to(fit)':>10s} {'SN ratio':>9s}"
        f" {'top-node w':>11s}"
    )
    for s, r in baseline.iterrows():
        ratio = (
            "   n/a" if not np.isfinite(r.sn_ratio_fitted_over_ums)
            else f"{r.sn_ratio_fitted_over_ums:9.2f}"
        )
        print(
            f"{s:<12s} {r.wp4_ums_map_Myr:8.3f} {r.wp5_fitted_posterior_mean_Myr:8.3f}"
            f" {r.age_shift_Myr:+7.3f} {r.turnoff_at_ums_map_Msun:10.1f}"
            f" {r.turnoff_at_fitted_Msun:10.1f} {ratio}"
            f" {r.top_node_posterior_weight:11.3f}"
        )
    print(
        f"\nCygOB2-B rails: {b.top_two_nodes_posterior_weight:.1%} of its fitted "
        f"posterior sits on the top two of nine nodes.  If its true age were the "
        f"top node ({b.top_node_Myr:.2f} Myr) its supernova count would be "
        f"{railing_headroom:+.1%}."
    )
    print("wrote tables/wp4_wp5_age_reconciliation.csv")
    print("wrote provenance/wp4_wp5_age_reconciliation_execution.json")


if __name__ == "__main__":
    main()
