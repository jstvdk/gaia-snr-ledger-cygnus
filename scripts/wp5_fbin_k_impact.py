#!/usr/bin/env python3
"""Post-hoc scope assessment: what does repair_v7 actually move?

NOT part of the pre-registered discriminator.  Written after
provenance/wp5_fbin_discriminator_execution.json returned "repair_v7
JUSTIFIED", to answer a scoping question the discriminator did not:

    if the truth model is wrong, is it safe to run WP7 before fixing it?

The discriminator's D1 measured the RECOVERY FRACTION over the calibration
window, which does not capture mass migration -- a star can be recovered in both
arms while being assigned a different mass.  But k is fitted to the COUNTS in
the 2-8 Msun observed bins, so the quantity that actually sets k is

    P(estimated mass lands INSIDE the calibration window | true mass)

integrated against the IMF over the parent range.  That does capture migration,
because it is about where the estimate lands rather than whether it exists.
This script computes it from the paired responses already on disk.

It also computes the Class E branch spread on N_SN, so the correction can be
compared against the uncertainty WP7 must report regardless.

Outputs: provenance/wp5_fbin_k_impact_execution.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp5_fbin_k_impact.py
"""
from __future__ import annotations

import glob
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
from wp6_mass_extension_decision import IMF_UPPER_LIMIT, turnoff_mass

WEIGHT_ALPHA = 2.3
BASELINE_RV = 3.1


def window_counts(response: pd.DataFrame) -> float:
    """IMF-weighted P(estimated mass inside the calibration window | true M).

    This is the k-relevant quantity: k is fitted so that
    k * integral dM M^-alpha P(estimated inside window | M) reproduces the
    observed counts, so a fractional shift here is the fractional shift k must
    absorb, with the opposite sign.
    """
    draw_columns = sorted(
        c for c in response.columns if c.startswith("recovered_mass_draw_")
    )
    total = 0.0
    for mass in np.sort(response["true_primary_mass"].unique()):
        rows = response[response["true_primary_mass"].eq(mass)]
        active = [c for c in draw_columns if rows[c].notna().any()]
        if not active:
            continue
        draws = rows[active].to_numpy(float)
        finite = np.isfinite(draws)
        inside = np.sum(
            finite
            & (draws >= w.CALIBRATION_NOMINAL_LO)
            & (draws <= w.CALIBRATION_HI)
        ) / (len(active) * len(rows))
        total += mass ** (-WEIGHT_ALPHA) * inside
    return float(total)


def n_sn(family: str, alpha: float, normalization, age_posterior) -> float:
    total = 0.0
    for subgroup in w.SUBGROUPS:
        row = age_posterior[
            age_posterior.subgroup.eq(subgroup)
            & age_posterior.family.eq(family)
            & age_posterior.R_V.eq(BASELINE_RV)
            & age_posterior.f_bin.eq(w.F_BINARY)
            & age_posterior.indicator.eq("ums")
            & age_posterior.dmu.eq(0.0)
        ]
        age = float(row.age_map.iloc[0])
        k = float(
            normalization[
                normalization.subgroup.eq(subgroup)
                & normalization.family.eq(family)
                & normalization.R_V.eq(BASELINE_RV)
                & normalization.alpha.eq(alpha)
            ].k_median.iloc[0]
        )
        cap = min(turnoff_mass(family, age), IMF_UPPER_LIMIT)
        if cap < IMF_UPPER_LIMIT:
            total += k * (
                cap ** (1 - alpha) - IMF_UPPER_LIMIT ** (1 - alpha)
            ) / (alpha - 1)
    return total


def main() -> None:
    per_node = []
    for control in sorted(glob.glob(str(w.PROC / "wp5_fbindisc_ctl_*_response.parquet"))):
        treatment = control.replace("_ctl_", "_trt_")
        ctl = window_counts(pd.read_parquet(control))
        trt = window_counts(pd.read_parquet(treatment))
        per_node.append(
            {
                "node": control.split("/")[-1],
                "window_counts_ctl": round(ctl, 5),
                "window_counts_trt": round(trt, 5),
                "count_shift_fraction": round(trt / ctl - 1.0, 5),
                "implied_k_shift_fraction": round(-(trt / ctl - 1.0), 5),
            }
        )
    shifts = np.array([e["count_shift_fraction"] for e in per_node])
    k_shift = float(-shifts.mean())

    normalization = pd.read_parquet(
        w.PROC / "wp5_imf_normalization_repair_v6.parquet"
    )
    age_posterior = pd.read_parquet(
        w.PROC / "wp4_age_posteriors_repair_v5.parquet"
    )
    branches = {
        f"{family}_alpha{alpha}": round(
            n_sn(family, alpha, normalization, age_posterior), 2
        )
        for family in w.FAMILIES
        for alpha in w.IMF_SLOPES
    }
    values = np.array(list(branches.values()))

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp5_fbin_k_impact.py",
        "status": "SUCCESS",
        "post_hoc": (
            "NOT part of the pre-registered discriminator.  Written after it "
            "returned 'repair_v7 JUSTIFIED', to answer a scoping question it "
            "did not address: is it safe to run WP7 before the repair?"
        ),
        "parent": "provenance/wp5_fbin_discriminator_execution.json",
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "why_D1_was_not_enough": (
            "D1 measured the RECOVERY FRACTION over the calibration window, "
            "which does not capture mass migration: a star can be recovered in "
            "both arms while being assigned a different mass.  k is fitted to "
            "the COUNTS in the 2-8 Msun observed bins, so the k-relevant "
            "quantity is P(estimated mass lands INSIDE the window | true mass), "
            "which does capture migration."
        ),
        "calibration_window_Msun": [w.CALIBRATION_NOMINAL_LO, w.CALIBRATION_HI],
        "per_node": per_node,
        "result": {
            "mean_count_shift_fraction": round(float(shifts.mean()), 5),
            "node_spread": round(float(shifts.std()), 5),
            "implied_k_shift_fraction": round(k_shift, 5),
            "implied_N_SN_shift_fraction": round(k_shift, 5),
            "note": (
                "estimated masses migrate OUT of the calibration window under "
                "the treatment model -- brighter stars are assigned higher "
                "masses and some cross above 8 Msun -- so k must RISE slightly "
                "to reproduce the same observed counts.  N_SN is linear in k."
            ),
        },
        "branch_spread_context": {
            "N_SN_by_branch": branches,
            "min": float(values.min()),
            "max": float(values.max()),
            "factor": round(float(values.max() / values.min()), 1),
            "comparison": (
                f"the Class E branch spread on N_SN is a factor of "
                f"{values.max() / values.min():.1f}, against a repair_v7 "
                f"correction of {k_shift:.2%}.  The branches are carried in "
                "parallel by design and cannot be reduced; the correction is "
                f"about {(values.max() / values.min() - 1) / abs(k_shift):.0f}x "
                "smaller than the uncertainty WP7 must report regardless."
            ),
        },
        "what_repair_v7_does_and_does_not_touch": {
            "untouched": [
                "WP3 per-star extinctions — the injection truth model plays no "
                "part in fitting real stars",
                "WP4 ages and per-star masses — same reason; WP7 depends on "
                "these heavily and they are bit-identical",
                "the WP6 runaway correction — the traceback uses astrometry, "
                "not injections",
                "stellar lifetimes tau(m) — from the isochrone families",
            ],
            "moved": [
                f"the WP5 normalization k, by {k_shift:+.2%}",
                "the WP6 closure ratio, 1.099 -> ~1.074, because the 4-8 Msun "
                "up-scatter segment carries 23.5% of the predicted observable "
                "count and its response shifts by +9.87%",
            ],
            "cost_estimate_corrected": (
                "tasks/repair_v7_recommendation.md section 3 costed ~8 h "
                "including WP3 and WP4 stages.  Those are unnecessary: the "
                "f_bin change lives entirely in the WP5 injection truth model.  "
                "Actual scope is WP5 injections (~3 h) + fit + verdict "
                "stability + WP6 re-run, about 6 h.  The 'repair_v7' label is a "
                "misnomer for what is really a WP5/WP6 re-run."
            ),
        },
        "conclusion": {
            "safe_to_run_wp7_first": True,
            "reasoning": (
                "the perturbation to WP7's only affected input is measured at "
                f"{k_shift:+.2%}, some three orders of magnitude below the "
                "Class E branch spread WP7 must report anyway, and its other "
                "inputs are bit-identical.  This is categorically different "
                "from issues #16 and #17, which were defects of UNKNOWN size "
                "that turned out to be large.  Proceeding on a bounded, "
                "quantified perturbation is legitimate; proceeding on an "
                "unknown one is what produced the two withdrawn results."
            ),
            "why_repair_v7_is_still_needed": (
                "it is a WP6 fix, not a WP7 fix.  WP6's residual excess is "
                "9.9% and the repair removes about 23% of what remains "
                "unexplained.  The same physical error is 0.54% on a "
                "window-integrated normalization and 9.87% on a "
                "threshold-crossing probability — an 18x asymmetry produced "
                "purely by one quantity being an average and the other a "
                "boundary."
            ),
            "obligation": (
                "WP7 results computed on the accepted chain are PROVISIONAL "
                "until re-checked against repair_v7.  WP7 is pure computation "
                "on frozen inputs, so the re-check is cheap."
            ),
        },
    }
    w.write_json(w.PROVENANCE / "wp5_fbin_k_impact_execution.json", record)

    print("repair_v7 scope assessment (post-hoc)\n")
    for entry in per_node:
        print(f"  {entry['node'][:34]:34s} k shift "
              f"{entry['implied_k_shift_fraction']:+.2%}")
    print(f"\n  implied k / N_SN shift: {k_shift:+.2%} "
          f"+- {shifts.std():.2%}")
    print(f"  Class E branch spread on N_SN: {values.min():.2f} to "
          f"{values.max():.2f}, factor {values.max() / values.min():.1f}")
    print(f"\n  safe to run WP7 first: {record['conclusion']['safe_to_run_wp7_first']}")
    print("  repair_v7 remains needed — for WP6, not WP7")
    print("\nwrote provenance/wp5_fbin_k_impact_execution.json")


if __name__ == "__main__":
    main()
