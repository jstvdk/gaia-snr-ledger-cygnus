#!/usr/bin/env python3
"""WP6 step 0a: record the decision to extend the injection response above 18 Msun.

Written BEFORE the extension injections are run, so the design and its
predicted downstream effects cannot be reconstructed after seeing the closure
ratio -- which is the load-bearing number of the paper and exactly the kind of
quantity where knowing the answer first is corrosive.

THE PROBLEM.  WP5's injection response was measured on
``wp5_common.MASS_GRID``, which stops at PARENT_MASS_HI = 18 Msun.  That
ceiling was derived and correct for WP5: the >8 Msun contribution to the
top calibration bin is 100.0% recovered at 18 Msun for every alpha.  But WP6's
closure test integrates from 8 Msun to the isochrone turnoff, which is 64 Msun
(PARSEC) or 72 Msun (MIST) at CygOB2-A and B's age and far higher at C's.
Observed member masses already reach 61 Msun.  The predicted and observed sides
of the closure test therefore cover different mass ranges and are not
comparable as things stand.

THE DECISION (project owner, 2026-07-28).  Extend the injections rather than
assume the response stays flat above 18 Msun.  The closure ratio is the number
the supernova count rests on, so it is measured, not extrapolated.

WHAT IS DELIBERATELY NOT CHANGED.  ``wp5_common.MASS_GRID`` and
``PARENT_MASS_HI`` are frozen.  The accepted repair_v6 fit filters its response
to that grid (wp5_joint_age_fit._category_counts), so widening it globally
would silently change an accepted result on the next re-run.  The extension is
injected on its own grid, written to its own files, and consumed only by WP6.

Output: provenance/wp6_mass_extension_decision.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp6_mass_extension_decision.py
"""
from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
import wp5_joint_age_fit as J

# Log-spaced continuation of MASS_GRID from its 18 Msun ceiling to the physical
# IMF limit.  Spacing coarsens upward because the IMF integrand there is small
# and smooth and the response integral uses non-uniform trapezoid weights, the
# same argument MASS_GRID itself already uses above 8 Msun.
WP6_MASS_EXTENSION = np.array(
    [20.0, 23.0, 27.0, 31.0, 36.0, 42.0, 48.0, 56.0, 65.0, 75.0, 87.0, 100.0, 115.0]
)
IMF_UPPER_LIMIT = 120.0


_TURNOFF_CACHE: dict[str, tuple[np.ndarray, np.ndarray]] = {}


def turnoff_sequence(family: str) -> tuple[np.ndarray, np.ndarray]:
    """Native (age, most massive star still on the PMS/MS locus) pairs.

    Enforced monotone non-increasing in age.  Both grids have one inversion
    from table sampling at the very top of the locus -- PARSEC 331.4 -> 336.0
    Msun between 1.995 and 2.239 Myr, MIST 69.0 -> 74.9 between 3.181 and 3.571
    -- and a turnoff that rises with age is unphysical, so a running minimum is
    applied and the correction is recorded.
    """
    if family in _TURNOFF_CACHE:
        return _TURNOFF_CACHE[family]
    frame = pd.read_parquet(w.PROC / f"wp3_isochrones_{family.lower()}.parquet")
    frame = frame[frame["label"] <= 1] if family == "PARSEC" else frame[frame["phase"] <= 1]
    ages = np.sort(frame["age_Myr"].unique())
    raw = np.array(
        [float(frame[np.isclose(frame["age_Myr"], age)]["Mini"].max()) for age in ages]
    )
    monotone = np.minimum.accumulate(raw)
    _TURNOFF_CACHE[family] = (ages, monotone)
    return ages, monotone


def turnoff_mass(family: str, age_myr: float) -> float:
    """Most massive star still alive at this age, continuous in age.

    NOT ``load_isochrone_between_ages(...)["Mini"].max()``.  That loader takes
    the *intersection* of the two bracketing tables' mass ranges, so its ceiling
    is always the older bracket's turnoff and is therefore a STEP function of
    age -- 48.0 Msun at both 4.00 and 4.20 Myr, 69.0 at both 3.00 and 3.50.
    Harmless for WP5, where every injected mass is far below the turnoff, but
    not for WP6: N_SN is essentially the IMF integral above this value, so a
    stepped turnoff makes N_SN a step function of age.  That is the same defect
    class as issue #13.

    Instead the turnoff is interpolated log-log between native tables, which is
    continuous and uses only tabulated values.
    """
    ages, masses = turnoff_sequence(family)
    return float(
        10.0
        ** np.interp(
            np.log10(float(age_myr)), np.log10(ages), np.log10(masses)
        )
    )


def main() -> None:
    age_posterior = pd.read_parquet(w.PROC / "wp4_age_posteriors_repair_v5.parquet")
    native = {family: J.native_isochrone_ages(family) for family in w.FAMILIES}

    plan, total_injections = [], 0
    for family in w.FAMILIES:
        for rv in w.R_V_BRANCHES:
            for subgroup in w.SUBGROUPS:
                nodes = J.truth_age_nodes(
                    age_posterior, subgroup, family, rv, native[family], snap=False
                )
                for age in nodes:
                    cap = min(turnoff_mass(family, age), IMF_UPPER_LIMIT)
                    masses = WP6_MASS_EXTENSION[WP6_MASS_EXTENSION <= cap]
                    total_injections += len(masses) * w.N_INJECT_PER_MASS
                    plan.append(
                        {
                            "subgroup": subgroup, "family": family, "R_V": float(rv),
                            "truth_age_Myr": round(float(age), 4),
                            "turnoff_Msun": round(cap, 1),
                            "extension_masses": [float(m) for m in masses],
                            "n_extension_masses": int(len(masses)),
                        }
                    )
    frame = pd.DataFrame(plan)

    turnoff_by_subgroup = {
        subgroup: {
            "min_turnoff_Msun": round(float(block.turnoff_Msun.min()), 1),
            "median_turnoff_Msun": round(float(block.turnoff_Msun.median()), 1),
            "max_turnoff_Msun": round(float(block.turnoff_Msun.max()), 1),
        }
        for subgroup, block in frame.groupby("subgroup")
    }

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp6_mass_extension_decision.py",
        "status": "DECISION_RECORDED",
        "work_package": "WP6 step 0a",
        "decision": (
            "Extend the injection response above 18 Msun by direct injection "
            "rather than assume the recovery fraction stays flat."
        ),
        "decided_by": "project owner, 2026-07-28",
        "alternative_rejected": (
            "Assert flatness above 18 Msun.  The response does plateau by ~5 "
            "Msun and is flat 8-18 Msun, and the loss is the WP2 quality filter "
            "rather than a magnitude limit, so flatness is plausible.  It was "
            "rejected because the closure ratio is the quantity N_SN rests on, "
            "and replacing a measurement with a plausible assumption at exactly "
            "the load-bearing step is the wrong trade."
        ),
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "design": {
            "frozen_and_not_touched": [
                "wp5_common.MASS_GRID",
                "wp5_common.PARENT_MASS_HI",
                "every repair_v1..v6 artifact",
            ],
            "why_frozen": (
                "the accepted repair_v6 fit filters its response to MASS_GRID's "
                "range in wp5_joint_age_fit._category_counts, so widening the "
                "grid globally would silently change an accepted result the "
                "next time the fit is re-run.  WP5's acceptance must not be a "
                "function of WP6's needs."
            ),
            "extension_grid_Msun": [float(m) for m in WP6_MASS_EXTENSION],
            "imf_upper_limit_Msun": IMF_UPPER_LIMIT,
            "per_node_cap": (
                "each node's extension is truncated at the isochrone turnoff at "
                "that node's truth age.  Above the turnoff the PMS/MS locus does "
                "not exist, so the injector would clamp and place a very massive "
                "star at the turnoff star's magnitude.  The cap is also "
                "physically right: stars above the turnoff are dead, and the "
                "closure test's living-star integral stops there."
            ),
            "outputs_separate": (
                "data/processed/wp6_massext_*_response.parquet, consumed only "
                "by WP6; no WP5 script globs that prefix"
            ),
            "seed_recipe": (
                "fresh default_rng(SEED) per node, matching the repair_v6 node "
                "recipe, so the extension shares donor, binary and extinction "
                "realizations with its parent node"
            ),
        },
        "cost": {
            "nodes": int(len(frame)),
            "total_injected_stars": int(total_injections),
            "extension_masses_per_node_min": int(frame.n_extension_masses.min()),
            "extension_masses_per_node_max": int(frame.n_extension_masses.max()),
            "comparison": (
                f"repair_v6 injected {162 * 45 * w.N_INJECT_PER_MASS:,} stars "
                f"over 45 masses per node; the extension adds "
                f"{total_injections:,}"
            ),
        },
        "turnoff_by_subgroup": turnoff_by_subgroup,
        "defect_found_and_fixed_while_designing_this": {
            "issue": "#14 — the isochrone turnoff was a step function of age",
            "found": (
                "The obvious way to read the turnoff — "
                "load_isochrone_between_ages(family, age)['Mini'].max() — is "
                "wrong.  That loader takes the INTERSECTION of the two "
                "bracketing native tables' mass ranges, so its ceiling is "
                "always the older bracket's turnoff and is constant across each "
                "native interval."
            ),
            "measured": {
                "PARSEC_4.00_Myr_stepped": 48.0,
                "PARSEC_4.20_Myr_stepped": 48.0,
                "MIST_3.00_Myr_stepped": 69.0,
                "MIST_3.50_Myr_stepped": 69.0,
                "PARSEC_4.00_Myr_interpolated": 57.9,
                "MIST_4.00_Myr_interpolated": 58.7,
            },
            "why_it_matters": (
                "N_SN is essentially the IMF integral above the turnoff, so a "
                "stepped turnoff makes N_SN a step function of age — the same "
                "defect class as issue #13, in the quantity the paper is about. "
                "It also manufactured a spurious 22% PARSEC/MIST disagreement "
                "at 4 Myr (48.0 vs 58.6) that vanishes once interpolated "
                "(57.9 vs 58.7)."
            ),
            "fix": (
                "turnoff_mass interpolates log-log between the native tables' "
                "own maxima, using only tabulated values."
            ),
            "monotonicity_correction": (
                "both grids have one unphysical inversion from table sampling "
                "at the top of the locus (PARSEC 331.4 -> 336.0 Msun between "
                "1.995 and 2.239 Myr; MIST 69.0 -> 74.9 between 3.181 and "
                "3.571).  A running minimum is applied because a turnoff cannot "
                "rise with age."
            ),
            "table_ceiling_caveat": (
                "below 3.16 Myr (PARSEC) and 2.52 Myr (MIST) the tabulated "
                "maximum is the table's own mass ceiling (300.0 and 210.2 Msun) "
                "rather than an evolutionary turnoff.  Both are far above the "
                "120 Msun IMF limit, so the value is capped and never used; "
                "CygOB2-C lives entirely in this regime, which is another way "
                "of saying C has lost no stars yet."
            ),
            "wp5_unaffected": (
                "WP5 injects nothing above 18 Msun, far below every turnoff, so "
                "no accepted result depends on this"
            ),
        },
        "predicted_downstream_effects": [
            {
                "id": "D1",
                "target": "WP5 (repair_v6) — the accepted normalization",
                "effect": "NONE, by construction",
                "why": (
                    "MASS_GRID is unchanged and the extension is written to a "
                    "separate file prefix, so no WP5 input changes.  This is "
                    "checkable and will be checked: the repair_v6 gate record "
                    "and mass-function bins must hash identically after the "
                    "extension run."
                ),
            },
            {
                "id": "D2",
                "target": "WP6 closure ratio",
                "effect": "the reason for the work; direction unknown",
                "why": (
                    "the predicted-living integral currently truncates at 18 "
                    "Msun while the observed count reaches 61 Msun, so the two "
                    "sides are not comparable.  Extending makes them comparable. "
                    "Whether the ratio then lands above or below 1 is NOT "
                    "predicted here, and no target value is written down."
                ),
            },
            {
                "id": "D3",
                "target": "N_SN (WP7)",
                "effect": "potentially large, and this is the point",
                "why": (
                    "N_SN is essentially the IMF integral above the turnoff.  "
                    "The closure test is what licenses extrapolating the WP5 "
                    "normalization into that range.  If the measured response "
                    "above 18 Msun differs from the flat assumption, the "
                    "licence changes with it."
                ),
            },
            {
                "id": "D4",
                "target": "CygOB2-C specifically",
                "effect": "largest extension, smallest expected consequence",
                "why": (
                    "C's turnoff sits far above the IMF limit at its 2.5 Myr "
                    "age, so C gets the full 13-point extension while A and B "
                    "are capped near 64-72 Msun.  But C has lost essentially no "
                    "stars yet, so its closure test is the cleanest test of the "
                    "method and the weakest constraint on N_SN."
                ),
            },
            {
                "id": "D5",
                "target": "runtime and reproducibility",
                "effect": "one extra chain step, permanently",
                "why": (
                    "any future WP5 version that WP6 consumes must re-run the "
                    "extension too, because the extension inherits that "
                    "version's extinction and age posteriors.  Recorded here so "
                    "it is not forgotten."
                ),
            },
        ],
        "verification_obligations": [
            "V1 — repair_v6 mass-function bins and gate record must hash "
            "identically after the extension run (D1 is checkable, so check it)",
            "V2 — report the measured recovery fraction above 18 Msun against "
            "the flat-extrapolation assumption that was rejected, so the cost "
            "of the alternative is on the record either way",
            "V3 — no extension mass may exceed its node's turnoff",
        ],
        "node_plan": plan,
    }
    w.write_json(w.PROVENANCE / "wp6_mass_extension_decision.json", record)

    print("WP6 step 0a — mass-extension decision recorded\n")
    print(f"  nodes                 {len(frame)}")
    print(f"  extension masses/node {frame.n_extension_masses.min()}"
          f"-{frame.n_extension_masses.max()}")
    print(f"  injected stars        {total_injections:,}")
    print("\n  isochrone turnoff by subgroup (caps the extension):")
    for subgroup, block in turnoff_by_subgroup.items():
        print(f"    {subgroup}  {block['min_turnoff_Msun']:6.1f} - "
              f"{block['max_turnoff_Msun']:6.1f} Msun  "
              f"(median {block['median_turnoff_Msun']:.1f})")
    print("\nwrote provenance/wp6_mass_extension_decision.json")


if __name__ == "__main__":
    main()
