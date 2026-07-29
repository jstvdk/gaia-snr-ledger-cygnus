#!/usr/bin/env python3
"""Pre-registration: WP7, the supernova ledger.

Written BEFORE the Monte Carlo engine is run.  Predictions, thresholds and
decision rules are fixed here and are not amended afterwards.

WHAT WP7 IS
-----------
N_SN posterior, explosion timeline R_SN(t), and time-since-last-supernova
posterior, per branch and as explicit functions of assumed age.  It is pure
computation on frozen inputs: nothing is refitted, no new data is read.

TWO ROUTES TO N_SN, AND WHY ONLY ONE OF THEM WORKS
--------------------------------------------------
The execution plan specifies a "missing = dead" bookkeeping with the recovered
runaways subtracted.  That route is NOT well-posed with the normalization this
project actually has, and the arithmetic was done at design time rather than
discovered after a run:

    turnoff route,   sum_sg k_sg * integral[turnoff, 120] M^-alpha dM =    8.46
    census route vs the LABELLED population (246.4)               =   17.4
    census route vs the FULL ledger        (380.6)                = -116.8

A negative death count is impossible, and the reason is structural.  Each k_sg
is fitted from that subgroup's own 2-8 Msun members, so sum_sg k_sg * integral
predicts the LABELLED, clustered population only.  The ledger's 380.6 adds 49
unlabelled members, 27 orphan anchors and 54.9 runaways that no per-subgroup
normalization ever predicted.  WP5 never fitted an association-wide k, so the
subtraction has no consistent left-hand side.

Consequences, fixed here:

  * the TURNOFF ROUTE is the measurement;
  * the census route is reported as a CONSISTENCY CHECK on the labelled
    population only, where it reduces to the WP6 closure test plus a small dead
    sliver, and where its residual (17.4 vs 8.46) is dominated by the 6.7%
    closure excess rather than by anything about supernovae;
  * the runaway count does NOT enter N_SN.  Runaways are living stars below the
    turnoff; they change neither k (fitted below 8 Msun) nor the turnoff.  What
    they bound is the fraction of supernovae that occurred OUTSIDE the
    association, which is a WP8 question about location, not a WP7 question
    about number.

This supersedes execution-plan WP7 step 2's runaway subtraction.  The plan's
instruction is not wrong in principle -- it is right for an association-wide
normalization that does not exist.

Output: provenance/wp7_ledger_prereg.json

Run:
  WP_REPAIR_VERSION=repair_v5 WP3_ANCHOR_PRIOR_MODE=kriging \
  PYTHONPATH=scripts python3 scripts/wp7_ledger_prereg.py
"""
from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w
from wp6_mass_extension_decision import IMF_UPPER_LIMIT

# ---- branch grid (all Class E, carried in parallel, never averaged) --------
SF_DURATIONS_MYR = (0.0, 1.0, 2.0)          # CUTS master-table row 27
EXPLODABILITY = ("all_explode", "islands")
BH_THRESHOLD_MSUN = 25.0                    # the "islands" branch cut
BH_THRESHOLD_SCAN = (20.0, 25.0, 30.0, 40.0, 60.0, 80.0, 120.0)

SN_THRESHOLD_MSUN = 8.0
N_ITERATIONS = 40000
CONVERGENCE_SPLIT = 20000
AGE_SCAN_MYR = tuple(np.round(np.arange(2.0, 6.01, 0.25), 2))
RECENT_WINDOW_MYR = 0.1                     # the "last SN < 100 kyr" question

# Ledger channel weights, from provenance/wp6_ledger_execution.json.
LEDGER_LABELLED = 246.4
LEDGER_UNLABELLED_MEMBERS = 49.0
LEDGER_ORPHAN_ANCHORS = 27.0
LEDGER_RUNAWAY_UNCLIPPED = 54.9
LEDGER_RUNAWAY_BINNED = 58.18


def main() -> None:
    living_excl_runaway = (
        LEDGER_LABELLED + LEDGER_UNLABELLED_MEMBERS + LEDGER_ORPHAN_ANCHORS
    )
    runaway_fraction = LEDGER_RUNAWAY_UNCLIPPED / (
        living_excl_runaway + LEDGER_RUNAWAY_UNCLIPPED
    )

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "script": "scripts/wp7_ledger_prereg.py",
        "status": "PREREGISTERED",
        "work_package": "WP7",
        "objective": (
            "N_SN posterior, R_SN(t) explosion timeline, and time-since-last-SN "
            "posterior, per branch and as explicit functions of assumed age"
        ),
        "environment": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
        "inputs_are_frozen": {
            "normalization": "data/processed/wp5_imf_normalization_repair_v7.parquet",
            "paired_posterior_draws": (
                "data/processed/wp5_imf_posterior_draws_repair_v7.npz -- 10000 "
                "PAIRED (k, truth_age) draws per (subgroup, family, R_V, alpha).  "
                "Pairing matters: k and the age are correlated through the WP5 "
                "fit, and drawing them independently would break that "
                "correlation and misstate the N_SN interval."
            ),
            "turnoff": (
                "wp6_mass_extension_decision.turnoff_mass, the SAME relation WP4 "
                "and WP6 use.  Lifetimes tau(m) are obtained by inverting it "
                "rather than by importing an external lifetime table, so the "
                "ledger cannot disagree with the closure test for the trivial "
                "reason of a different tau(m)."
            ),
            "nothing_is_refitted": (
                "WP7 reads frozen artifacts only.  In particular alpha is NOT "
                "refitted -- the closure test is the analysis's only "
                "out-of-sample check and converting it into a fitted slope "
                "would destroy it (PROJECT_TRACE section 10c item 2)."
            ),
        },
        "design_finding_measured_before_preregistration": {
            "what": (
                "the execution plan's 'missing = dead' route is not well-posed "
                "with the normalization this project has"
            ),
            "arithmetic": {
                "turnoff_route_N_SN": 8.46,
                "census_route_vs_labelled_246p4": 17.4,
                "census_route_vs_full_ledger_380p6": -116.8,
                "branch": "PARSEC, R_V=3.1, alpha=2.3, k_median, posterior-mean ages",
            },
            "why": (
                "each k_sg is fitted from that subgroup's own 2-8 Msun members, "
                "so sum_sg k_sg * integral predicts the LABELLED clustered "
                "population (246.4).  The ledger's 380.6 adds 49 unlabelled "
                "members, 27 orphan anchors and 54.9 runaways that no "
                "per-subgroup normalization predicted.  WP5 never fitted an "
                "association-wide k."
            ),
            "honesty_note": (
                "this was computed BEFORE this pre-registration was written and "
                "is recorded as a measured design finding, NOT dressed up as a "
                "prediction.  The predictions below are restricted to things the "
                "Monte Carlo has not yet been asked."
            ),
            "reconciliation_of_the_two_living_totals": (
                "246.4 (WP6 closure census, per subgroup) + 49.0 (members with "
                "no subgroup label) = 295.4 (ledger member channel).  The "
                "closure test filters on subgroup membership; the ledger filters "
                "on having a finite P(M>8).  Both are correct for their own "
                "purpose and they are not interchangeable."
            ),
        },
        "decisions_fixed_here": {
            "measurement_route": "turnoff",
            "census_route_status": (
                "reported as a consistency check on the LABELLED population "
                "only, with its residual attributed to the WP6 closure excess "
                "rather than to supernovae"
            ),
            "runaways_do_not_enter_N_SN": (
                "runaways are living stars below the turnoff.  They change "
                "neither k (fitted from 2-8 Msun) nor the turnoff, so they "
                "cannot change a number that is an integral of the IMF above "
                "the turnoff.  This SUPERSEDES execution-plan WP7 step 2, whose "
                "subtraction is correct only for an association-wide "
                "normalization that does not exist."
            ),
            "what_runaways_DO_bound": (
                f"the fraction of supernovae that occurred OUTSIDE the "
                f"association.  Stars ejected before their death exploded "
                f"elsewhere.  With {LEDGER_RUNAWAY_UNCLIPPED} runaways against "
                f"{living_excl_runaway:.1f} retained living stars above 8 Msun, "
                f"at most {runaway_fraction:.1%} of the ledger's supernovae "
                f"were not in-situ.  Reported as a bound on location, carried "
                f"into WP8's cavity argument; NOT subtracted from N_SN."
            ),
            "runaway_systematic_still_carried": (
                f"the clipped/unclipped runaway totals "
                f"({LEDGER_RUNAWAY_BINNED} vs {LEDGER_RUNAWAY_UNCLIPPED}) "
                f"differ by 3.3, which moves the out-of-association bound by "
                f"under a percentage point.  Carried, but immaterial."
            ),
            "birth_time_convention": (
                "the WP4 age is treated as the MIDPOINT of the formation "
                "interval, so stars are born uniformly in "
                "[t_age - delta/2, t_age + delta/2].  This keeps the mean age "
                "equal to the fitted age and makes the star-formation duration "
                "a pure spread rather than a systematic shift.  The alternative "
                "convention (t_age = formation start) would make every "
                "non-zero duration reduce N_SN by construction, confounding the "
                "branch with a bias."
            ),
            "death_test_avoids_inverting_tau": (
                "a star of mass m born at look-back t_b is dead iff "
                "m > turnoff(t_b).  This needs no lifetime inversion and is "
                "exactly consistent with WP4/WP6 by construction.  tau(m) is "
                "inverted ONLY to date the explosion of stars already known to "
                "have died, where the turnoff relation is single-valued and "
                "well inside its validity range."
            ),
            "sampling_is_discrete": (
                "N_born above the minimum turnoff is drawn as Poisson, and "
                "individual masses are drawn from the IMF.  At N ~ 8 the "
                "posterior is asymmetric and expectation values alone would "
                "misrepresent it.  Only stars above the minimum turnoff over "
                "the birth window are drawn, which is exact -- no star below it "
                "can have died -- and keeps the engine cheap."
            ),
        },
        "branch_grid": {
            "families": list(w.FAMILIES),
            "R_V": list(w.R_V_BRANCHES),
            "alpha": list(w.IMF_SLOPES),
            "sf_duration_Myr": list(SF_DURATIONS_MYR),
            "explodability": list(EXPLODABILITY),
            "n_branches": (
                len(w.FAMILIES) * len(w.R_V_BRANCHES) * len(w.IMF_SLOPES)
                * len(SF_DURATIONS_MYR) * len(EXPLODABILITY)
            ),
            "subgroups": list(w.SUBGROUPS),
            "alpha_2p6_is_disfavoured": (
                "carried, but the WP6 census disfavours it.  Reported per "
                "branch and never averaged across alpha."
            ),
        },
        "explodability": {
            "all_explode": "every star that has died and exceeds 8 Msun explodes",
            "islands": (
                f"no supernova above {BH_THRESHOLD_MSUN} Msun (direct collapse "
                f"to a black hole); explodes below it"
            ),
            "why_the_fine_island_structure_is_omitted": (
                "Sukhbold et al. 2016 and Ertl et al. 2016 give an interleaved "
                "pattern of explodable and non-explodable ZAMS masses between "
                "roughly 15 and 25 Msun.  That structure is IRRELEVANT here: "
                "the smallest turnoff anywhere on the branch grid is about "
                "52 Msun, so every star that has died in Cyg OB2 is far above "
                "the island region.  Rendering the islands in detail would be "
                "fake precision.  The branch reduces to a single question -- do "
                "stars well above 40 Msun explode -- and that is what is varied."
            ),
            "sensitivity_is_scanned_not_hidden": (
                f"N_SN is reported against a black-hole threshold scanned over "
                f"{list(BH_THRESHOLD_SCAN)} Msun, so the reader sees the "
                f"dependence directly instead of a binary branch label"
            ),
        },
        "monte_carlo": {
            "iterations": N_ITERATIONS,
            "convergence_split": CONVERGENCE_SPLIT,
            "seed": int(w.SEED),
            "posterior_draws_available": 10000,
            "resampling": (
                "iterations beyond the 10000 available posterior draws resample "
                "the (k, age) pairs with replacement.  This adds no posterior "
                "information -- it only reduces Monte Carlo noise from the "
                "Poisson and IMF sampling, which is the intended effect, and it "
                "is stated so no one reads 40000 as 40000 independent posterior "
                "samples."
            ),
        },
        "predictions": [
            {
                "id": "L1",
                "statement": (
                    "CygOB2-C contributes exactly zero supernovae on every "
                    "PARSEC branch and a non-zero number on every MIST branch.  "
                    "The two families disagree about C's age (2.51 vs 3.16 Myr) "
                    "and that disagreement straddles the age at which the "
                    "PARSEC turnoff falls below the 120 Msun IMF limit."
                ),
                "falsifies_if": (
                    "any PARSEC-C iteration yields N_SN > 0, or the MIST-C "
                    "posterior median is zero"
                ),
                "if_falsified": (
                    "the family branch is behaving differently from what the "
                    "point-estimate turnoffs imply and the age-to-turnoff "
                    "mapping must be re-examined before N_SN is reported"
                ),
                "why_it_matters": (
                    "for CygOB2-C the isochrone family is not a minor "
                    "systematic -- it is the difference between 'C has produced "
                    "no supernovae' and 'C has produced some'.  This must be "
                    "reported as a branch disagreement, never averaged away."
                ),
            },
            {
                "id": "L2",
                "statement": (
                    f"the islands branch yields N_SN = 0 in EVERY iteration of "
                    f"every branch, because the smallest turnoff on the grid is "
                    f"about 52 Msun and nothing below the "
                    f"{BH_THRESHOLD_MSUN} Msun black-hole threshold has had time "
                    f"to die"
                ),
                "falsifies_if": "any islands iteration yields N_SN > 0",
                "if_falsified": (
                    "some branch is reaching ages far older than the WP4 "
                    "posteriors support, which would indicate an error in the "
                    "star-formation-duration handling"
                ),
                "why_it_matters": (
                    "it means the entire supernova budget of Cyg OB2 is "
                    "conditional on whether stars above 40 Msun explode at all.  "
                    "That is a genuine, publishable conditional -- and it sets "
                    "up WP8's strongest test, because PSR J2032+4127 is a "
                    "neutron star, and a neutron star requires a SUCCESSFUL "
                    "explosion.  Its existence is evidence against the pure "
                    "islands branch, or evidence that its progenitor was less "
                    "massive than the present turnoff implies."
                ),
            },
            {
                "id": "L3",
                "statement": (
                    "N_SN INCREASES with star-formation duration.  N_SN is a "
                    "convex function of age near the turnoff, so spreading "
                    "births symmetrically about a fixed mean age raises the "
                    "mean count by Jensen's inequality."
                ),
                "falsifies_if": (
                    "the baseline-branch N_SN median at delta = 2 Myr is below "
                    "the value at delta = 0"
                ),
                "if_falsified": (
                    "either the convexity argument is wrong in this regime or "
                    "the midpoint birth-time convention was not implemented as "
                    "specified; both require diagnosis before reporting"
                ),
            },
            {
                "id": "L4",
                "statement": (
                    "the N_SN posterior is right-skewed: mean > median at the "
                    "baseline branch, from Poisson sampling plus turnoff "
                    "convexity"
                ),
                "falsifies_if": "median >= mean at the baseline branch",
                "if_falsified": (
                    "the small-number stochasticity is behaving unexpectedly "
                    "and the interval reporting must be re-derived"
                ),
            },
            {
                "id": "L5",
                "statement": (
                    f"P(last SN < {RECENT_WINDOW_MYR * 1000:.0f} kyr) lies "
                    f"between 0.30 and 0.70 for the association total on the "
                    f"all-explode baseline branch.  Order-of-magnitude: about "
                    f"8.5 supernovae spread over the roughly 1 Myr since the "
                    f"turnoff first fell below 120 Msun gives a rate near "
                    f"8/Myr, so 1 - exp(-0.85) ~ 0.57."
                ),
                "falsifies_if": "the value falls outside [0.30, 0.70]",
                "if_falsified": (
                    "the explosion-epoch distribution is not approximately "
                    "uniform over the recent past as the rate argument assumes, "
                    "which is itself worth reporting"
                ),
            },
            {
                "id": "L6",
                "statement": (
                    "the Monte Carlo is converged: every reported N_SN quantile "
                    "moves by less than 1% between the first "
                    f"{CONVERGENCE_SPLIT} iterations and all {N_ITERATIONS}"
                ),
                "falsifies_if": "any reported quantile moves by 1% or more",
                "if_falsified": "iterations are increased until it passes",
            },
        ],
        "gate_criteria": {
            "G7a": (
                "Monte Carlo converged -- posteriors stable under doubling "
                "iterations (prediction L6)"
            ),
            "G7b": (
                "the branch spread is documented and no branch is averaged "
                "away, with the family disagreement on CygOB2-C reported "
                "explicitly"
            ),
            "G7c": (
                "the Knoedlseder+2002 / Martin+2010 comparison is written, with "
                "any factor-level difference explained by input differences "
                "(census size, ages, mass range) rather than hand-waved"
            ),
            "G7d": (
                "the census route is reported with its failure diagnosed, not "
                "silently dropped"
            ),
        },
        "external_anchors_for_G7c": {
            "Knoedlseder_2002": (
                "models Cyg OB2 as 1-4 Myr old with the first supernovae near "
                "4 Myr in the coeval picture.  Our A and B sit at 4.00 and "
                "4.09 Myr and produce a small non-zero count; C at 2.52 Myr "
                "(PARSEC) produces none.  Qualitative agreement is expected and "
                "its absence would be a tension."
            ),
            "Martin_2010": (
                "uses 120 stars in 20-120 Msun at 1584 pc and 2.5 Myr for the "
                "Cyg OB2 row, and reports 10-20 supernovae over the last Myr "
                "for the WHOLE Cygnus complex.  Not directly comparable: that "
                "is a complex-wide figure over a different volume.  The "
                "comparison must be made on rate per association, and the "
                "distance and age differences stated."
            ),
        },
        "what_wp7_cannot_do": [
            "it cannot measure N_SN independently of the WP4 ages -- everything "
            "inherits them, which is why the age-sensitivity scan is mandatory "
            "rather than optional",
            "it cannot locate the explosions; the turnoff route counts "
            "supernovae from the association's stellar population wherever they "
            "occurred, and the in-situ fraction is bounded, not computed",
            "it does not model binary mass transfer, which can both strip "
            "envelopes and rejuvenate accretors; acknowledged as unmodelled",
            "it cannot distinguish a failed supernova from a successful one -- "
            "that is the explodability branch, and the branch spread is the "
            "honest answer",
        ],
        "age_sensitivity": {
            "scan_Myr": list(AGE_SCAN_MYR),
            "mandatory": (
                "N_SN and P(last SN < 100 kyr) are reported as explicit "
                "FUNCTIONS of assumed age in addition to the marginalized "
                "values.  This is the honesty plot for the age sensitivity and "
                "is not optional."
            ),
        },
    }
    w.write_json(w.PROVENANCE / "wp7_ledger_prereg.json", record)

    print("WP7 — the supernova ledger — pre-registered\n")
    print("  measurement route : turnoff")
    print("  census route      : consistency check only (not well-posed at")
    print("                      association level -- see the design finding)")
    print(f"  runaways          : bound the OUT-OF-ASSOCIATION fraction at "
          f"{runaway_fraction:.1%}, not subtracted from N_SN")
    print(f"\n  branches          : {record['branch_grid']['n_branches']} "
          f"x {len(w.SUBGROUPS)} subgroups")
    print(f"  iterations        : {N_ITERATIONS}")
    print("\n  predictions:")
    for entry in record["predictions"]:
        head = entry["statement"].split(".")[0]
        print(f"    {entry['id']}  {head[:66]}")
    print("\nwrote provenance/wp7_ledger_prereg.json")


if __name__ == "__main__":
    main()
