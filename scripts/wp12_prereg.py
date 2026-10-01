#!/usr/bin/env python3
"""WP12 -- pre-registration for the manuscript-revision analysis.

WP12 is the quantitative work required by
`tasks/manuscript_reframing_revision_brief.md`.  It adds no new observation and
refits nothing upstream: every input is a frozen `repair_v7`-chain product,
consumed read-only.  What it adds is *presentation-grade* structure that the
accepted chain always contained but never exposed --- which branches rest on
mass-function cells that passed their own gate, how the closure test behaves
when the three subgroups are not forced onto one slope, and what the cocoon
scenario score does when its two weakest terms are allowed to vary instead of
being fixed at a bound.

This file is the specification.  It is written and executed BEFORE any WP12
result exists, it records the SHA-256 of every input it will consume, and it
fixes each definition and pass condition in advance.  WP11's preregistration
was criticized for carrying no input hashes (brief item 13); this one does, and
the manuscript describes its status accurately rather than generously.

Prior knowledge is disclosed explicitly in DISCLOSED_PRIOR_KNOWLEDGE below.
Some of these numbers were computed while reading the brief, before this file
was written.  Saying so is the point: a preregistration that hides what its
author already knew is worth less than one that lists it.

Run:
  PYTHONPATH=scripts python3 scripts/wp12_prereg.py
"""
from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone

import wp5_common as w

# --------------------------------------------------------------------- scope
WP5_VERSION = "repair_v7"
WP3_WP4_VERSION = "repair_v5"

# Every input WP12 is allowed to read.  Hashed at preregistration time and
# re-hashed at execution time; a mismatch fails the run.
FROZEN_INPUTS: dict[str, str] = {
    "wp5_normalization": f"data/processed/wp5_imf_normalization_{WP5_VERSION}.parquet",
    "wp5_posterior_draws": f"data/processed/wp5_imf_posterior_draws_{WP5_VERSION}.npz",
    "wp5_fit_execution": f"provenance/wp5_imf_fit_execution_{WP5_VERSION}.json",
    "wp5_gate_record_v6": "provenance/wp5_repair_v6_gate.json",
    "wp6_closure": "tables/wp6_closure_repair_v7.csv",
    "wp6_ledger_execution": "provenance/wp6_ledger_execution.json",
    "wp7_ledger": "tables/wp7_ledger.csv",
    "wp7_branch_sets": "tables/wp7_alpha_headline_branch_sets.csv",
    "wp7_bh_threshold_scan": "tables/wp7_bh_threshold_scan.csv",
    "wp8_crosschecks": "tables/wp8_crosschecks.csv",
    "wp9_verdict": "tables/wp9_verdict.csv",
    "wp9_sensitivity": "tables/wp9_sensitivity.csv",
    "wp11_isotope_forecast": "tables/wp11_isotope_forecast.csv",
    "wp4_age_posteriors": f"data/processed/wp4_age_posteriors_{WP3_WP4_VERSION}.parquet",
    "wp3_extinction": f"data/processed/wp3_extinction_{WP3_WP4_VERSION}.parquet",
    "wp3_isochrones_parsec": "data/processed/wp3_isochrones_parsec.parquet",
    "wp3_isochrones_mist": "data/processed/wp3_isochrones_mist.parquet",
    "wp2_members": "data/processed/wp2_members.parquet",
    "wp2_subgroup_labels": "tables/wp2_subgroup_labels.parquet",
    "wp5_mass_function_bins": f"data/processed/wp5_mass_function_bins_{WP5_VERSION}.parquet",
}

# WP12 must not write to any of these.  Checked by the execution scripts.
FORBIDDEN_TO_WRITE = sorted(FROZEN_INPUTS.values())

# ------------------------------------------------------------------- WP12.1
# The residual-gate landscape.  Nothing is recomputed: the per-cell verdicts are
# read from the repair_v7 normalization, which is the version the ledger
# consumes.  The manuscript has been quoting the repair_v6 *report* table for
# the failure breakdown while the ledger ran on repair_v7; the two agree on
# 40/54 and disagree on which 14 fail.  repair_v7 is authoritative here because
# it is the chain that produced every downstream number.
GATE_LANDSCAPE = {
    "definition": (
        "A subgroup fit cell is (subgroup, family, R_V, alpha).  There are "
        "3 x 2 x 3 x 3 = 54.  A cell passes the WP5 residual gate iff "
        "chi_p >= 0.01 AND trend_p >= 0.05 AND max_abs_pearson_residual <= 3.0 "
        "(scripts/wp5_fit_imf.py).  A family-R_V-alpha COMBINATION passes "
        "all-subgroup iff all three of its cells pass."
    ),
    "quantities": [
        "cells_passing / cells_total, and the failure breakdown by subgroup, "
        "family, R_V and alpha, on repair_v7",
        "combinations passing in all three subgroups, of 18",
        "headline (alpha != 2.6) cells passing, of 36",
        "headline combinations passing all-subgroup, of 12",
        "ledger branches (family x R_V x alpha x sf_duration) built entirely "
        "on passing cells, of 36 headline branches",
    ],
    "why": (
        "The manuscript currently reports '40 of 54 mass-function fits pass' "
        "beside a headline branch set of 36, which invites the reading that "
        "the headline set is validated.  It is not: the headline set is "
        "selected on alpha, not on gate status."
    ),
}

# ------------------------------------------------------------------- WP12.2
# Closure by subgroup and slope, and the mixed-slope sensitivity.
#
# The carried grid forces one alpha on all three subgroups.  That is a modelling
# choice, not a measurement, and the closure test says the subgroups do not want
# the same slope.  The mixed-slope branch relaxes it in the only way the frozen
# products allow: each subgroup keeps its OWN (k, age) posterior draws at its
# OWN closing slope, and the ledger is re-run on that combination.
CLOSURE_AND_MIXED_SLOPE = {
    "closing_alpha_definition": (
        "For each (subgroup, family, R_V) cell, the slope at which the "
        "predicted/observed ratio of living stars above 8 Msun equals 1, by "
        "linear interpolation of log(closure_ratio) against alpha over the "
        "three carried slopes.  Reported as INSIDE or OUTSIDE the carried "
        "[2.0, 2.6] range; an outside value is an extrapolation and is "
        "labelled as one."
    ),
    "mixed_slope_definition": (
        "A mixed-slope branch assigns each subgroup the slope, among the three "
        "carried, whose closure ratio is closest to unity in that subgroup for "
        "the given (family, R_V).  The ledger is then re-run drawing each "
        "subgroup's k and truth-age from the frozen posterior draws at ITS "
        "assigned slope, and sampling its IMF at that same slope.  No new fit "
        "is performed and no draw file is modified."
    ),
    "iterations": 200_000,
    "seed_recipe": "np.random.default_rng(20260803) master, per-subgroup spawn",
    "reproduction_check": (
        "Running the same engine with a single global slope must reproduce the "
        "stored WP7 N_SN_mean for that branch to within Monte-Carlo tolerance. "
        "Declared tolerance: 3% relative on N_SN_mean for every one of the 36 "
        "headline branches, and 0.02 absolute on P(last SN < 100 kyr)."
    ),
}

# ------------------------------------------------------------------- WP12.3
# The cocoon quantity, renamed and made a function of its weakest term.
#
# P_verdict = C1 * C3 * C4 is arithmetically correct and epistemically
# mislabelled.  C4 is an upper bound from a 2-D, footprint-limited runaway
# census applied to already-dead, more massive progenitors; C3 is a
# deterministic mass-to-type map; independence is asserted.  WP12 renames the
# product and reports it against C4 rather than at a fixed C4.
SCENARIO_SCORE = {
    "rename": (
        "conditional scenario-availability score S = C1 x C3 x C4.  NOT a "
        "calibrated probability.  The manuscript must not present it as one."
    ),
    "c4_grid": "np.arange(0.30, 1.0001, 0.005)",
    "c4_thresholds": (
        "C4_any = max over alpha=2.0 headline branches of 0.5 / (C1 * C3): the "
        "C4 below which AT LEAST ONE alpha=2.0 branch falls under 0.5.  "
        "C4_all = min over the same set: the C4 below which EVERY alpha=2.0 "
        "branch falls under 0.5.  The manuscript's current sentence conflates "
        "these two."
    ),
    "also_report": (
        "the C4 at which the first alpha=2.3 branch would RISE above 0.5, "
        "which is > 1 for most branches and is the honest statement that the "
        "alpha split is not an artefact of the adopted C4."
    ),
}

# ------------------------------------------------------------------- WP12.4
# C3 as a model mapping, with a bounded alternative.
#
# The implemented C3 is an indicator: 1 if the progenitor exceeds 30 Msun.  It
# returns exactly 1 because every simulated progenitor does.  That is a
# statement about the mapping, not about nature.  WP12 adds two alternative
# mappings, both bounded and both sourced, and reports C3 as a range.
C3_SUBTYPE = {
    "baseline_map": (
        "deterministic step at 30 Msun (wp9_verdict_prereg.STRIPPED_PROGENITOR"
        "_MSUN).  Retained unchanged as the baseline so the stored WP9 number "
        "is reproduced, not replaced."
    ),
    "alternatives": {
        "threshold_scan": (
            "the same step function with the threshold moved over "
            "[20, 25, 30, 35, 40, 45, 50, 60, 70] Msun.  Bounds how much of C3 "
            "is the 30 Msun choice."
        ),
        "logistic": (
            "a smooth stripping probability p(M) = 1 / (1 + exp(-(M - M0)/dM)) "
            "with M0 in {30, 45, 60} Msun and dM in {5, 10} Msun.  A "
            "single-parameter stand-in for the fact that envelope removal "
            "depends on wind mass loss, metallicity, rotation and binarity "
            "rather than on ZAMS mass alone."
        ),
        "pessimistic_cap": (
            "an explicit worst case in which only progenitors above 60 Msun "
            "are assumed to strip.  This is not a preferred model; it is the "
            "lower edge of the bracket the manuscript must quote."
        ),
    },
    "iterations": 200_000,
    "seed_recipe": "np.random.default_rng(20260804) master, per-subgroup spawn",
    "reproduction_check": (
        "the 30 Msun step must return C3 = 1.000 on every headline branch, "
        "reproducing the stored WP9 value exactly (it is an indicator over a "
        "population whose minimum lies above the threshold, so Monte-Carlo "
        "noise cannot move it)."
    ),
}

# ------------------------------------------------------------------- WP12.5
# The neighbouring-association budget the original plan called WP8.5.
#
# WP8 declined to claim a quantitative budget.  The brief requires either
# completing one or scoping the claim out.  We complete it, coarsely and from a
# single primary tabulation, and we validate the coarse method against our own
# measured ledger before applying it to the neighbours.
NEIGHBOUR_BUDGET = {
    "source": (
        "Martin et al. 2010, A&A 511, A86, Table 1 -- 'Characteristics of the "
        "Cygnus OB associations and stellar clusters collectively referred to "
        "as the Cygnus complex', reproduced there from Knodlseder et al. 2002, "
        "A&A 390, 945.  Each row gives an OBSERVED massive-star count within a "
        "stated ZAMS mass interval, a distance and an age."
    ),
    "method": (
        "For each population: normalize a single power-law IMF of slope alpha "
        "to the observed count over the quoted mass interval, integrate above "
        "the turnoff mass at the quoted age using THIS PROJECT'S OWN turnoff "
        "relation (the same wp6_mass_extension_decision.turnoff_mass used by "
        "the ledger), and call the result the cumulative number of stars that "
        "have already died.  Carry alpha in {2.0, 2.3, 2.6} and both isochrone "
        "families, exactly as the main ledger does."
    ),
    "declared_coarseness": [
        "the counts are OBSERVED, not completeness-corrected, so every derived "
        "supernova number is a lower bound",
        "the ages and distances are heterogeneous literature values from 2002, "
        "not measured here",
        "a single power law is imposed with no measured normalization "
        "uncertainty",
        "coevality is assumed within each population (no star-formation "
        "duration branch)",
        "all deaths are counted as successful explosions, matching the main "
        "ledger's all-explode branch and inheriting its conditionality",
    ],
    "validation": (
        "The identical coarse method is applied to Martin+2010's own Cyg OB2 "
        "row (120 stars in [20, 120] Msun at 2.5 Myr, 1584 pc) and the answer "
        "is compared with this paper's measured ledger.  The coarse estimator "
        "is only usable for the neighbours if it lands within an order of "
        "magnitude of the measured value for the one population we have "
        "measured."
    ),
    "propagated_quantity": (
        "f_OB2 = N_SN(Cyg OB2, measured, in window W) / [ N_SN(Cyg OB2) + "
        "sum over neighbours N_SN(neighbour, coarse, in window W) ], evaluated "
        "as a RATE ratio over the recent window.  This is the fraction of "
        "cavity-region explosions attributable to Cyg OB2 under the coarse "
        "budget.  It is an additional conditional factor on the scenario "
        "score, reported separately and NOT silently multiplied into the "
        "headline number."
    ),
    "membership_of_the_cavity": (
        "Which populations share the Cygnus X cavity is a projection "
        "judgement, not a measurement.  Two nested sets are reported: "
        "STRICT = {Cyg OB1, Cyg OB9} (the two the Cygnus-X literature places "
        "with Cyg OB2 in the same complex) and WIDE = STRICT + {Cyg OB3, "
        "Cyg OB8, Ber 87, NGC 6913, NGC 6910, NGC 6871, Ber 86, IC 4996}.  "
        "Cyg OB7 is excluded from both: at 832 pc it is a foreground object."
    ),
}

# --------------------------------------------------------------- predictions
# Falsifiable statements, fixed here, scored in the execution scripts.  Each
# must be recorded PASS or FAIL and, if it fails, reported as failed.
PREDICTIONS = [
    {
        "id": "R1",
        "statement": (
            "Fewer than half of the 36 headline ledger branches rest entirely "
            "on mass-function cells that pass the WP5 residual gate."
        ),
        "pass_if": "all_pass_headline_branches < 18",
        "falsifies_if": "all_pass_headline_branches >= 18",
        "why_it_matters": (
            "if it fails, the manuscript's current framing is fine and the "
            "brief's item 4 is unnecessary"
        ),
    },
    {
        "id": "R2",
        "statement": (
            "The qualitative alpha split in the cocoon score survives "
            "restriction to all-subgroup-pass branches: at least one alpha=2.0 "
            "branch remains above 0.5 and no alpha=2.3 branch rises above it."
        ),
        "pass_if": (
            "n(alpha=2.0, all-pass, S > 0.5) >= 1 and "
            "n(alpha=2.3, all-pass, S > 0.5) == 0"
        ),
        "falsifies_if": "either arm violates the above",
    },
    {
        "id": "R3",
        "statement": (
            "The three subgroups do not share a closing slope: the spread in "
            "closing alpha across subgroups exceeds 0.15 at the baseline "
            "family and R_V."
        ),
        "pass_if": "max(closing_alpha) - min(closing_alpha) > 0.15 at PARSEC, R_V=3.1",
        "falsifies_if": "the spread is <= 0.15, i.e. one slope does describe all three",
    },
    {
        "id": "R4",
        "statement": (
            "The mixed-slope ledger lies inside the carried headline branch "
            "range for N_SN, i.e. relaxing the single-slope assumption does "
            "not open a new range beyond the one already reported."
        ),
        "pass_if": "NSN(mixed) within [min, max] of the 36 headline branches",
        "falsifies_if": "the mixed-slope value falls outside the carried range",
        "why_it_matters": (
            "if it fails, the paper's quoted range understates its own model "
            "uncertainty and must be widened"
        ),
    },
    {
        "id": "R5",
        "statement": (
            "C3 is not robust to its mapping: under the pessimistic 60 Msun "
            "stripping assumption C3 falls below 0.95 on at least one headline "
            "branch."
        ),
        "pass_if": "min over headline branches of C3(60 Msun step) < 0.95",
        "falsifies_if": "C3 stays above 0.95 even at a 60 Msun threshold",
        "why_it_matters": (
            "if it fails, C3 = 1 is robust to the mapping after all and the "
            "brief's item 9 is a labelling issue only, not a numerical one"
        ),
    },
    {
        "id": "R6",
        "statement": (
            "The coarse neighbour estimator, applied to Martin+2010's own "
            "Cyg OB2 row, lands within a factor of 10 of this paper's measured "
            "baseline N_SN."
        ),
        "pass_if": "1/10 <= N_coarse(CygOB2) / N_measured(baseline) <= 10",
        "falsifies_if": "the coarse estimator misses by more than an order of magnitude",
        "why_it_matters": (
            "this is the gate on WP12.5.  If it FAILS, the neighbour budget is "
            "not reported as a number and the paper scopes the claim out "
            "instead, per the brief's second option."
        ),
    },
    {
        "id": "R7",
        "statement": (
            "Cyg OB2 does not dominate the recent-explosion budget of the "
            "Cygnus complex: f_OB2 < 0.9 on the WIDE cavity set at the "
            "baseline slope."
        ),
        "pass_if": "f_OB2(WIDE, alpha=2.3) < 0.9",
        "falsifies_if": "f_OB2 >= 0.9, i.e. the neighbours are negligible",
    },
]

# ------------------------------------------------------- disclosed priors
DISCLOSED_PRIOR_KNOWLEDGE = [
    "The failure breakdown of the repair_v7 gate (14 failures: B 7, A 4, C 3; "
    "MIST 10, PARSEC 4) was computed while reading the revision brief, before "
    "this preregistration was written.  R1 is therefore NOT blind: the value "
    "15/36 was known.  R1 is retained as a recorded expectation, not as a "
    "test, and is labelled as such in the outcome.",
    "The C4 thresholds 0.7198 and 0.5806 were likewise computed before this "
    "file existed.  They are arithmetic consequences of the stored WP9 table, "
    "not a new measurement, and no prediction is registered on them.",
    "The revision brief itself asserts several of these numbers.  Where WP12's "
    "recomputation disagrees with the brief, WP12's value from the frozen "
    "artifact wins and the disagreement is recorded.",
    "R3-R7 were written before the corresponding quantities were computed.",
]

ANTI_TUNING_RULE = (
    "No WP12 script may alter any upstream artifact, threshold, gate "
    "definition, branch definition or headline branch set.  If a WP12 result "
    "is unwelcome it is reported, not repaired.  Failed predictions stay "
    "failed."
)


def build() -> dict:
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp12_prereg.py",
        "work_package": "WP12 -- manuscript revision analysis",
        "brief": "tasks/manuscript_reframing_revision_brief.md",
        "status": "PREREGISTRATION -- written before any WP12 result existed",
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "scope": (
            "Read-only over the frozen repair_v7 chain.  WP12 adds no "
            "observation, refits nothing upstream, and changes no accepted "
            "number.  It exposes structure the accepted chain already had."
        ),
        "versions": {"wp5": WP5_VERSION, "wp3_wp4": WP3_WP4_VERSION},
        "frozen_inputs": {
            name: {
                "path": rel,
                "sha256": (
                    w.sha256(w.ROOT / rel) if (w.ROOT / rel).exists() else None
                ),
                "exists": (w.ROOT / rel).exists(),
            }
            for name, rel in sorted(FROZEN_INPUTS.items())
        },
        "forbidden_to_write": FORBIDDEN_TO_WRITE,
        "wp12_1_gate_landscape": GATE_LANDSCAPE,
        "wp12_2_closure_and_mixed_slope": CLOSURE_AND_MIXED_SLOPE,
        "wp12_3_scenario_score": SCENARIO_SCORE,
        "wp12_4_c3_subtype": C3_SUBTYPE,
        "wp12_5_neighbour_budget": NEIGHBOUR_BUDGET,
        "predictions": PREDICTIONS,
        "disclosed_prior_knowledge": DISCLOSED_PRIOR_KNOWLEDGE,
        "anti_tuning_rule": ANTI_TUNING_RULE,
    }


# Products repair_v7 wrote without a suffix; a later chain reads its tagged copy.
CHAIN_UNVERSIONED = {
    "wp6_ledger_execution", "wp7_ledger", "wp7_branch_sets",
    "wp7_bh_threshold_scan", "wp8_crosschecks", "wp9_verdict",
    "wp9_sensitivity", "wp11_isotope_forecast",
}


def chain_record() -> None:
    """Hash record for a later chain (issue #19), from the SAME specification.

    The 2026-08-03 preregistration is not edited.  This copies it verbatim --
    definitions, predictions R1-R7, disclosed priors, anti-tuning rule -- and
    replaces only ``frozen_inputs``: each path re-pointed to the chain's own
    product and re-hashed now.  WP12 on that chain then verifies against this
    record exactly as WP12 on repair_v7 verifies against the original.
    """
    import copy
    import json

    import chain as C

    if C.CHAIN == C.LEGACY:
        raise SystemExit("--chain-record is for a chain other than repair_v7")
    original_path = w.PROVENANCE / "wp12_revision_prereg.json"
    out = w.PROVENANCE / f"wp12_revision_prereg_{C.CHAIN}.json"
    if out.exists():
        raise SystemExit(f"{out.relative_to(w.ROOT)} exists; a hash record is written once")
    original = json.loads(original_path.read_text())
    record = copy.deepcopy(original)
    inputs = {}
    for name, entry in sorted(original["frozen_inputs"].items()):
        rel = entry["path"].replace(f"_{WP5_VERSION}", f"_{C.V['wp5']}")
        if rel == entry["path"] and name in CHAIN_UNVERSIONED:
            rel = C.tag_rel(rel)
        path = w.ROOT / rel
        inputs[name] = {
            "path": rel,
            "sha256": w.sha256(path) if path.exists() else None,
            "exists": path.exists(),
            "repair_v7_path": entry["path"],
        }
    record.update({
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp12_prereg.py --chain-record",
        "status": (
            f"HASH RECORD for {C.CHAIN} -- the specification and predictions are "
            "the 2026-08-03 preregistration's, copied verbatim; only "
            "frozen_inputs is re-pointed and re-hashed"
        ),
        "chain": C.CHAIN,
        "derived_from": {
            "path": str(original_path.relative_to(w.ROOT)),
            "sha256": w.sha256(original_path),
            "created_utc": original["created_utc"],
        },
        "reason": (
            "issue #19: the repair_v7 chain consumed anchor masses read at "
            "pre-repair ages; see provenance/issue19_repair_v8_prereg.json"
        ),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "versions": {**original["versions"], "chain": C.CHAIN, "wp5": C.V["wp5"],
                     "wp4_masses": C.V["wp4_masses"]},
        "frozen_inputs": inputs,
        "forbidden_to_write": sorted(e["path"] for e in inputs.values()),
        "predictions_note": (
            "R1-R7 were written before any repair_v7 WP12 result existed.  On "
            f"{C.CHAIN} the repair_v7 WP12 results ARE known; re-scoring R1-R7 "
            "here is a re-execution of a fixed specification, not a fresh test."
        ),
    })
    voided = sorted(w.PROVENANCE.glob(f"wp12_revision_prereg_{C.CHAIN}_void_*.json"))
    if voided:
        record["voided_predecessors"] = {
            str(p.relative_to(w.ROOT)): w.sha256(p) for p in voided
        }
        record["voided_note"] = (
            "issue #19, 2026-10-01: the first repair_v8 hash record pinned a "
            "WP11 forecast run at 200,000 iterations instead of the published "
            "500,000.  WP11 was re-run at 500,000; the first record was renamed "
            "(not deleted) and this one written in its place."
        )
    missing = [n for n, e in inputs.items() if not e["exists"]]
    if missing:
        raise SystemExit(f"chain inputs missing: {missing}")
    w.write_json(out, record)
    print(f"wrote {out.relative_to(w.ROOT)} ({len(inputs)} inputs hashed)")


def main() -> None:
    if "--chain-record" in sys.argv[1:]:
        chain_record()
        return
    record = build()
    out = w.PROVENANCE / "wp12_revision_prereg.json"
    if out.exists():
        raise SystemExit(
            f"{out} already exists.  A preregistration is written once.  "
            "Delete it deliberately if you really mean to re-specify."
        )
    w.write_json(out, record)
    missing = [
        name for name, entry in record["frozen_inputs"].items()
        if not entry["exists"]
    ]
    print("WP12 preregistration written to provenance/wp12_revision_prereg.json")
    print(f"  frozen inputs hashed : {len(record['frozen_inputs'])}")
    if missing:
        print(f"  MISSING              : {', '.join(missing)}")
    print(f"  predictions declared : {len(PREDICTIONS)}")
    print(f"  priors disclosed     : {len(DISCLOSED_PRIOR_KNOWLEDGE)}")
    if missing:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
