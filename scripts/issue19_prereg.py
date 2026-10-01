#!/usr/bin/env python3
"""Issue #19 -- pre-registration of the repair_v8 rerun.

repair_v8 regenerates the products that consumed the two pre-repair WP4
artifacts found in issue #19 (tasks/issue19_stale_wp4_inputs_brief.md):

  F2  the 150 spectroscopic-HRD anchor masses, read off isochrones at the
      PRE-REPAIR ages and copied unchanged onto every R_V branch;
  F1  the 2.25-5.67 Myr "two-indicator" age envelope (narrative only).

This file is written and executed BEFORE any repair_v8 product exists.  It
records the SHA-256 of every input the rerun will consume and of every
repair_v7 comparison artifact the predictions are scored against, fixes the
predictions P1-P6 and the exact rule that scores each, and fixes the adoption
rule.  It refuses to overwrite itself: a preregistration is written once.

Decisions D1-D3 of the brief (section 4) were taken at the brief's
recommendation when the user instructed the agent to proceed with it; they are
recorded here so that they are on the record before the rerun, not after.

Run:
  PYTHONPATH=scripts python3 scripts/issue19_prereg.py
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import wp5_common as w

OUT = w.PROVENANCE / "issue19_repair_v8_prereg.json"

NEW_VERSION = "repair_v8"
WP3_WP4_VERSION = "repair_v5"   # WP3 extinction and WP4 ages: unchanged
RESPONSE_VERSION = "repair_v7"  # injection responses: reused, see REUSE below
COMPARE_VERSION = "repair_v7"   # the chain every prediction is scored against

# ----------------------------------------------------------------- decisions
DECISIONS = {
    "D1": {
        "question": "fix F3 now or separately?",
        "taken": (
            "F3(b) now: anchor masses become R_V-dependent (one mass per "
            "family x R_V branch, read at that branch's own repair_v5 age with "
            "that branch's own de-reddened M_G0).  F3(a) -- initial mass Mini "
            "instead of present-day Mass -- is NOT part of repair_v8; it is "
            "deferred to a separately pre-registered, separately scored change "
            "so the two effects are attributable."
        ),
        "source": "brief section 4, recommended option",
    },
    "D2": {
        "question": "order relative to WP13",
        "taken": (
            "WP13 M0/M1 are read only on the repaired chain.  Freezing the "
            "WP13 materiality thresholds is a WP13 task and is not done here."
        ),
        "source": "brief section 4, recommended option",
    },
    "D3": {
        "question": "version name",
        "taken": (
            "repair_v8 for every regenerated product.  Every repair_v7-and-"
            "earlier artifact is left byte-identical (rule 3); unversioned "
            "repair_v7 tables are not rewritten -- repair_v8 writes "
            "_repair_v8-suffixed siblings."
        ),
        "source": "brief section 4, recommended option",
    },
    "taken_by": (
        "the agent, on the user's instruction to proceed with the brief "
        "(2026-10-01); flagged to the user in the completion report"
    ),
}

# ---------------------------------------------------------------- what moves
CHANGE = {
    "wp4_anchors": (
        "scripts/wp4_anchors_hrd.py gains --age-version (no default; the "
        "unversioned pre-repair posterior is refused) and --extinction-version. "
        "For each anchor and each (family, R_V) branch: age = the repair_v5 "
        "upper-MS MAP of that subgroup/family/R_V at f_bin = 0.4, dmu = 0, "
        "snapped to the nearest native isochrone age (the frozen procedure's "
        "rule, unchanged); M_G0 = G0_abs_rv{R_V} from wp3_extinction_repair_v5; "
        "mass = present-day Mass of the nearest-chi isochrone point (unchanged "
        "procedure, SIG_LOGTE = 0.03, SIG_MG0 = 0.40).  Unlabelled anchors use "
        "the median of the three subgroup MAPs for that family and R_V, as "
        "before.  Output data/processed/wp4_anchor_hrd_repair_v8.parquet."
    ),
    "wp4_masses": (
        "scripts/wp4_mass_posteriors_repair.py reads the versioned anchor "
        "file and uses each branch's own anchor mass.  Photometric posteriors "
        "are recomputed by the identical code, seeds and inputs as repair_v5. "
        "Stale label photometric_posterior_repair_v1 -> photometric_posterior "
        "(cosmetic).  Output wp4_mass_posteriors_repair_v8.parquet + samples."
    ),
    "downstream": (
        "WP5 joint age-k fit, WP6 census / closure / attribution / living "
        "ledger, WP7 ledger + scans + headline branch sets + binary bound, "
        "WP8 cross-checks, WP9 verdict, WP4-WP5 age and association-mass "
        "reconciliations, WP12 analyses: re-run unchanged in method on the "
        "repair_v8 masses.  Nothing is retuned; no threshold moves."
    ),
}

REUSE = {
    "injection_responses": (
        "The repair_v7 WP5 node responses (wp5_agenode_*_repair_v7_*) and WP6 "
        "mass-extension responses (wp6_massext_*_repair_v7_*) are reused "
        "unchanged.  Reason, verified by code inspection on 2026-10-01: "
        "wp5_injections_agenodes.py / wp5_injections_repair.py / "
        "wp6_massive_injections.py never open wp4_mass_posteriors_* -- it "
        "appears only in their provenance input list -- and the truth-age "
        "nodes come from wp4_age_posteriors_repair_v5, which repair_v8 does "
        "not change.  Recovered masses for injected stars come from "
        "wp4_repair_common.infer_mass_samples, the photometric path, which "
        "never uses the anchor table."
    ),
    "runaways": (
        "wp6_runaways / wp6_runaway_crossmatch read Gaia astrometry and "
        "membership only; not re-run.  Their tables feed the repair_v8 living "
        "ledger unchanged."
    ),
}

# ------------------------------------------------------------- predictions
PREDICTIONS = [
    {
        "id": "P1",
        "statement": (
            "PARSEC, subgroup A: the WP5 calibration window loses 8.0 +/- 0.5 "
            "membership-weighted stars on every PARSEC branch."
        ),
        "measured_quantity": (
            "membership_weighted_calibration_sources(repair_v7) - "
            "membership_weighted_calibration_sources(repair_v8) for "
            "subgroup CygOB2-A, family PARSEC, at each R_V in {3.0, 3.1, 3.5} "
            "(alpha = 2.3 row; the window count does not depend on alpha)"
        ),
        "pass_if": "7.5 <= loss <= 8.5 at all three PARSEC R_V branches",
    },
    {
        "id": "P2",
        "statement": (
            "MIST branches, and subgroups B and C on every branch: k and "
            "N_death within 1 %."
        ),
        "measured_quantity": (
            "(a) |k_median(v8)/k_median(v7) - 1| for every WP5 cell with "
            "subgroup in {B, C} (36 cells) or family MIST (27 cells; 45 "
            "distinct cells in total); (b) N_SN_mean, all_explode, every "
            "sf_duration, for subgroup-scope rows of B and C on every branch "
            "and for every MIST row (subgroup and association scope)"
        ),
        "pass_if": (
            "every cell in (a) within 1 % relative, AND every row in (b) "
            "within 1 % relative or within 0.01 deaths absolute where the "
            "repair_v7 value is below 1 (C is 0.00 at baseline, where a "
            "relative tolerance is undefined)"
        ),
    },
    {
        "id": "P3",
        "statement": "Baseline association N_death lies in [8.20, 8.43].",
        "measured_quantity": (
            "N_SN_mean of the association-scope, all_explode ledger row at "
            "PARSEC / R_V 3.1 / alpha 2.3 / sf_duration 0, repair_v8"
        ),
        "pass_if": "8.20 <= value <= 8.43 (the repair_v7 value is 8.4349)",
    },
    {
        "id": "P4",
        "statement": "A's PARSEC closure ratios rise by 8-20 %.",
        "measured_quantity": (
            "closure_ratio(v8)/closure_ratio(v7) - 1 for the nine CygOB2-A "
            "PARSEC cells (3 R_V x 3 alpha) of wp6_closure"
        ),
        "pass_if": (
            "0.08 <= rise <= 0.20 in all nine cells; the median rise is "
            "reported but does not gate"
        ),
    },
    {
        "id": "P5",
        "statement": (
            "The association closing slope stays within [2.10, 2.30]."
        ),
        "measured_quantity": (
            "closing_slope_grid_median_by_subgroup['association'] from the "
            "repair_v8 WP12 closure-slope run (repair_v7: 2.1966)"
        ),
        "pass_if": "2.10 <= value <= 2.30",
    },
    {
        "id": "P6",
        "statement": (
            "The repair_v5 retained age envelope is 2.00-4.01 Myr and no PMS "
            "row is retained."
        ),
        "measured_quantity": (
            "rows of wp4_age_posteriors_repair_v5 with measurable and not "
            "grid_railed: min and max age_map rounded to 0.01 Myr, and the "
            "count with indicator == 'pms'"
        ),
        "pass_if": "min == 2.00 AND max == 4.01 AND pms_count == 0",
        "note": (
            "Not a forecast: the diagnostic already measured it (disclosed "
            "below).  Registered so that the macro which replaces the "
            "hard-coded envelope is checked against a value fixed in advance."
        ),
    },
]

# Expectations from brief section 3.  Reported before/after; they do NOT gate
# anything and are not scored PASS/FAIL.
REPORT_ONLY_EXPECTATIONS = {
    "k_A_PARSEC_relative_change": "-1 % to -3 %",
    "A_N_death_baseline": "about 4.05-4.15 (repair_v7 4.17)",
    "A_census_above_8_PARSEC_rv3.1": "about 65 (repair_v7 57 thresholded)",
    "A_closing_slope_grid_median": "lower than 2.34, about 2.2-2.3",
    "WP5_gate": "40/54 cells and 15/36 all-passing headline branches may change",
}

# -------------------------------------------------------- integrity checks
# These are code-correctness checks, not science predictions.  A failure stops
# the chain until the code is fixed; it is never resolved by moving a number.
INTEGRITY_CHECKS = [
    {
        "id": "I1",
        "check": (
            "Every non-anchor row of wp4_mass_posteriors_repair_v8 (and its "
            "posterior sample cube) is identical to repair_v5 on all six "
            "branches.  The photometric path is untouched by repair_v8."
        ),
    },
    {
        "id": "I2",
        "check": (
            "The repair_v8 R_V = 3.1 anchor masses reproduce "
            "tables/issue19_anchor_mass_rederivation.csv (the scope "
            "diagnostic, same procedure at the same ages) exactly."
        ),
    },
    {
        "id": "I3",
        "check": (
            "Replay: the modified wp5_fit_imf_joint.py, fed the repair_v5 "
            "masses and the repair_v7 node responses, reproduces "
            "wp5_imf_normalization_repair_v7 (k quantiles, gate columns) "
            "exactly.  Guards the version plumbing, not the science."
        ),
    },
    {
        "id": "I4",
        "check": (
            "Every repair_v7-chain artifact hashed below under "
            "compare_artifacts still has the same SHA-256 after the rerun."
        ),
    },
]

ADOPTION_RULE = (
    "repair_v8 corrects a defect; it is not a model choice.  It becomes the "
    "chain the manuscript quotes if and only if I1-I4 pass.  Adoption does "
    "NOT depend on P1-P6: a failed prediction is recorded as failed and means "
    "the defect's effect was misunderstood, not that the fix is wrong.  If a "
    "downstream gate (WP5 G3, WP6, WP7, WP9, WP12 R1-R7) changes verdict on "
    "repair_v8, the new verdict is reported as it stands and opened as an "
    "issue in PROJECT_TRACE section 9; repair_v7 is not restored and no "
    "threshold is moved to recover the old verdict."
)

DISCLOSED_PRIOR_KNOWLEDGE = [
    "tables/issue19_anchor_mass_rederivation.csv (2026-09-30/10-01): at the "
    "R_V = 3.1 repair_v5 ages, 8 A anchors (sum of membership probability "
    "8.0) cross 2-8 -> >8 Msun on PARSEC; no other crossing of 2 or 8 Msun.  "
    "P1 is built on this.  The per-R_V (F3(b)) masses have NOT been computed; "
    "P1 and P2 at R_V = 3.0 and 3.5 are therefore genuine forecasts.",
    "tables/issue19_envelope_by_version.csv: repair_v5 retains 66 rows, 0 PMS, "
    "MAP span 2.00-4.01 Myr.  P6 restates it.",
    "The repair_v7 values quoted in each pass_if were read from the repair_v7 "
    "products hashed below.",
    "No repair_v8 product, and no WP5/WP6/WP7 quantity computed from "
    "corrected anchor masses, existed when this file was written.",
]


def _node_files() -> list[str]:
    pattern = [
        f"wp5_agenode_*_{RESPONSE_VERSION}_response.parquet",
        f"wp5_agenode_*_{RESPONSE_VERSION}_curve.parquet",
        f"wp6_massext_*_{RESPONSE_VERSION}_response.parquet",
        f"wp6_massext_*_{RESPONSE_VERSION}_curve.parquet",
    ]
    files = []
    for glob in pattern:
        files += sorted(str(p.relative_to(w.ROOT)) for p in w.PROC.glob(glob))
    return files


CONSUMED_INPUTS = [
    "data/processed/wp1_spectroscopic_anchors.parquet",
    "data/processed/wp2_members.parquet",
    "tables/wp2_subgroup_labels.parquet",
    f"data/processed/wp3_extinction_{WP3_WP4_VERSION}.parquet",
    f"data/processed/wp3_extinction_posterior_{WP3_WP4_VERSION}.npz",
    f"data/processed/wp4_age_posteriors_{WP3_WP4_VERSION}.parquet",
    "data/processed/wp3_isochrones_parsec.parquet",
    "data/processed/wp3_isochrones_mist.parquet",
    f"provenance/wp6_massive_injections_execution_{RESPONSE_VERSION}.json",
    "tables/wp6_runaways.csv",
    "tables/wp6_runaway_crossmatch.csv",
    "provenance/wp6_runaways_execution.json",
]

COMPARE_ARTIFACTS = [
    f"data/processed/wp4_mass_posteriors_{WP3_WP4_VERSION}.parquet",
    f"data/processed/wp4_mass_posterior_samples_{WP3_WP4_VERSION}.npz",
    "data/processed/wp4_anchor_hrd.parquet",
    "data/processed/wp4_age_posteriors.parquet",
    f"data/processed/wp5_imf_normalization_{COMPARE_VERSION}.parquet",
    f"data/processed/wp5_imf_posterior_draws_{COMPARE_VERSION}.npz",
    f"data/processed/wp5_association_mass_{COMPARE_VERSION}.parquet",
    f"data/processed/wp5_mass_function_bins_{COMPARE_VERSION}.parquet",
    f"data/processed/wp5_completeness_curves_{COMPARE_VERSION}.parquet",
    f"provenance/wp5_imf_fit_execution_{COMPARE_VERSION}.json",
    "tables/wp6_massive_census.csv",
    "tables/wp6_orphan_anchors.csv",
    f"tables/wp6_closure_{COMPARE_VERSION}.csv",
    f"tables/wp6_closure_attribution_{COMPARE_VERSION}.csv",
    "tables/wp6_massive_census.cat",
    "tables/wp7_ledger.csv",
    "tables/wp7_rsn_curves.csv",
    "tables/wp7_age_sensitivity.csv",
    "tables/wp7_bh_threshold_scan.csv",
    "tables/wp7_convergence.csv",
    "tables/wp7_alpha_headline_branch_sets.csv",
    "tables/wp7_binary_bound.csv",
    "tables/wp7_binary_bound_branches.csv",
    "tables/wp8_crosschecks.csv",
    "tables/wp8_tension_list.csv",
    "tables/wp9_verdict.csv",
    "tables/wp9_sensitivity.csv",
    "tables/wp12_wp5_gate_map.csv",
    "tables/wp12_combination_gate.csv",
    "tables/wp12_branch_gate_table.csv",
    "tables/wp12_closure_by_alpha.csv",
    "tables/wp12_closing_slopes.csv",
    "tables/wp12_mixed_slope_ledger.csv",
    "tables/wp12_c4_scan.csv",
    "tables/wp12_c3_subtype.csv",
    "tables/wp12_scenario_score.csv",
    "tables/wp4_wp5_age_reconciliation.csv",
    "tables/wp5_association_mass_reconciliation.csv",
    "provenance/wp12_revision_prereg.json",
    "provenance/wp12_closure_slopes_execution.json",
    "provenance/wp7_ledger_execution.json",
    "provenance/wp9_verdict_execution.json",
    "provenance/wp6_closure_test_execution_repair_v7.json",
    "provenance/wp6_ledger_execution.json",
    "provenance/wp6_massive_census_execution.json",
    "provenance/wp4_anchors_hrd_execution.json",
    "provenance/wp4_mass_repair_execution.json",
    "tables/issue19_anchor_mass_rederivation.csv",
    "tables/issue19_envelope_by_version.csv",
]


def hashed(paths: list[str]) -> dict:
    out = {}
    for rel in paths:
        path = w.ROOT / rel
        out[rel] = w.sha256(path) if path.exists() else None
    return out


def wp12_frozen_still_intact() -> dict:
    """The WP12 preregistration's own hashes, re-checked now."""
    record = json.loads((w.PROVENANCE / "wp12_revision_prereg.json").read_text())
    moved = {}
    for name, entry in record["frozen_inputs"].items():
        path = w.ROOT / entry["path"]
        digest = w.sha256(path) if path.exists() else None
        if digest != entry["sha256"]:
            moved[name] = {"path": entry["path"], "now": digest}
    return {"checked": len(record["frozen_inputs"]), "moved": moved,
            "intact": not moved}


def main() -> None:
    if OUT.exists():
        raise SystemExit(
            f"{OUT.relative_to(w.ROOT)} already exists -- a preregistration is "
            "written once and never regenerated."
        )
    existing_v8 = sorted(
        str(p.relative_to(w.ROOT))
        for p in list(w.PROC.glob(f"*{NEW_VERSION}*"))
        + list(w.TABLES.glob(f"*{NEW_VERSION}*"))
    )
    if existing_v8:
        raise SystemExit(
            "repair_v8 products already exist; a preregistration written after "
            f"them would be worthless: {existing_v8[:5]}"
        )
    nodes = _node_files()
    if len([n for n in nodes if "agenode" in n and n.endswith("_response.parquet")]) != 162:
        raise SystemExit("expected 162 repair_v7 WP5 node responses")
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/issue19_prereg.py",
        "issue": "#19 -- stale pre-repair WP4 inputs",
        "brief": "tasks/issue19_stale_wp4_inputs_brief.md",
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "versions": {
            "new": NEW_VERSION,
            "wp3_extinction_and_wp4_ages": WP3_WP4_VERSION,
            "injection_responses_reused": RESPONSE_VERSION,
            "compared_against": COMPARE_VERSION,
        },
        "decisions": DECISIONS,
        "change": CHANGE,
        "reuse": REUSE,
        "predictions": PREDICTIONS,
        "report_only_expectations": REPORT_ONLY_EXPECTATIONS,
        "integrity_checks": INTEGRITY_CHECKS,
        "adoption_rule": ADOPTION_RULE,
        "failed_predictions_rule": (
            "Failed predictions stay recorded as failed.  This file is not "
            "amended after the rerun; corrections, if any, go in the outcome "
            "record with the reason."
        ),
        "disclosed_prior_knowledge": DISCLOSED_PRIOR_KNOWLEDGE,
        "consumed_inputs": hashed(CONSUMED_INPUTS),
        "reused_injection_responses": hashed(nodes),
        "compare_artifacts": hashed(COMPARE_ARTIFACTS),
        "wp12_frozen_inputs_recheck": wp12_frozen_still_intact(),
    }
    missing = [k for k, v in {**record["consumed_inputs"],
                              **record["compare_artifacts"]}.items() if v is None]
    if missing:
        raise SystemExit(f"inputs missing at preregistration: {missing}")
    w.write_json(OUT, record)
    print(f"wrote {OUT.relative_to(w.ROOT)}")
    print(f"  consumed inputs hashed : {len(record['consumed_inputs'])}")
    print(f"  reused node files      : {len(nodes)}")
    print(f"  compare artifacts      : {len(record['compare_artifacts'])}")
    print(f"  WP12 frozen inputs intact: {record['wp12_frozen_inputs_recheck']['intact']}")


if __name__ == "__main__":
    main()
