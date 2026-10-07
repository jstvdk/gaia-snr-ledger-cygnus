#!/usr/bin/env python3
"""repair_v9 -- pre-registration of the chain and the age scan.

Written and executed BEFORE any repair_v9 chain product exists.  It records
the SHA-256 of every input the chain consumes, of every repair_v8 artifact
the checks and predictions are scored against, and of every script that runs
or scores the chain; it fixes the integrity checks I1-I5 (with tolerances),
the predictions V1-V5 (with the rule that scores each), the age scan and the
adoption rule.  Tolerances and thresholds are imported from the scoring code
itself (repair_v9_integrity.TOL, repair_v9_score constants), so the record and
the code cannot disagree; both files are hashed here.

It refuses to overwrite itself and refuses to run if any repair_v9 chain
output already exists.

Basis: tasks/repair_v9_prereg_proposal.md (signed off by the owner on
2026-10-07, see owner_signoff), brief §4 and §4.3-4.4
(tasks/stage_repair_v9_wp4_age_redesign_and_wp13.md), owner decisions 6-8
(provenance/decisions_2026_10_07.json).

Run:  PYTHONPATH=scripts python3 scripts/repair_v9_prereg.py
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import repair_v9_integrity as RI
import repair_v9_score as RS
import wp5_common as w
from wp6_mass_extension_decision import IMF_UPPER_LIMIT, turnoff_mass

OUT = w.PROVENANCE / "repair_v9_prereg.json"

INPUTS = {
    # ages
    "headline_ages": "data/processed/wp4_age_posteriors_repair_v9_headline.parquet",
    "stage1_ages": "data/processed/wp4_age_posteriors_repair_v9.parquet",
    "replay_ages_v5": "data/processed/wp4_age_posteriors_repair_v5.parquet",
    "headline_ages_record": "provenance/repair_v9_headline_ages_execution.json",
    "stage1_outcome": "provenance/wp4v9_age_outcome.json",
    "decisions_2026_10_07": "provenance/decisions_2026_10_07.json",
    # unchanged upstream
    "wp1_gaia_narrow": "data/processed/wp1_gaia_narrow.parquet",
    "wp1_2mass_join": "data/processed/wp1_2mass_join.parquet",
    "wp1_spectroscopic_anchors": "data/processed/wp1_spectroscopic_anchors.parquet",
    "wp2_members": "data/processed/wp2_members.parquet",
    "wp2_subgroup_labels": "tables/wp2_subgroup_labels.parquet",
    "wp2_membership_manifest": "provenance/wp2_membership_manifest.json",
    "wp3_extinction_v5": "data/processed/wp3_extinction_repair_v5.parquet",
    "wp3_extinction_posterior_v5": "data/processed/wp3_extinction_posterior_repair_v5.npz",
    "wp3_repair_execution": "provenance/wp3_repair_execution.json",
    "isochrones_parsec": "data/processed/wp3_isochrones_parsec.parquet",
    "isochrones_mist": "data/processed/wp3_isochrones_mist.parquet",
    # repair_v8 comparison artifacts (I1, V1-V2, before/after)
    "v8_anchors": "data/processed/wp4_anchor_hrd_repair_v8.parquet",
    "v8_masses": "data/processed/wp4_mass_posteriors_repair_v8.parquet",
    "v8_mass_samples": "data/processed/wp4_mass_posterior_samples_repair_v8.npz",
    "v8_normalization": "data/processed/wp5_imf_normalization_repair_v8.parquet",
    "v8_ledger": "tables/wp7_ledger_repair_v8.csv",
    "v8_closure": "tables/wp6_closure_repair_v8.csv",
    "wp12_original_prereg": "provenance/wp12_revision_prereg.json",
}

CODE = [
    # new
    "scripts/repair_v9_headline_ages.py", "scripts/repair_v9_injections.py",
    "scripts/repair_v9_age_scan.py", "scripts/repair_v9_integrity.py",
    "scripts/repair_v9_score.py", "scripts/repair_v9_prereg.py",
    "scripts/run_repair_v9_chain.sh",
    # changed for repair_v9 (each change inert for every earlier version; I1 tests it)
    "scripts/chain.py", "scripts/wp5_joint_age_fit.py",
    "scripts/wp4_mass_posteriors_repair.py", "scripts/wp4_anchors_hrd.py",
    "scripts/wp12_prereg.py",
    # run unchanged
    "scripts/wp4v9_common.py", "scripts/wp4_repair_common.py", "scripts/wp5_common.py",
    "scripts/wp5_injections_repair.py", "scripts/wp5_fit_imf_joint.py",
    "scripts/wp5_fbin_discriminator_prereg.py", "scripts/wp6_mass_extension_decision.py",
    "scripts/wp6_massive_injections.py", "scripts/wp6_massive_census.py",
    "scripts/wp6_closure_test.py", "scripts/wp6_closure_attribution.py", "scripts/wp6_ledger.py",
    "scripts/wp7_ledger.py", "scripts/wp7_ledger_prereg.py", "scripts/wp7_convergence_scan.py",
    "scripts/wp7_alpha_headline_adopt.py", "scripts/wp7_binary_bound.py",
    "scripts/wp8_crosschecks.py", "scripts/wp9_verdict.py", "scripts/wp11_isotope_forecast.py",
    "scripts/wp4_wp5_age_reconciliation.py", "scripts/wp5_association_mass_reconciliation.py",
    "scripts/wp5_alpha_plausibility.py", "scripts/wp12_common.py", "scripts/wp12_gate_landscape.py",
    "scripts/wp12_closure_slopes.py", "scripts/wp12_scenario_score.py",
    "scripts/wp12_neighbour_budget.py", "scripts/wp10_inputs.py",
]

OUTPUTS_THAT_MUST_NOT_EXIST = [
    "data/processed/wp4_anchor_hrd_repair_v9.parquet",
    "data/processed/wp4_anchor_hrd_repair_v9_replay.parquet",
    "data/processed/wp4_mass_posteriors_repair_v9.parquet",
    "data/processed/wp4_mass_posteriors_repair_v9_replay.parquet",
    "data/processed/wp5_imf_normalization_repair_v9.parquet",
    "data/processed/wp5_imf_normalization_repair_v9_replay.parquet",
    "provenance/wp5_injections_agenodes_execution_repair_v9.json",
    "provenance/wp6_massive_injections_execution_repair_v9.json",
    "tables/wp6_closure_repair_v9.csv", "tables/wp7_ledger_repair_v9.csv",
    "tables/wp7_ledger_repair_v9_replay.csv", "tables/repair_v9_age_scan.csv",
    "data/processed/wp4_age_posteriors_repair_v9_scan00.parquet",
    "provenance/repair_v9_integrity.json", "provenance/repair_v9_outcome.json",
]


def v2_turnoff_terms() -> dict:
    """The sign V2 predicts, derived from wp6_closure_test before any run.

    closure ratio = observed / (k * integral[4, M_to] M^-alpha R(M) dM).  An age
    change reaches it three ways: M_to (the cap of the integral), the refitted
    k, and 'observed' (the census, from photometric masses at the new age).
    Only the first has a sign fixed by the code: a younger age raises M_to,
    enlarges the predicted count and LOWERS the ratio.  Its size is computed
    here with R = 1 above 8 Msun (step response) at alpha 2.3, from the repair_v8
    WP5 truth-age posterior means to the repair_v9 headline MAPs."""
    n8 = pd.read_parquet(w.ROOT / INPUTS["v8_normalization"])
    head = pd.read_parquet(w.ROOT / INPUTS["headline_ages"])
    a = 2.3
    integral = lambda m: (8 ** (1 - a) - min(m, IMF_UPPER_LIMIT) ** (1 - a)) / (a - 1)  # noqa: E731
    out = {}
    for fam in w.FAMILIES:
        for sg in w.SUBGROUPS:
            t8 = float(n8[n8.subgroup.eq(sg) & n8.family.eq(fam) & np.isclose(n8.R_V, 3.1)
                          & np.isclose(n8.alpha, a)].truth_age_posterior_mean_Myr.iloc[0])
            t9 = float(head[head.subgroup.eq(sg) & head.family.eq(fam) & np.isclose(head.R_V, 3.1)
                            & np.isclose(head.f_bin, 0.4) & np.isclose(head.dmu, 0.0)].age_map.iloc[0])
            m8, m9 = min(turnoff_mass(fam, t8), IMF_UPPER_LIMIT), min(turnoff_mass(fam, t9), IMF_UPPER_LIMIT)
            out[f"{sg}|{fam}"] = {"age_v8_Myr": round(t8, 3), "age_v9_Myr": round(t9, 3),
                                  "turnoff_v8_Msun": round(m8, 1), "turnoff_v9_Msun": round(m9, 1),
                                  "ratio_factor_turnoff_only": round(integral(m8) / integral(m9), 4),
                                  "predicted_direction": "down" if m9 > m8 else "up"}
    return out


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"{OUT.relative_to(w.ROOT)} exists; a pre-registration is written once")
    present = [p for p in OUTPUTS_THAT_MUST_NOT_EXIST if (w.ROOT / p).exists()]
    if present:
        raise SystemExit(f"repair_v9 outputs already exist, cannot pre-register: {present}")
    head_rec = json.loads((w.ROOT / INPUTS["headline_ages_record"]).read_text())
    v2 = v2_turnoff_terms()

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/repair_v9_prereg.py",
        "issue": "#21 (WP4 extinction-error model) and #20 (is C young?)",
        "brief": "tasks/stage_repair_v9_wp4_age_redesign_and_wp13.md §4",
        "proposal": {"path": "tasks/repair_v9_prereg_proposal.md",
                     "sha256": w.sha256(w.ROOT / "tasks/repair_v9_prereg_proposal.md")},
        "owner_signoff": {
            "date": "2026-10-07 (about 18:39 CEST)",
            "verbatim": "proceed with execution and document all decisions and results",
            "logged_in": "AI_PROMPTS.md",
            "what_was_signed": (
                "the proposal as written, plus two additions shown to the owner before sign-off: "
                "(1) V2's expected signs derived from wp6_closure_test (below), with the "
                "disclosure that the turnoff term is only 2-5 %; (2) I1 compares across "
                "platforms, so its tolerance is Monte Carlo / floating-point level, not exact"),
        },
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "platform_note": (
            "repair_v8 and Stage 1 were produced on macOS arm64; repair_v9 runs on the owner's "
            "AlmaLinux 9.6 x86_64 machine with the same pinned package versions (python 3.11.15, "
            "numpy 2.4.6, pandas 3.0.3, scipy 1.17.1, scikit-learn 1.9.0).  Floating-point results "
            "can differ in the last bits; INSTALL.md: such differences are not findings."),
        "versions": {
            "chain": "repair_v9 (scripts/chain.py; ADOPTED stays repair_v8 until the rule below passes)",
            "wp3_extinction": "repair_v5 (unchanged)",
            "wp4_ages": "repair_v9_headline (new)",
            "wp4_masses": "repair_v9 (regenerated at the new ages)",
            "wp5": "repair_v9 (regenerated)",
            "responses": "repair_v9 (WP5 nodes and WP6 extension regenerated; repair_v7 truth model: "
                         "interpolated nodes, mass-dependent f_bin 'extended')",
        },
        "headline_ages": {
            "rule": head_rec["rule"],
            "decision_cell_R_V_3p1": head_rec["decision_cell_R_V_3p1"],
            "B_flag": "not measured for B: borrows the A + C combined spectroscopic posterior",
            "built_before_this_record": (
                "scripts/repair_v9_headline_ages.py ran before this pre-registration, as the handoff "
                "orders (step 3).  It is an INPUT (a deterministic combination of frozen Stage 1 "
                "curves), contains no chain result, and is hashed here."),
        },
        "code_changes": {
            "scripts/chain.py": "repair_v9, repair_v9_replay and repair_v9_scan00..10 registered; ADOPTED unchanged",
            "scripts/wp5_joint_age_fit.py": "repair_v9 and repair_v9_replay added to AGE_INTERPOLATED_VERSIONS",
            "scripts/wp4_mass_posteriors_repair.py": (
                "--age-version (default WP_REPAIR_VERSION, i.e. unchanged).  The 'stale mass_method "
                "label' of the proposal was already fixed in repair_v8 (issue19_prereg); no change."),
            "scripts/wp4_anchors_hrd.py": "--preregistration (provenance field only; default unchanged)",
            "scripts/wp12_prereg.py": (
                "--chain-record re-points wp4_age_posteriors to the chain-declared age product "
                "(for repair_v8 that is repair_v5, i.e. unchanged)"),
            "scripts/repair_v9_injections.py": (
                "new parallel, resumable driver around the unchanged inject_curve; fresh "
                "default_rng(SEED) per node, so order and worker do not matter (tested by I1c)"),
        },
        "steps": "scripts/run_repair_v9_chain.sh, phases A-D (anchors+masses, injections, fits+ledgers, WP6-WP12)",
        "iteration_counts": {"wp7_ledger": 2_000_000, "wp11": 500_000, "age_scan_ledger": 200_000},
        "age_scan": {
            "decision": "owner decision 7 of 2026-10-07",
            "grid_Myr": {"PARSEC": [2.00, 2.24, 2.51, 2.82, 3.16, 3.55, 3.98, 4.47, 5.01, 5.62, 6.31],
                         "MIST": [2.00, 2.25, 2.52, 2.83, 3.18, 3.57, 4.01, 4.50, 5.05, 5.67, 6.37]},
            "rule": ("at scan point i all three subgroups are set to the i-th native age of their family "
                     "(point-mass age table); anchors and photometric masses at t; ONE native node per "
                     "subgroup x family x R_V (snapped, no interpolation, never a reused gate-G2 snapshot); "
                     "WP5 fit (k refitted at t); WP7 engine per subgroup and all 54 branches at 200,000 "
                     "iterations, all-explode, storing the full N_death distribution"),
            "combination": ("subgroup populations are independent Poisson draws, so N_death for any "
                            "(t_A, t_B, t_C) is the convolution of the stored per-subgroup distributions"),
            "status": "sensitivity, never the headline; no WP6, no WP8-12",
        },
        "integrity_checks": {
            "tolerances": RI.TOL,
            "I1": ("replay -- the repair_v9 code path fed the repair_v5 ages reproduces repair_v8.  Scored "
                   "per part as PASS_EXACT, PASS_WITHIN_TOLERANCE (tolerances above) or FAIL: "
                   "I1a anchors (max |diff| any numeric column <= I1a_anchor_abs, same rows and NaNs); "
                   "I1b masses (<= I1b_mass_frac_beyond of entries beyond I1b_mass_rel relative; "
                   "baseline sum P(M>8) within I1b_expected_gt8_abs); "
                   "I1c injector (the 27 PARSEC R_V 3.1 repair_v7 nodes regenerated: >= I1c_draw_frac_equal "
                   "of recovered-mass draws equal to I1c_draw_rel per node, isotonic curve within "
                   "I1c_curve_abs); I1d WP5 (k_median within I1d_k_rel in all 54 cells, residual gate "
                   "differing in <= I1d_gate_flips_max cells); I1e WP7 (baseline N_death mean, association "
                   "and each subgroup, within I1e_ndeath_abs of repair_v8 at 2,000,000 iterations)."),
            "I2": ("anchor masses are per family x R_V branch and each anchor's age_used equals the native "
                   "age nearest its subgroup's repair_v9 headline MAP (unlabelled anchors: the median of "
                   "the three subgroup MAPs) -- zero mismatches"),
            "I3": ("every consumed input's SHA-256 matches this record (no exception); every frozen_code "
                   "hash matches, or the change is listed in provenance/repair_v9_deviations.json with its "
                   "new hash and reason"),
            "I4": ("every repair_v9 execution record references only the headline age table (replay: "
                   "repair_v5; scan: its own scan table) and no wp10_inputs.FORBIDDEN artifact; the new "
                   "scripts contain no repair_v5, repair_v8 or unversioned age-table literal"),
            "I5": ("for each subgroup, PARSEC and MIST, R_V 3.1, alpha 2.3, coeval: the scan's mean N_death "
                   "interpolated linearly in age and averaged over the repair_v9 WP5 truth-age draws "
                   "reproduces the repair_v9 ledger mean within max(I5_abs, I5_rel x headline)"),
        },
        "predictions": {
            "basis": ("estimates made before the run from repair_v8's k at fixed MAP ages (PARSEC, R_V 3.1, "
                      "alpha 2.3): A ~0.9 at 3.16 Myr, B ~2.2, C ~2.6, total ~5.7; averaging over the "
                      "posterior widths raises these because N_death is convex in age"),
            "V1a": "PARSEC baseline association N_death mean < repair_v8's (read from tables/wp7_ledger_repair_v8.csv)",
            "V1b": f"PARSEC baseline association N_death mean in {list(RS.V1B_RANGE)}",
            "V1c": "baseline subgroup means: A falls, B falls, C rises (vs repair_v8)",
            "V2": {
                "statement": ("baseline-cell (PARSEC, R_V 3.1, alpha 2.3) closure ratio moves A down, B down, "
                              "C up vs repair_v8; MIST reported, not scored"),
                "derivation": v2_turnoff_terms.__doc__,
                "turnoff_term_by_cell": v2,
                "disclosure": ("the turnoff term alone is a 2-5 % effect; the refitted k and the new census "
                               "masses can override it, so V2 is a weak prediction"),
            },
            "V3": f"C baseline N_death mean > {RS.V3_MIN} (C is NOT below the first-death boundary)",
            "V4": (f"scan: mean N_death(t) non-decreasing in t for every subgroup x 54 branches (a drop "
                   f"larger than max({RS.V4_DROP_ABS}, {RS.V4_DROP_REL:.0%} of the previous point) is a "
                   f"violation), and < {RS.V4_ZERO_MAX} at t <= {RS.V4_YOUNG_MAX_MYR} Myr on the coeval "
                   "branches (operationalises the proposal's 'is 0 for t <= 2.5 Myr': with a 1-2 Myr "
                   "birth spread, or MIST's 108 Msun turnoff at 2.52 Myr, it is not exactly 0)"),
            "V5": (f"scan: at every t >= {RS.V5_FROM_MYR} Myr, alpha 2.3, coeval, each family x R_V, the "
                   f"three subgroups' mean N_death agree within max/min <= {RS.V5_RATIO_MAX}"),
            "scoring": "scripts/repair_v9_score.py, as written; predictions never gate adoption",
        },
        "adoption_rule": {
            "rule": "ADOPTED = 'repair_v9' in scripts/chain.py iff I1a-I1e, I2, I3, I4 and I5 all pass",
            "then": ["reports/repair_v9_completion_report.md with a before/after table for every headline "
                     "number", "close issue #21, update issue #20", "update CLAUDE.md headline numbers and "
                     "the withdrawn table",
                     "manuscript plumbing (wp10_inputs, numbers.tex, tables, figures, wp10_validate)"],
            "ordering_note": ("the proposal lists manuscript plumbing as chain step 7; it is executed only "
                              "AFTER adoption because it changes what the manuscript quotes, and no repair_v9 "
                              "number may be quoted before adoption (handoff §7)"),
            "if_not_adopted": "repair_v8 stays the quoted chain; the failure is reported and goes to the owner",
        },
        "failed_predictions_rule": "a failed prediction stays failed and is reported as such",
        "disclosed_prior_knowledge": [
            "Stage 1 (reports/wp4v9_age_redesign.md): GA1 failed in all six cells; the Gaia CMD does not "
            "measure these ages; A and C carry test c, C coeval with A (W1, W2 true)",
            "the headline ages were computed before this record (A 3.16, B 3.55, C 3.55 Myr PARSEC; "
            "A 3.57, B 3.57, C 4.01 MIST)",
            "the rough per-subgroup estimates in predictions.basis, and that deaths per subgroup go from 0 "
            "at <= 2.5 Myr to ~1 at 3.2, ~4 at 4.0 and ~8 at 5.0 Myr (handoff §4)",
            "a 1-node timing test of inject_curve at the headline ages (A, PARSEC, R_V 3.1, 3.162 Myr) was "
            "run to size the worker pool; its output was not written to disk and not inspected beyond its "
            "shape and runtime",
            "no repair_v9 chain product (anchors, masses, nodes, fits, ledgers, scan) existed when this was written",
        ],
        "consumed_inputs": {name: {"path": rel, "sha256": w.sha256(w.ROOT / rel)}
                            for name, rel in INPUTS.items()},
        "frozen_code": {rel: w.sha256(w.ROOT / rel) for rel in CODE},
    }
    w.write_json(OUT, record)
    print(f"wrote {OUT.relative_to(w.ROOT)}  sha256 {w.sha256(OUT)}")
    for k, v in v2.items():
        print(f"  V2 {k}: {v['age_v8_Myr']} -> {v['age_v9_Myr']} Myr, M_to {v['turnoff_v8_Msun']} -> "
              f"{v['turnoff_v9_Msun']}, ratio x{v['ratio_factor_turnoff_only']} ({v['predicted_direction']})")


if __name__ == "__main__":
    main()
