#!/usr/bin/env python3
"""WP13 -- pre-registration of the pooled-versus-resolved ablation (M0 vs M1).

Transcribes tasks/wp13_prereg_proposal.md as signed off by the owner on
2026-10-08 (provenance/decisions_2026_10_08.json: D1(a), D2, D3 = 18 scored /
36 reported, D4, D5, D6), the brief's §3-§6 and §9
(tasks/wp13_pooled_vs_resolved_ablation_brief.md), and the clarifications
A1-A6 below, which the proposal left open and which were put to the owner
before this record was written (ADDENDUM_SIGNOFF).

Written and executed BEFORE any M0-side number exists -- T1 included, which is
computable in seconds from stored curves and deliberately is not computed
here.  It records the SHA-256 of every consumed input, of every script that
runs or scores WP13, and of wp13_rules.py, which holds every threshold and
scoring rule as code (the scorer imports them from there).

Refuses to overwrite itself, refuses to run if any WP13 output exists, and
refuses to run before the WP13 driver and scorer exist (they are hashed here).

Run:  PYTHONPATH=scripts python3 scripts/wp13_prereg.py
"""
from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone

import wp13_rules as R
import wp5_common as w

OUT = w.PROVENANCE / "wp13_prereg.json"

ADDENDUM_SIGNOFF = (
    "2026-10-08 (~17:30 CEST), four multiple-choice answers in the session: A1 'A+C on both sides'; "
    "A2 'expected count mu'; A3/A4 'accept as proposed'; A5/A6 'yes, that order' (WP13 code written "
    "and plan-tested before this record, then the owner approves the commit, then M0 runs)")

INPUTS = {
    # M1 (repair_v9, hash-verified, reproduced not recomputed)
    "headline_ages": "data/processed/wp4_age_posteriors_repair_v9_headline.parquet",
    "headline_ages_record": "provenance/repair_v9_headline_ages_execution.json",
    "stage1_ages": "data/processed/wp4_age_posteriors_repair_v9.parquet",
    "testc_curves": "provenance/issue20b_run_execution.json",
    "testc_table": "tables/issue20b_age_tests.csv",
    "native_ages": "provenance/wp4v9_fit_execution.json",
    "m1_anchors": "data/processed/wp4_anchor_hrd_repair_v9.parquet",
    "m1_masses": "data/processed/wp4_mass_posteriors_repair_v9.parquet",
    "m1_mass_samples": "data/processed/wp4_mass_posterior_samples_repair_v9.npz",
    "m1_injections_record": "provenance/wp5_injections_agenodes_execution_repair_v9.json",
    "m1_normalization": "data/processed/wp5_imf_normalization_repair_v9.parquet",
    "m1_draws": "data/processed/wp5_imf_posterior_draws_repair_v9.npz",
    "m1_ledger": "tables/wp7_ledger_repair_v9.csv",
    "m1_age_scan": "tables/wp7_age_sensitivity_repair_v9.csv",
    "m1_prereg": "provenance/repair_v9_prereg.json",
    "m1_integrity": "provenance/repair_v9_integrity.json",
    # held fixed between M0 and M1
    "wp1_spectroscopic_anchors": "data/processed/wp1_spectroscopic_anchors.parquet",
    "wp2_members": "data/processed/wp2_members.parquet",
    "wp2_subgroup_labels": "tables/wp2_subgroup_labels.parquet",
    "wp3_extinction_v5": "data/processed/wp3_extinction_repair_v5.parquet",
    "wp3_extinction_posterior_v5": "data/processed/wp3_extinction_posterior_repair_v5.npz",
    "isochrones_parsec": "data/processed/wp3_isochrones_parsec.parquet",
    "isochrones_mist": "data/processed/wp3_isochrones_mist.parquet",
    # basis and disclosure
    "brief": "tasks/wp13_pooled_vs_resolved_ablation_brief.md",
    "proposal": "tasks/wp13_prereg_proposal.md",
    "decisions_2026_10_08": "provenance/decisions_2026_10_08.json",
    "brief_rebase_v9": "tables/wp13_brief_rebase_repair_v9.csv",
}

CODE = [
    # WP13 (new)
    "scripts/wp13_prereg.py", "scripts/wp13_rules.py", "scripts/wp13_m0.py",
    "scripts/wp13_score.py", "scripts/wp13_t7.py",
    # engines, run unchanged apart from the inert pooled-label hook (W13-I1 tests it)
    "scripts/chain.py", "scripts/wp4_common.py", "scripts/wp5_common.py",
    "scripts/wp4v9_common.py", "scripts/issue20b_common.py", "scripts/issue20b_run.py",
    "scripts/wp4_anchors_hrd.py", "scripts/wp4_mass_posteriors_repair.py",
    "scripts/wp4_repair_common.py", "scripts/wp5_injections_repair.py",
    "scripts/repair_v9_injections.py", "scripts/wp5_joint_age_fit.py",
    "scripts/wp5_fit_imf_joint.py", "scripts/wp7_ledger.py", "scripts/wp7_ledger_prereg.py",
    "scripts/wp6_mass_extension_decision.py", "scripts/repair_v9_integrity.py",
]

OUTPUTS_THAT_MUST_NOT_EXIST = [
    "data/processed/wp4_age_posteriors_wp13_m0.parquet",
    "data/processed/wp4_anchor_hrd_wp13_m0.parquet",
    "data/processed/wp4_mass_posteriors_wp13_m0.parquet",
    "provenance/wp5_injections_agenodes_execution_wp13_m0.json",
    "data/processed/wp5_imf_normalization_wp13_m0.parquet",
    "data/processed/wp5_imf_posterior_draws_wp13_m0.npz",
    "provenance/wp13_m0_ages_execution.json",
    "data/processed/wp4_anchor_hrd_wp13_identity.parquet",
    "data/processed/wp4_mass_posteriors_wp13_identity.parquet",
    "data/processed/wp5_imf_normalization_wp13_identity.parquet",
    "tables/wp13_t1.csv", "tables/wp13_ablation.csv", "tables/wp13_t7.csv", "tables/wp13_d6_grid.csv",
    "data/processed/wp13_baseline_curves.npz",
    "provenance/wp13_ledger_execution.json", "provenance/wp13_t7_execution.json",
    "provenance/wp13_integrity.json", "provenance/wp13_outcome.json",
    "reports/wp13_ablation.md",
]

RUN_ORDER = [
    "PYTHONPATH=scripts python3 scripts/wp13_score.py I3",
    "PYTHONPATH=scripts python3 scripts/wp13_m0.py identity_anchors",
    "PYTHONPATH=scripts python3 scripts/wp13_m0.py identity_masses",
    "PYTHONPATH=scripts python3 scripts/wp13_m0.py identity_inject",
    "PYTHONPATH=scripts python3 scripts/wp13_m0.py identity_fit",
    "PYTHONPATH=scripts python3 scripts/wp13_score.py I1",
    "PYTHONPATH=scripts python3 scripts/wp13_m0.py ages        # W13-I4",
    "PYTHONPATH=scripts python3 scripts/wp13_m0.py anchors",
    "PYTHONPATH=scripts python3 scripts/wp13_m0.py masses",
    "PYTHONPATH=scripts python3 scripts/wp13_m0.py inject",
    "PYTHONPATH=scripts python3 scripts/wp13_m0.py fit",
    "PYTHONPATH=scripts python3 scripts/wp13_score.py t1",
    "PYTHONPATH=scripts python3 scripts/wp13_score.py ledger    # W13-I2",
    "PYTHONPATH=scripts python3 scripts/wp13_t7.py",
    "PYTHONPATH=scripts python3 scripts/wp13_score.py score",
]

CLARIFICATIONS = {
    "A1_T1_same_data": (
        "A Bayes factor needs both models to see the same data.  T1 is therefore computed on the A and "
        "C test-c stars only (59 + 43 = 102), on both sides: ln BF = ln Z(A) + ln Z(C) - ln Z(A+C), "
        "where Z(A+C) uses the summed curve.  B's 4 stars enter neither side of T1 (they have no M1 "
        "age of their own, D2).  M0's age for the ledger still uses all 106 stars (D1(a)).  The "
        "test-c ln L is a membership-weighted sum over stars (scripts/issue20b_run.py), so summing the "
        "stored per-subgroup curves is exact; the likelihood is membership-weighted, so the BF is a "
        "weighted-likelihood BF and is reported as such.  Prior: uniform in log10 age over the native "
        "grid (WP4's convention), wp13_rules.log_evidence."),
    "A2_T2_T3_intervals": (
        "T2's paired 95 % interval is taken on the per-iteration EXPECTED count mu = k int M^-alpha "
        "above the turnoff (parameter uncertainty), paired by iteration index -- not on the Poisson "
        "realisation, whose 95 % interval (about +-7 at N ~ 7) could never exclude a +-10 % band, "
        "which would make T2 unpassable by construction.  T3's interval is the Monte Carlo interval "
        "of the paired indicator difference.  N_death and P(last < 100 kyr) are reported as in WP7 "
        "(realised counts), unchanged."),
    "A3_T7_design": (
        "Ledger-level injections (brief §5), 200 realisations per design, 20,000 ledger iterations per "
        "model per realisation, PARSEC R_V 3.1 alpha 2.3 delta 0, all-explode.  Truth k_s = the M1 "
        "k_median of each subgroup.  Fitted ages: t_hat ~ N(t_true, sigma) and posterior draws "
        "~ N(t_hat, sigma); fitted k: k_true x lognormal with the M1 relative k width; sigma_s = the "
        "M1 WP5 truth-age posterior sd of subgroup s, sigma_0 = M0's measured pooled sd.  "
        "T7a (coeval): every subgroup at t = 3.516 Myr (the disclosed k-weighted mean), M0 fitted "
        "at the same truth; T2 or T3 'passes' (wp13_rules.t2/t3) in < 10 % of realisations; the "
        "same realisations give T4's D null (95th percentile).  T7b (multi-age): true ages = the M1 "
        "WP5 posterior means (A 3.48, B 3.70, C 3.38 Myr); M0's fitted truth = the k-weighted mean "
        "true age; pass iff |bias_N(M1)| < |bias_N(M0)| and |bias_P(M1)| < |bias_P(M0)| and M1's "
        "68 % interval of mu covers mu_true in 58-78 % of realisations.  Disclosed weakness: at "
        "ledger level the M1 coverage is near-tautological (its posterior is centred on a draw with "
        "the true width); it tests propagation, not the age fit."),
    "A4_D6_grid": (
        f"(age gap x young fraction): a fraction f of the association's k (sum of M1 k_median) sits "
        f"at 3.516 - gap Myr, the rest at 3.516 Myr; gap in {list(R.D6_AGE_GAP_MYR)} Myr, f in "
        f"{list(R.D6_YOUNG_FRACTION)}; 200 realisations each; M1 = the two components fitted "
        f"separately, M0 = one age at the k-weighted mean.  The one-age estimate is called adequate "
        f"where |bias N| <= {R.D6_ADEQUATE_N_REL:.0%} of the truth and |bias P| <= "
        f"{R.D6_ADEQUATE_P_ABS}.  Descriptive: the map frames the paper; it gates nothing."),
    "A5_M0_implementation": (
        "M0 is the repair_v9 chain re-run with ONE label, POOLED_LABEL, substituted in memory for "
        "A, B and C (the 61 unlabelled members stay excluded).  No engine file is edited: "
        "scripts/wp13_m0.py runs each engine in-process (runpy) after two in-memory patches -- the "
        "labels table is read back pooled and the shared SUBGROUPS list becomes [POOLED_LABEL]; the "
        "driver's versions join wp5_joint_age_fit.AGE_INTERPOLATED_VERSIONS (repair_v9's node rule).  "
        "The engines it runs are byte-identical to repair_v9's frozen_code (checked when this record "
        "was written; chain.py differs only by the adoption flag).  With pooling off the same driver "
        "must reproduce repair_v9 (W13-I1).  Steps: (1) pooled age "
        "table: the summed A + B + C test-c curves (eps 0.05) per family x R_V, summarised by "
        "wp4v9_common.posterior, one 'ums' row per key, flagged 'pooled_A_B_C'; (2) anchors and "
        "mass posteriors at the pooled age (wp4_anchors_hrd, wp4_mass_posteriors_repair); (3) WP5 "
        "truth-age node injections for the pooled footprint (donors = union of the three subgroups' "
        "donors), 9 nodes x 6 family/R_V cells, via repair_v9_injections.py; (4) WP5 joint fit "
        "(wp5_fit_imf_joint.py) -> pooled (k, truth age) draws per family x R_V x alpha; (5) ledger: "
        "wp7_ledger.run_population on the 54 branches, 2,000,000 iterations, own seeds; M1 replayed "
        "per subgroup with the stored seeding, paired by iteration index.  WP6 is not re-run (M0 "
        "needs no closure).  M0-lite (k_ALL = k_A + k_B + k_C) is computed as a labelled fallback and "
        "may not be quoted as the best one-population model."),
    "A6_verdict_gap": (
        "The brief's verdict rule leaves two combinations unnamed: T1 fails but an effect test passes, "
        "and T1 and an effect pass but T6 fails.  They are reported as 'unclassified' with the "
        "individual results, never mapped to a named outcome afterwards (wp13_rules.verdict)."),
}

INTEGRITY = {
    "W13-I1": ("driver inert: with pooling off, scripts/wp13_m0.py reproduces repair_v9's anchors, mass "
               "posteriors, the 27 PARSEC R_V 3.1 WP5 nodes, and (identity masses + repair_v9 responses) "
               "the WP5 k_median in all 54 cells, within repair_v9_integrity.TOL (I1a-I1d tolerances); "
               "the ledger part is W13-I2"),
    "W13-I2": ("M1 replay: re-running the WP7 engine with the stored seeding reproduces "
               "tables/wp7_ledger_repair_v9.csv N_SN_mean and P_last_SN_within_100kyr on every scored "
               "and reported branch within I1e_ndeath_abs and 0.005"),
    "W13-I3": "every consumed input and frozen-code hash matches this record (no exception)",
    "W13-I4": ("pooled-curve construction: the A + C sum equals the stored B-borrowed curve in "
               "provenance/repair_v9_headline_ages_execution.json (b_combined_curves) to 1e-9; the "
               "pooled table has n_stars = 106 at every family x R_V"),
    "W13-I5": ("held fixed: same members and membership probabilities, 1,331 labelled stars pooled, "
               "61 unlabelled excluded from both, extinction repair_v5, the same isochrone tables, IMF "
               "slopes and the 120 Msun ceiling, the 2-8 Msun calibration window, f_bin, the formation "
               "durations, all-explode primary, turnoff relation and RNG recipe -- asserted in the "
               "execution record"),
}


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"{OUT.relative_to(w.ROOT)} exists; a pre-registration is written once")
    if ADDENDUM_SIGNOFF is None:
        raise SystemExit("the owner's sign-off of clarifications A1-A6 is not recorded yet")
    missing = [p for p in CODE if not (w.ROOT / p).exists()]
    if missing:
        raise SystemExit(f"WP13 code not written yet, cannot be hashed: {missing}")
    present = [p for p in OUTPUTS_THAT_MUST_NOT_EXIST if (w.ROOT / p).exists()]
    if present:
        raise SystemExit(f"WP13 outputs already exist, cannot pre-register blind: {present}")

    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp13_prereg.py",
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "chain": "repair_v9 (M1, adopted 2026-10-07)",
        "basis": {"brief": INPUTS["brief"], "proposal": INPUTS["proposal"],
                  "decisions": INPUTS["decisions_2026_10_08"]},
        "owner_signoff": {
            "proposal": ("2026-10-08: issue #22 accepted as a caveat; D1(a); D2; D3 = T6 scored on the 18 "
                         "alpha = 2.3 branches, the 36 reported; D4; D5; D6 -- "
                         "provenance/decisions_2026_10_08.json"),
            "clarifications_A1_A6": ADDENDUM_SIGNOFF,
            "logged_in": "AI_PROMPTS.md",
        },
        "falsifiable_sentence": (
            "Relative to the best single-age, single-normalization model of the same Gaia-selected stars "
            "(M0), the subgroup-resolved model (M1) changes the inferred number of stellar deaths by at "
            "least T_N and the probability of an event within the last 100 kyr by at least T_P, because "
            "Cyg OB2-C lies below the first-death boundary while Cyg OB2-A and -B lie above it.  (Kept as "
            "written; its 'because' clause is already false on repair_v9.)"),
        "models": {
            "M1": "the adopted repair_v9 ledger, unchanged, hash-verified, replayed for pairing (W13-I2)",
            "M0": "one age, one k per family x R_V x alpha cell, on the 1,331 labelled members (A5)",
            "M0_lite": "M0's pooled age with k_ALL = k_A + k_B + k_C draw-wise; fallback, never 'the best'",
        },
        "branches": {"scored_T6": "18: PARSEC/MIST x R_V 3.0/3.1/3.5 x delta 0/1/2 Myr at alpha 2.3",
                     "reported": "36: the same at alpha 2.0 and 2.3", "all_54": "reported in the table"},
        "tests": {
            "T1": f"ln BF(M1:M0) >= {R.T1_LNBF_MIN} on the baseline cell and on >= {R.T1_MIN_CELLS} of 6 "
                  f"family x R_V cells (A1)",
            "T2": f"|dN| >= {R.T2_ABS} and |dN|/N_M1 >= {R.T2_REL}, paired 95 % interval excluding "
                  f"+-{R.T2_BAND_REL:.0%} of N_M1 (A2)",
            "T3": f"|dP(last < 100 kyr)| >= {R.T3_ABS}, paired 95 % interval excluding 0 (A2)",
            "T4": f"D >= {R.T4_D_MIN} and D > the {R.T4_NULL_QUANTILE:g}th percentile of the T7a null, or "
                  f"the first-death epoch shifts by >= {R.T4_FIRST_DEATH_SHIFT_MYR} Myr",
            "T5": f"for some subgroup |N_M1,sub - N_M0 k_sub / sum k| >= {R.T5_ABS} (baseline)",
            "T6": f"sign of dN and dP, and whichever of T2/T3 passed on the baseline, hold on >= "
                  f"{R.T6_FRACTION:.0%} of the 18 scored branches and on both families",
            "T7a": f"coeval injections: T2 or T3 passes in < {R.T7A_MAX_FALSE_PASS:.0%} of "
                   f"{R.T7_REALISATIONS} realisations (A3)",
            "T7b": f"multi-age injections: M1 bias smaller than M0's in N and in P; M1 68 % coverage in "
                   f"{list(R.T7B_COVERAGE)} (A3)",
            "baseline_cell": R.BASELINE,
            "explodability": "all-explode primary; hard cutoff reported as a secondary row",
        },
        "verdict_rule": {
            "material improvement": "T1 and (T2 or T3 or T4) and T6 and T7a and T7b",
            "supported but immaterial": "T1 and T7 and not (T2 or T3 or T4)",
            "equivalent": "not T1 and not (T2 or T3 or T4)",
            "overfit": "not T7a or not T7b, whatever else passes",
            "unclassified": "any other combination (A6)",
            "code": "wp13_rules.verdict",
        },
        "framing_by_outcome": "brief §6, unchanged (D6); 'equivalent' adds the D6 grid (A4)",
        "clarifications": CLARIFICATIONS,
        "run_order": RUN_ORDER,
        "plan_tests_before_this_record": (
            "wp13_m0.py identity_plan / plan (dry runs: node plans only, nothing written; the pooled plan "
            "used the headline table's subgroup median, not an M0 age) and a smoke test of wp13_score "
            "ledger and wp13_t7 on a synthetic M0 stand-in (subgroup A's repair_v9 draws relabelled), "
            "2,000 iterations, 3 realisations, outputs in the session scratchpad -- no M0 number existed"),
        "integrity_checks": INTEGRITY,
        "thresholds": {k: getattr(R, k) for k in dir(R) if k.isupper()},
        "disclosed_prior_knowledge": [
            "brief §2 as re-based on repair_v9 (tables/wp13_brief_rebase_repair_v9.csv): a single "
            "k-weighted common age of 3.516 Myr gives N_death 6.644 against the resolved 6.657 and "
            "P(last < 100 kyr) 0.697 against 0.697; T2 and T3 are expected to fail and the expected "
            "outcome is 'equivalent' (or 'supported but immaterial' if T1 passes)",
            "M0's age is close to known: the A + C test-c product is B's borrowed headline posterior "
            "(PARSEC 3.55, MIST 3.57 Myr at R_V 3.1; provenance/repair_v9_headline_ages_execution.json), "
            "and M0's 106-star curve differs from it only by B's 4 stars",
            "at equal ages the three subgroups' death curves agree within 7 % (repair_v9 age scan)",
            "T1 was not computed; the stored curves were opened only to list their keys "
            "(c|0.05|CygOB2-B exists), no value was read",
            "issue #22 (WP5 counts vs spectroscopic ages) is open and is a caveat shared by M0 and M1 (D4)",
        ],
        "anti_tuning": (
            "brief §9: no WP13 script alters any upstream artifact, threshold, gate or branch set; every "
            "input is hash-verified; M1's numbers are the stored ones, replayed for pairing, never moved; "
            "an unwelcome M0 is reported; failed criteria stay failed; a defect in a check is fixed and "
            "recorded in provenance/wp13_deviations.json, a threshold is never loosened"),
        "consumed_inputs": {name: {"path": rel, "sha256": w.sha256(w.ROOT / rel)}
                            for name, rel in INPUTS.items()},
        "frozen_code": {rel: w.sha256(w.ROOT / rel) for rel in CODE},
    }
    w.write_json(OUT, record)
    print(f"wrote {OUT.relative_to(w.ROOT)}  sha256 {w.sha256(OUT)}")


if __name__ == "__main__":
    main()
