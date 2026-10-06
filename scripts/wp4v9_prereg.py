#!/usr/bin/env python3
"""repair_v9 Stage 1 -- pre-registration of the WP4 age redesign (issue #21).

Brief: tasks/stage_repair_v9_wp4_age_redesign_and_wp13.md, section 3.
Owner decision: provenance/decisions_2026_10_06.json (option b).

Written before any real-data posterior or injection of the redesign exists.
Refuses to overwrite itself and refuses to run if any Stage 1 output exists.
`--dry-run` writes the record to a scratch path for owner review instead.

Run:
  PYTHONPATH=scripts python3 scripts/wp4v9_prereg.py [--dry-run PATH]
"""
from __future__ import annotations

import argparse
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

import wp4v9_common as V
import wp5_common as w

OUTPUTS = [
    "data/processed/wp4_age_posteriors_repair_v9.parquet",
    "tables/wp4v9_injection_validation.csv",
    "tables/wp4v9_faint_end_study.csv",
    "tables/wp4v9_real_fits.csv",
    "provenance/wp4v9_injection_execution.json",
    "provenance/wp4v9_fit_execution.json",
    "provenance/wp4v9_age_outcome.json",
    "reports/wp4v9_age_redesign.md",
]

METHOD = {
    "M1_extinction": (
        "wp4v9_common.m1_loglike.  Each star's likelihood is integrated over its own WP3 A_V "
        "posterior (wp3_extinction_posterior_repair_v5.npz, at the branch's R_V), represented by "
        "K = 32 equal-probability quantile nodes.  For each node the star is de-reddened along "
        "r = (k_BP - k_RP, k_G) from band_coefficients(R_V) and compared with WP4's particle "
        "mixture (wp4_common.build_model_particles, unchanged: 400 log-mass points, 24 q values, "
        "IMF 2.3, f_bin of the branch) under photometric errors + MAG_FLOOR 0.02 + SIGMA_INT 0.03 "
        "only.  Anchors carry their narrow intrinsic-colour posterior automatically.  Stars with "
        "no finite posterior (9 of 1,392) are excluded."
    ),
    "M2_window": (
        "Primary window M_G0 <= -1.0; sensitivities -0.5 and -1.5.  Membership is decided on the "
        "star's posterior-MEDIAN M_G0.  CHOSEN RULE: the window is modelled as a SELECTION "
        "FUNCTION inside the likelihood.  For star i the cut is a cut on its observed g = G - mu at "
        "the fixed edge + k_G med_i(A_V); the likelihood is the untruncated population density "
        "(particles to edge + 3 mag) marginalised over the A_V nodes, divided by the probability "
        "that a model star with the same A_V nodes and photometric error passes the same cut "
        "(Gaussian CDF).  The alternative the brief offered -- the same median cut applied to "
        "the model particles (truncation) -- was tested on synthetic data before registration "
        "and recovered 1.5 Myr for a true 2.51 and 1.1 Myr for a true 3.98 (stars scatter into "
        "the window from the much more numerous fainter population because the A_V posteriors "
        "are 0.75 mag wide), so it is not used."
    ),
    "M3_unchanged": (
        "forward model and particles of wp4_common; f_bin in {0.3, 0.4, 0.5}; native grid 1-10 Myr; "
        "PARSEC and MIST; R_V in {3.0, 3.1, 3.5}; membership-probability weights; distance fixed at "
        "1.6245 kpc with +/- 0.060 mag refits at the decision cell (model magnitudes shifted, as "
        "WP4); measurability gate (>= 15 stars, not grid-railed) and posterior convention "
        "(wp4_common.posterior_from_loglike) with 95 % bounds added."
    ),
    "secondary_error_model": (
        "the Gaussian correlated-error model (wp4v9_common.gauss_loglike, correlated=True; A_V at "
        "the posterior median, sigma the posterior 68 % half-width, window truncation as WP4) "
        "on every cell -- sensitivity only.  The old WP4 model (independent errors, WP3 A_V and "
        "av_err, truncation) is refitted on every cell for comparison."
    ),
    "double_counting_sensitivity": (
        "Added at the owner's request on 2026-10-06, before registration.  WP3's A_V posterior was "
        "fitted to the same G, BP, RP (with an age-agnostic template prior: every native age of "
        "both families weighted equally), so M1 uses the optical colours twice and its intervals "
        "may be too narrow.  Synthetic tests cannot detect this by construction.  Check: refit M1 "
        "in the primary window with every star's A_V nodes widened x2 about the median "
        "(wp4v9_common.widened; clipped at 0; window membership unchanged), every subgroup x "
        "family x R_V at f_bin 0.4.  Reported: the shift of the posterior median and the change "
        "of the 68 % width.  Report only; it gates nothing.  A shift larger than the primary 68 % "
        "half-width is flagged 'double-counting sensitive' wherever the age is quoted."
    ),
    "M4_faint_end_study": {
        "base": "faint window -1 < M_G0 <= 1.5 alone, M1, decision cell, both families",
        "causes_one_at_a_time": {
            "i_pms": "MIST only: model particles without the PMS phase (phase -1); PARSEC's "
                     "label column is not a reliable phase flag in these tables (MS stars of "
                     "14-53 Msun carry label 0) and is not used",
            "ii_binaries": "f_bin 0.6 and 0.7",
            "iii_contamination": "membership probability P > 0.9 only",
            "iv_av_outliers": "outlier mixture eps = 0.05 with a uniform density 0.2 mag^-2 "
                              "(colour range 2.0 x magnitude band 2.5)",
        },
        "reconciles_if": "the faint-window 68 % interval overlaps the bright-window (primary) 68 % interval",
        "gates": "nothing (diagnostic, never adopted)",
    },
}

INJECTIONS = {
    "design": (
        "wp4v9_common.synthetic per subgroup x family at true ages 2.5, 3.2, 4.0 and 5.0 Myr "
        "(snapped to native: PARSEC 2.51 3.16 3.98 5.01; MIST 2.52 3.18 4.01 5.05), R_V 3.1, "
        "f_bin 0.4.  Each synthetic star takes a random real star's photometric errors and A_V "
        "posterior (the subgroup's labelled members); its TRUE A_V is a draw from that posterior, "
        "its intrinsic photometry a draw from the model particles; the window is applied on the "
        "posterior-median M_G0 exactly as for real data; the star count equals the real count in "
        "the window (subgroup, R_V 3.1)."
    ),
    "realisations": "100 per cell for the primary (M1, bright window); 20 per cell for M1 and the "
                    "old WP4 model in the faint and full windows and for the old model in the "
                    "bright window (report only, and W3)",
    "point_estimate": "posterior median",
    "seed": 20261006,
    "GA1": (
        "PASS for a subgroup x family iff, for M1 in the bright window, |mean(median - truth)| < "
        "0.2 Myr at every one of the four true ages AND the 68 % coverage pooled over the 400 "
        "realisations lies in [0.55, 0.80]"
    ),
}

ACCEPTANCE = {
    "evaluated_on": "each subgroup x family at the decision cell (R_V 3.1, f_bin 0.4); robustness "
                    "over the other R_V x f_bin cells is a qualifier (robust / fragile)",
    "GA1": "injection validation as above",
    "GA2": (
        "A and C only: the M1 bright-window 68 % interval overlaps the Phase A' test-c spectroscopic "
        "68 % interval (tables/issue20b_age_tests.csv, test c, eps 0.05, same family, R_V 3.1)"
    ),
    "GA3": (
        "less than 5 % of the age posterior lies at ages whose turnoff "
        "(wp6_mass_extension_decision.turnoff_mass) is below m_lo, the largest lower bound on "
        "initial mass among the subgroup's spectroscopic members of class III-V (or "
        "unclassified), non-composite, O/B primary, not extreme-hot.  Each star's bound is "
        "wp4v9_common.minimum_initial_mass: the smallest initial mass of any tabulated point of "
        "the family's isochrones (any native age, any phase) within chi <= 2 of the star "
        "(sigma 0.03 dex, 0.40 mag; table T_eff; M_G0 of the branch).  A bound that assumes no age."
    ),
    "B": "no spectroscopic check (4 stars): accepted on GA1 and GA3 only and flagged "
         "'photometry-only' wherever quoted",
    "adopted_if_pass": "the M1 bright-window posterior, for every R_V x f_bin cell of that subgroup x family",
    "fallback_if_GA2_fails": (
        "the JOINT likelihood: M1 bright window over the subgroup's NON-anchor stars (disjoint "
        "data) x the test-c spectroscopic HRD curve of Phase A' (eps 0.05, same family and R_V).  "
        "Used for every R_V cell of that subgroup x family (test c has no f_bin variants beyond "
        "0.4, so the joint uses the f_bin 0.4 spectroscopic curve throughout)."
    ),
    "branches_if": (
        "the joint posterior is bimodal (wp4v9_common.bimodal: two maxima, separating minimum < "
        "0.5 x the smaller peak, smaller mode >= 10 % of the probability) or fails GA3: carry the "
        "M1 bright-window and the test-c posteriors as two separate branches through repair_v9; "
        "no claim may then rest on the difference between subgroup ages"
    ),
    "if_GA1_or_GA3_fails": (
        "the subgroup x family is NOT ADOPTED from the redesign; repair_v9 then carries the test-c "
        "spectroscopic posterior (A, C) as the only age, or, for B, the M1 bright-window posterior "
        "flagged 'failed GA1/GA3'.  Stated in the outcome; the owner decides before repair_v9"
    ),
}

PREDICTIONS = {
    "W1": {"statement": "C's adopted age is >= 3.2 Myr in both families",
           "scored_on": "adopted posterior MAP, decision cell"},
    "W2": {"statement": "|A - C| adopted age difference < 1 Myr in both families",
           "scored_on": "adopted posterior MAPs, decision cell"},
    "W3": {"statement": "the old independent-error model is biased young in the faint window on injections",
           "scored_on": "mean(median - truth) of the old model, faint window, averaged over the 12 "
                        "subgroup x true-age cells, is < -0.2 Myr in both families"},
    "W4": {"statement": "the faint-window study does not reconcile the faint window under any single cause",
           "scored_on": "a cause reconciles if it brings the faint 68 % into overlap with the bright "
                        "68 % for every subgroup with a measurable fit in every family it applies to; "
                        "W4 PASS iff no cause reconciles (if the M1 faint window already overlaps "
                        "without any cause, W4 is recorded as moot)"},
}

OUTPUT_SCHEMA = (
    "wp4_age_posteriors_repair_v9.parquet: the repair_v5 columns (subgroup, family, R_V, f_bin, "
    "indicator, dmu, n_stars, measurable, grid_railed, exclusion_reason, age_map, age_lo68, "
    "age_hi68, age_lo90, age_hi90, age_mean) plus window, error_model, adopted, fallback, "
    "branch, photometry_only.  Exactly one row per (subgroup, family, R_V, f_bin, dmu) carries "
    "indicator == 'ums' and adopted == True -- the row downstream code selects; every other fit "
    "is stored with indicator == 'ums_sensitivity'.  If branches are triggered, the second "
    "branch is stored with indicator 'ums_branch_spectroscopic' and repair_v9's preregistration "
    "must add the branch dimension before WP5."
)

DISCLOSED = [
    "Issue #21 diagnostic (tables/issue21_error_model_diagnostic.csv, PARSEC / MIST, R_V 3.1, "
    "f_bin 0.4): WP4 full window A 3.98 / 4.01, B 3.55 / 4.01, C 2.51 / 3.18; correlated full "
    "window A 5.62 / 6.37, B 5.62 / 8.03, C 5.01 / 4.50; bright window (M_G0 <= -1, either model) "
    "A 2.82 / 4.01 (broad, 68 % 1.1-2.8 PARSEC), B 3.98 / 2.83 (correlated 2.82 / 4.01), C 3.98 / 3.18.",
    "Phase A / A' (provenance/issue20_phase_a_outcome.json, issue20b_outcome.json): spectroscopic "
    "HRD A 3.16 / 3.57, C 3.55 / 4.01 (test c); near-IR free-extinction CMD A 5.0 / 7.2, B 7.9 / "
    "8.0, C 4.5 / 4.5 (full window, considered mis-specified); ESP-HS invalid.",
    "The previous agent expects C at ~3.2-4 Myr.  A's bright window was broad in the diagnostic, "
    "so GA2 for A-PARSEC may fail (brief 3.4).",
    "Before registration the present agent ran the method on SYNTHETIC data only: the median-"
    "truncation window rule failed (1.5 / 1.1 Myr for true 2.51 / 3.98, M2 above); the "
    "selection-function rule recovered 2.55 [2.43, 2.75], 3.90 [2.73, 4.33], 4.65 [3.85, 5.20] for "
    "true 2.51, 3.98, 5.01 (C-like, PARSEC, one realisation each); the old model on those same "
    "synthetic sets gave 2.48 and 2.19 for true 2.51 and 3.98 (median-truncation selection).  "
    "Real stars were used only for their photometric errors, A_V posteriors and the window "
    "counts already in the brief; no real-data M1 likelihood has been computed.",
    "The A_V posteriors of photometric stars have a median 68 % half-width of 0.75 mag (anchors "
    "0.03); 9 stars have no posterior.",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", type=Path, default=None)
    args = ap.parse_args()
    out = args.dry_run or V.PREREG_PATH
    if args.dry_run is None:
        if V.PREREG_PATH.exists():
            raise SystemExit("provenance/wp4v9_age_prereg.json exists -- written once")
        existing = [p for p in OUTPUTS if (w.ROOT / p).exists()]
        if existing:
            raise SystemExit(f"Stage 1 outputs already exist: {existing}")
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp4v9_prereg.py",
        "dry_run": args.dry_run is not None,
        "brief": "tasks/stage_repair_v9_wp4_age_redesign_and_wp13.md section 3",
        "decision_record": "provenance/decisions_2026_10_06.json",
        "python": sys.version.split()[0], "platform": platform.platform(),
        "packages": {"numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
        "method": METHOD,
        "constants": {"K_nodes": V.K_NODES, "windows": {k: list(v) for k, v in V.WINDOWS.items()},
                      "primary_window": V.PRIMARY_WINDOW, "selection_margin_mag": V.SELECTION_MARGIN,
                      "decision_cell": V.DECISION, "widen_factor": V.WIDEN_FACTOR, "f_bins": list(V.F_BINS),
                      "ga3": {"sig_logte": V.SIG_LOGTE_GA3, "sig_mg0": V.SIG_MG0_GA3,
                              "chi": V.CHI_GA3, "max_posterior_fraction": 0.05}},
        "injections": INJECTIONS,
        "acceptance": ACCEPTANCE,
        "predictions": PREDICTIONS,
        "output_schema": OUTPUT_SCHEMA,
        "disclosed_prior_knowledge": DISCLOSED,
        "method_code": {rel: w.sha256(w.ROOT / rel) for rel in (
            "scripts/wp4v9_common.py", "scripts/wp4_common.py", "scripts/wp4_fit_ages.py",
            "scripts/issue20_common.py", "scripts/wp6_mass_extension_decision.py")},
        "consumed_inputs": {n: {"path": rel, "sha256": w.sha256(w.ROOT / rel)}
                            for n, rel in V.INPUTS.items()},
        "planned_outputs": OUTPUTS,
        "failed_predictions_rule": "results stay as they come out; this file is never amended",
    }
    w.write_json(out, record)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
