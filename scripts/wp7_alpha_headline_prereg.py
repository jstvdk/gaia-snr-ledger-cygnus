#!/usr/bin/env python3
"""D1 pre-registration: restrict the HEADLINE branch set to alpha in {2.0, 2.3}.

[wp5_alpha_plausibility.md](../reports/wp5_alpha_plausibility.md) measured how
plausible each IMF slope is and explicitly refused to act: "This is a
plausibility measurement for a decision that has NOT been taken.  If it is ever
acted on, that adoption requires its own pre-registration."  This file is that
pre-registration.  It is written and committed BEFORE the adoption script runs.

WHAT IS BEING DECIDED
  Which branches the manuscript's headline N_SN statement ranges over.  Nothing
  else.  This is a reporting convention, not a measurement, and it is the one
  item in the pre-WP10 brief that is a judgement call rather than a defect fix.

WHAT MAKES IT LEGITIMATE, AND WHAT WOULD MAKE IT TUNING
  The criterion is E1 -- the WP5 calibration-window Poisson chi-square -- which
  is a gate statistic fixed long before WP7 existed and is internal to WP5.  It
  costs nothing to use.  The criterion is explicitly NOT E2, the >8 Msun census
  closure, which is this project's only out-of-sample validation of the IMF
  extrapolation; spending E2 to select alpha would make WP9's answer to
  devil's-advocate objection 1 circular.  That standing decision is
  PROJECT_TRACE section 10 item 2 and it is not being revisited.

  It would be tuning if: the criterion were chosen after seeing which branch set
  gave a preferred N_SN; or the excluded branch were deleted rather than
  reported.  Neither is the case, and the second is checked mechanically below.

HONESTY CLAUSE -- WHAT WAS ALREADY KNOWN WHEN THIS WAS WRITTEN
  Unlike the repair pre-registrations, this one is NOT blind.  The candidate
  branch sets were costed in wp5_alpha_plausibility_execution.json on
  2026-07-29, and while assessing the decision the executing agent additionally
  recomputed the per-set N_SN, P(>=1 SN) and P(last SN < 100 kyr) summaries from
  the stored ledger.  Every such quantity is listed verbatim under
  `known_at_declaration` so that it can never be mistaken for a blind
  prediction.  The predictions scored below are confined to consequences that
  had NOT been computed when this file was written.

Outputs:
  provenance/wp7_alpha_headline_adoption_prereg.json

Run:
  PYTHONPATH=scripts python3 scripts/wp7_alpha_headline_prereg.py
"""
from __future__ import annotations

import platform
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd

import wp5_common as w

LEDGER = w.TABLES / "wp7_ledger.csv"
PLAUSIBILITY = w.PROVENANCE / "wp5_alpha_plausibility_execution.json"


def main() -> None:
    record = {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/wp7_alpha_headline_prereg.py",
        "item": "D1 of tasks/pre_wp10_assessment_brief.md",
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "status": "PRE-REGISTERED, NOT YET EXECUTED",
        "decision": {
            "statement": (
                "The manuscript's HEADLINE branch set for N_SN, the explosion "
                "timeline and the last-supernova posterior is restricted to "
                "alpha in {2.0, 2.3} x {PARSEC, MIST} x R_V {3.0, 3.1, 3.5} x "
                "SF duration {0, 1, 2} Myr, all-explode -- 36 branches.  "
                "alpha = 2.6 is reported in the sensitivity table, never "
                "deleted."
            ),
            "evidence_basis": "E1 only",
            "E1": (
                "WP5 calibration-window Poisson chi-square over 2-8 Msun, an "
                "existing gate statistic.  alpha = 2.6 is the best fit in 1 of "
                "18 cells and has the worst median chi-square (10.80 against "
                "6.86 for 2.3 and 10.30 for 2.0)."
            ),
            "E2_explicitly_not_spent": (
                "the >8 Msun census closure independently agrees (alpha = 2.6 "
                "wins 0 of 18 cells) but is NOT used as the criterion, so it "
                "remains reportable as an out-of-sample validation"
            ),
            "not_adopted": (
                "restricting further to alpha = 2.3 alone.  E1 does not "
                "support it -- alpha = 2.0 wins 7 of 18 cells -- so it would "
                "require E2.  Standing decision, PROJECT_TRACE section 10 item 2."
            ),
        },
        "scope": {
            "changes": [
                "which branches the manuscript quotes as the headline range",
                "the accompanying sentence in PROJECT_TRACE and the WP7 report",
            ],
            "does_not_change": [
                "any stored table, figure or posterior draw",
                "any computed value of N_SN, P(last SN), R_SN(t) or P_verdict",
                "the WP9 framing verdict, whose headline set was ALREADY "
                "alpha in {2.0, 2.3} (36 branches) before this decision",
                "the WP5 branch grid, which keeps all 54 cells",
                "the explodability, R_V, family and SF-duration branch axes",
            ],
        },
        "predictions": [
            {
                "id": "D1-P1",
                "statement": (
                    "The structural finding survives the restriction: on every "
                    "one of the 36 retained branches the ledger returns "
                    "exactly zero supernovae for any black-hole threshold at "
                    "or below 40 Msun."
                ),
                "test": (
                    "tables/wp7_bh_threshold_scan.csv, retained branches only, "
                    "max N_SN over all thresholds <= 40 Msun must be 0.0"
                ),
                "threshold": "exactly 0.0",
                "computed_when_declared": False,
                "if_it_fails": (
                    "the restriction has changed a structural claim and must "
                    "not be adopted without re-deriving that claim"
                ),
            },
            {
                "id": "D1-P2",
                "statement": (
                    "The restriction does not weaken the WP8 pulsar agreement: "
                    "on the retained set P(last SN within the widened 151-401 "
                    "kyr pulsar window) is at least as high as on the full "
                    "54-branch set, because the dropped branches are the "
                    "lowest-rate ones."
                ),
                "test": (
                    "association-scope all-explode P(last SN <= 401 kyr) "
                    "recomputed from tables/wp7_rsn_curves.csv; retained-set "
                    "minimum >= full-set minimum"
                ),
                "threshold": "retained minimum >= full-set minimum",
                "computed_when_declared": False,
                "if_it_fails": (
                    "record as failed; the WP8 agreement would then be "
                    "partly carried by branches the headline no longer quotes"
                ),
            },
            {
                "id": "D1-P3",
                "statement": (
                    "Nothing is deleted.  After adoption, tables/wp7_ledger.csv "
                    "still contains all 54 all-explode association branches "
                    "including all 18 at alpha = 2.6, and its sha256 is "
                    "unchanged from before the adoption ran."
                ),
                "test": "row count and sha256 of tables/wp7_ledger.csv",
                "threshold": "sha256 identical, 18 alpha=2.6 branches present",
                "computed_when_declared": False,
                "if_it_fails": (
                    "the adoption has mutated a stored product and must be "
                    "reverted"
                ),
            },
            {
                "id": "D1-P4",
                "statement": (
                    "The baseline branch is unaffected.  PARSEC, R_V = 3.1, "
                    "alpha = 2.3, coeval, all-explode still gives N_SN = 8.43, "
                    "so the paper's central number does not move as a result "
                    "of a decision about which branches surround it."
                ),
                "test": "baseline row of tables/wp7_ledger.csv",
                "threshold": "N_SN_mean within 0.01 of 8.43",
                "computed_when_declared": True,
                "note": (
                    "declared as a prediction anyway because it is the check "
                    "that separates 'narrowing the reported range' from "
                    "'moving the result'"
                ),
            },
        ],
        "known_at_declaration": {
            "source_1": "provenance/wp5_alpha_plausibility_execution.json (2026-07-29)",
            "source_2": (
                "recomputed from tables/wp7_ledger.csv by the executing agent "
                "while assessing the decision, 2026-07-30, before this file "
                "was written"
            ),
            "association_all_explode_branch_summaries": {
                "all_54": {
                    "N_SN_mean_range": [1.93, 28.74],
                    "N_SN_ensemble_median": 8.79,
                    "factor": 14.93,
                    "P_last_SN_within_100kyr_range": [0.183, 0.889],
                    "P_at_least_one_range": [0.8393, 1.0],
                },
                "retained_36": {
                    "N_SN_mean_range": [5.63, 28.74],
                    "N_SN_ensemble_median": 13.29,
                    "factor": 5.11,
                    "P_last_SN_within_100kyr_range": [0.411, 0.889],
                    "P_at_least_one_range": [0.9935, 1.0],
                },
                "dropped_18_alpha_2p6": {
                    "N_SN_mean_range": [1.93, 4.24],
                    "N_SN_ensemble_median": 3.17,
                    "P_last_SN_within_100kyr_range": [0.183, 0.320],
                },
            },
            "correction_carried_into_the_adoption": (
                "wp5_alpha_plausibility.md section 8 and item D1 of the "
                "pre-WP10 brief both state the retained set as '5.63-28.74, "
                "median approximately 9'.  The ensemble median of the retained "
                "36 branches is 13.29, not 9; 8.79 is the median of the FULL "
                "54-branch set.  The prose picked up the wrong row of "
                "wp5_alpha_plausibility_execution.json, whose "
                "candidate_branch_sets block has the correct 13.29.  The "
                "adoption must not repeat the error, and must state that "
                "dropping alpha = 2.6 RAISES the ensemble median by about 50% "
                "rather than merely trimming a low tail."
            ),
        },
        "reporting_obligations_accepted_with_adoption": [
            "alpha = 2.6 appears in the sensitivity table with its own N_SN "
            "range (1.93-4.24) and its own P_verdict (0.142-0.250)",
            "the exclusion criterion (E1) is named in the text where the "
            "headline range is quoted",
            "the text states that E2 was deliberately not spent, and why",
            "the headline is led by the BASELINE branch value with the carried "
            "range beside it, never by the ensemble median across branches -- "
            "branches are carried, never averaged (plan section 1.4), and an "
            "ensemble median of a non-probability-weighted branch set is not a "
            "posterior summary",
        ],
        "falsification": (
            "This adoption is reversed if any of D1-P1, D1-P2 or D1-P3 fails, "
            "or if a future version of the WP5 fit makes alpha = 2.6 "
            "competitive on E1 (best in more than 6 of 18 cells, or median "
            "chi-square below alpha = 2.3's)."
        ),
        "inputs": {
            str(p.relative_to(w.ROOT)): w.sha256(p)
            for p in (LEDGER, PLAUSIBILITY)
            if p.exists()
        },
    }
    out = w.PROVENANCE / "wp7_alpha_headline_adoption_prereg.json"
    w.write_json(out, record)
    print(f"pre-registered D1 -> {out.relative_to(w.ROOT)}")
    print("  decision:", record["decision"]["statement"][:70], "...")
    print(f"  {len(record['predictions'])} predictions declared, "
          f"{sum(1 for p in record['predictions'] if not p['computed_when_declared'])}"
          f" of them not yet computed")


if __name__ == "__main__":
    main()
