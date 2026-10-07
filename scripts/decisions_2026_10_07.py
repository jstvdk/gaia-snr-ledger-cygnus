#!/usr/bin/env python3
"""Record the project owner's decisions of 2026-10-07 on the repair_v9 Stage 1
outcome (reports/wp4v9_age_redesign.md §8).

Written once; refuses to overwrite.  The decisions are stored with the evidence
they rest on (paths and SHA-256), what they change, and what they do NOT change.
Nothing is recomputed here.

Run:
  PYTHONPATH=scripts python3 scripts/decisions_2026_10_07.py
"""
from __future__ import annotations

from datetime import datetime, timezone

import wp5_common as w

OUT = w.PROVENANCE / "decisions_2026_10_07.json"

EVIDENCE = [
    "reports/wp4v9_age_redesign.md",
    "provenance/wp4v9_age_prereg.json",
    "provenance/wp4v9_age_outcome.json",
    "tables/wp4v9_injection_validation.csv",
    "tables/wp4v9_real_fits.csv",
    "tables/wp4v9_posthoc_refits.csv",
    "tables/wp4v9_posthoc_star_contributions.csv",
    "tables/issue20b_age_tests.csv",
]

DECISIONS = {
    "6_accept_stage1_outcome": {
        "question": "Accept the frozen Stage 1 outcome for A and C: GA1 failed in every cell, so "
                    "repair_v9 carries the Phase A' test-c spectroscopic-HRD posteriors (eps 0.05) "
                    "for A (PARSEC 3.16, MIST 3.57 Myr) and C (3.55, 4.01 Myr)?",
        "decision": "Accepted, with full documentation: the ages are spectroscopic-HRD ages (Gaia G, "
                    "Gaia-based distance and WP3 extinction on the luminosity axis; literature "
                    "spectral-type temperatures on the other), not Gaia-photometry-only ages, and the "
                    "failure of the Gaia optical CMD age (injections, real-data instability, one-star "
                    "dominance) is reported as a result in its own right.",
        "owner_concern_recorded": "whether adopting spectroscopic ages diminishes the Gaia work; the "
                                  "answer given: membership, subgroups, distance, extinction, the "
                                  "census that sets k, and the kinematics stay Gaia products; only "
                                  "the claim of a Gaia-photometry-only age is lost.",
    },
    "7_age_scan_all_subgroups": {
        "question": "B has no usable age (flagged M1 1.4-5.7 Myr across R_V; 4 spectroscopic stars). "
                    "Options offered: (a) borrow A+C's spectroscopic age, (b) fixed-age branches, "
                    "(c) keep the flagged M1 row.",
        "decision": "Every subgroup (A, B and C) is scanned over all plausible ages: the 11 native "
                    "isochrone ages from 2.0 to 6.3 Myr (PARSEC 2.00-6.31, MIST 2.00-6.37), each "
                    "subgroup independently, with the WP5 normalisation k refitted at every age. "
                    "Subgroup ledgers are independent Poisson populations summed per iteration, so any "
                    "(t_A, t_B, t_C) combination follows from the three one-dimensional scans. The "
                    "full scan is published.",
        "headline_weighting": "A and C: their own test-c spectroscopic posteriors. B: the A+C combined "
                              "spectroscopic posterior (product of the two test-c likelihood curves, "
                              "same family and R_V), flagged 'not measured for B' wherever quoted.",
        "alternatives_declined": ["full 1-10 Myr grid", "scan only, no headline",
                                  "B uniform over the scan"],
    },
}

CHANGES = ("the repair_v9 pre-registration (brief §4.1) adds the per-subgroup age-scan dimension and "
           "the headline weighting above; it supersedes the frozen flagged-M1 row as B's headline "
           "age (that row stays in wp4_age_posteriors_repair_v9.parquet, unedited)")
DOES_NOT_CHANGE = ("the Stage 1 scoring (GA1-GA3, W1-W4) and every Stage 1 output; repair_v8 stays the "
                   "quoted chain until repair_v9 passes its pre-registered adoption rule")


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"{OUT.relative_to(w.ROOT)} exists; a decision record is written once")
    w.write_json(OUT, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/decisions_2026_10_07.py",
        "decided_by": "project owner, 2026-10-07 ('I accept it, but with proper documentation'; "
                      "'All subgroups should test all possible ages'; then chose the 2.0-6.3 Myr "
                      "native grid and the A,C-spectroscopic / B-borrows headline)",
        "decisions": DECISIONS,
        "changes": CHANGES,
        "does_not_change": DOES_NOT_CHANGE,
        "evidence": {p: w.sha256(w.ROOT / p) for p in EVIDENCE},
    })
    print(f"wrote {OUT.relative_to(w.ROOT)}")


if __name__ == "__main__":
    main()
