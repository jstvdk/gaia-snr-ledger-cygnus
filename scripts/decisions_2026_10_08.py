#!/usr/bin/env python3
"""Record the project owner's decisions of 2026-10-08: issue #22 and the WP13
pre-registration proposal (tasks/wp13_prereg_proposal.md, decisions D1-D6).

Written once; refuses to overwrite.  The decisions are stored with the evidence
they rest on (paths and SHA-256), what they change, and what they do NOT change.
Nothing is recomputed here, and no M0-side number is read.

Run:
  PYTHONPATH=scripts python3 scripts/decisions_2026_10_08.py
"""
from __future__ import annotations

from datetime import datetime, timezone

import wp5_common as w

OUT = w.PROVENANCE / "decisions_2026_10_08.json"

EVIDENCE = [
    "reports/repair_v9_completion_report.md",
    "tasks/wp13_prereg_proposal.md",
    "tasks/wp13_pooled_vs_resolved_ablation_brief.md",
    "provenance/repair_v9_outcome.json",
    "tables/wp13_brief_rebase_repair_v9.csv",
]

DECISIONS = {
    "issue_22": {
        "question": "The WP5 low-mass counts disagree with the spectroscopic ages (gate 28/54 cells, "
                    "C 6/18; C passes only at <= 3.16 / 3.57 Myr).  (a) accept as a stated caveat, "
                    "or (b) investigate first under its own pre-registration?",
        "decision": "(a) accepted as a stated caveat; proceed to WP13.  #22 stays OPEN in "
                    "PROJECT_TRACE §9 and is quoted wherever a repair_v9 number is.",
    },
    "WP13_D1_age_method": {
        "question": "How is M0's pooled age measured?",
        "decision": "(a) the test-c spectroscopic-HRD fit (eps 0.05) on the pooled spectroscopic "
                    "sample A 59 + C 43 + B 4 = 106 stars, same family and R_V -- the method M1 uses "
                    "for A and C.",
        "declined": "(b) the WP4 photometric upper-MS fit on the pooled CMD (Stage 1: it does not "
                    "measure these ages, and it would compare two age methods).",
    },
    "WP13_D2_T1": {
        "question": "T1 as ln BF(M1 : M0) = ln Z_A + ln Z_C - ln Z_pooled on the same grid prior, "
                    "with B excluded on the M1 side (two measured ages); threshold ln BF >= 2.3 on "
                    "the baseline cell and on >= 4 of the 6 family x R_V cells?",
        "decision": "accepted as drafted",
    },
    "WP13_D3_T6_branches": {
        "question": "Which branch set is 'headline' for T6?",
        "decision": "T6 is scored on the 18 alpha = 2.3 branches (stage brief, decision 2 of "
                    "2026-10-02); the 36 retained branches (alpha 2.0 and 2.3) are reported "
                    "alongside, not scored.",
    },
    "WP13_D4_wait_for_22": {
        "question": "Does WP13 wait for issue #22?",
        "decision": "No.  M0 and M1 use identical WP5 counts, so #22 is a caveat shared by both "
                    "sides and is stated as such.",
    },
    "WP13_D5_thresholds": {
        "question": "Keep the brief's §5 thresholds T2-T7 unchanged?",
        "decision": "accepted: unchanged (T2's |dN| >= 3 is 45 % of the repair_v9 baseline 6.66; "
                    "changing it after the disclosure would be tuning).",
    },
    "WP13_D6_verdict_and_framing": {
        "question": "Keep the verdict rule and §6's framing; if 'equivalent', the title becomes "
                    "'When is a fixed-age supernova count adequate? Cygnus OB2 as a benchmark', "
                    "with T7 extended over an (age-spread x young-fraction) grid?",
        "decision": "accepted as drafted",
    },
}

CHANGES = ("tasks/wp13_prereg_proposal.md becomes the basis of scripts/wp13_prereg.py, which is "
           "written, run, hashed and committed before any M0-side number is computed (T1 included)")
DOES_NOT_CHANGE = ("the adopted chain (repair_v9), every repair_v9 product and number, issue #22's "
                   "status (open), and the brief's thresholds")


def main() -> None:
    if OUT.exists():
        raise SystemExit(f"{OUT.relative_to(w.ROOT)} exists; a decision record is written once")
    w.write_json(OUT, {
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "script": "scripts/decisions_2026_10_08.py",
        "decided_by": "project owner, 2026-10-08 (~17:00 CEST), answering four multiple-choice "
                      "questions in the session that followed the SSD data transfer; each "
                      "recommended option was chosen.  The prompt is logged in AI_PROMPTS.md.",
        "decisions": DECISIONS,
        "changes": CHANGES,
        "does_not_change": DOES_NOT_CHANGE,
        "evidence": {p: w.sha256(w.ROOT / p) for p in EVIDENCE},
    })
    print(f"wrote {OUT.relative_to(w.ROOT)}")


if __name__ == "__main__":
    main()
